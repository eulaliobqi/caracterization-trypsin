// modules/local/transdecoder.nf — TransDecoder 5.7 para predição de ORFs
// Usa hints de BLAST (UniProt) e Pfam para melhorar predição

process TRANSDECODER {
    tag        "${meta.id}"
    label      'process_high'
    publishDir "${params.outdir}/02_orfs/transdecoder", mode: params.publish_dir_mode
    conda      "${projectDir}/../envs/orf_prediction.yml"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta2), path(blast_hints)   // diamond/blast tsv
    tuple val(meta3), path(hmm_hints)     // hmmer domtblout

    output:
    tuple val(meta), path("${meta.id}.transdecoder.pep"),  emit: pep
    tuple val(meta), path("${meta.id}.transdecoder.cds"),  emit: cds
    tuple val(meta), path("${meta.id}.transdecoder.gff3"), emit: gff3
    path "versions.yml",                                    emit: versions

    script:
    """
    # Passo 1: Extrair ORFs longas candidatas
    TransDecoder.LongOrfs \\
        -t ${fasta} \\
        --min_protein_length ${params.min_aa_length}

    # Passo 2: Predição com hints (BLAST + Pfam melhoram F1)
    TransDecoder.Predict \\
        -t ${fasta} \\
        --retain_blastp_hits ${blast_hints} \\
        --retain_pfam_hits ${hmm_hints} \\
        --single_best_only \\
        --cpu ${task.cpus}

    # Renomear saídas com prefixo meta.id
    mv ${fasta}.transdecoder.pep  ${meta.id}.transdecoder.pep
    mv ${fasta}.transdecoder.cds  ${meta.id}.transdecoder.cds
    mv ${fasta}.transdecoder.gff3 ${meta.id}.transdecoder.gff3

    echo "ORFs preditas: \$(grep -c '>' ${meta.id}.transdecoder.pep)"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        transdecoder: \$(TransDecoder.LongOrfs --version 2>&1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.transdecoder.pep
    touch ${meta.id}.transdecoder.cds
    touch ${meta.id}.transdecoder.gff3
    touch versions.yml
    """
}
