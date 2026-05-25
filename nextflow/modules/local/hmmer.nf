// modules/local/hmmer.nf — HMMER 3.4 para busca de domínios Pfam
// Inclui HMMER (hints gerais) e HMMER_TRYPSIN (Pfam PF00089 Tryp_SPc)

process HMMER {
    tag        "${meta.id}"
    label      'process_medium'
    publishDir "${params.outdir}/02_orfs/hmmer", mode: params.publish_dir_mode
    conda      "${projectDir}/../envs/orf_prediction.yml"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}_pfam.domtblout"), emit: domtblout
    path "versions.yml",                                 emit: versions

    script:
    """
    # Busca contra Pfam completo para hints de TransDecoder
    hmmscan \\
        --domtblout ${meta.id}_pfam.domtblout \\
        --cpu ${task.cpus} \\
        -E ${params.hmmer_evalue} \\
        --domE ${params.hmmer_evalue} \\
        ${params.pfam_db} \\
        ${fasta} \\
        > /dev/null  # suprimir output verbose

    echo "Pfam hits: \$(grep -v '^#' ${meta.id}_pfam.domtblout | wc -l)"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hmmer: \$(hmmscan -h 2>&1 | head -2 | tail -1 | awk '{print \$3}')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_pfam.domtblout
    touch versions.yml
    """
}

process HMMER_TRYPSIN {
    tag        "${meta.id}"
    label      'process_medium'
    publishDir "${params.outdir}/03_trypsin_ids/hmmer", mode: params.publish_dir_mode
    conda      "${projectDir}/../envs/annotation.yml"

    input:
    tuple val(meta), path(pep)

    output:
    tuple val(meta), path("trypsin_hmm.domtblout"),  emit: domtblout
    tuple val(meta), path("trypsin_hmm_ids.txt"),    emit: trypsin_ids
    path "versions.yml",                              emit: versions

    script:
    """
    # Extrair apenas perfil PF00089 (Tryp_SPc) do Pfam
    hmmfetch ${params.pfam_db} Tryp_SPc > trypsin.hmm
    hmmpress trypsin.hmm

    # Busca específica por domínio de tripsina
    hmmsearch \\
        --domtblout trypsin_hmm.domtblout \\
        --cpu ${task.cpus} \\
        -E ${params.hmmer_evalue} \\
        trypsin.hmm \\
        ${pep} \\
        > trypsin_hmm.out

    # Extrair IDs que passaram no critério
    awk '!/^#/ && \$7 <= ${params.hmmer_evalue} {print \$1}' trypsin_hmm.domtblout \
        | sort -u > trypsin_hmm_ids.txt

    echo "HMMER Tryp_SPc hits: \$(wc -l < trypsin_hmm_ids.txt)"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hmmer: \$(hmmsearch -h 2>&1 | head -2 | tail -1 | awk '{print \$3}')
    END_VERSIONS
    """

    stub:
    """
    touch trypsin_hmm.domtblout
    echo "TRINITY_DN1_trypsin_stub" > trypsin_hmm_ids.txt
    touch versions.yml
    """
}
