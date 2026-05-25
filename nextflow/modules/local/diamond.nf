// modules/local/diamond.nf — DIAMOND BLASTp para identificação de tripsinas
// Inclui processo DIAMOND (ORF hints) e DIAMOND_TRYPSIN (identificação específica)

process DIAMOND {
    tag        "${meta.id}"
    label      'process_medium'
    publishDir "${params.outdir}/02_orfs/diamond", mode: params.publish_dir_mode
    conda      "${projectDir}/../envs/orf_prediction.yml"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}_diamond_hints.tsv"), emit: tsv
    path "versions.yml",                                    emit: versions

    script:
    """
    diamond blastp \\
        --query ${fasta} \\
        --db ${params.uniprot_db} \\
        --outfmt 6 qseqid sseqid stitle pident length qlen slen evalue bitscore \\
        --evalue ${params.diamond_evalue} \\
        --max-target-seqs 1 \\
        --sensitive \\
        --threads ${task.cpus} \\
        --out ${meta.id}_diamond_hints.tsv

    echo "DIAMOND hits: \$(wc -l < ${meta.id}_diamond_hints.tsv)"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diamond: \$(diamond --version 2>&1 | awk '{print \$3}')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_diamond_hints.tsv
    touch versions.yml
    """
}

// DIAMOND específico para identificação de tripsinas na Fase 3
process DIAMOND_TRYPSIN {
    tag        "${meta.id}"
    label      'process_medium'
    publishDir "${params.outdir}/03_trypsin_ids", mode: params.publish_dir_mode
    conda      "${projectDir}/../envs/annotation.yml"

    input:
    tuple val(meta), path(pep)

    output:
    tuple val(meta), path("trypsins_confident.fasta"),   emit: confident_fasta
    tuple val(meta), path("trypsins_suggestive.fasta"),  emit: suggestive_fasta
    tuple val(meta), path("identification_report.tsv"),  emit: report_tsv
    tuple val(meta), path("trypsin_ids_confident.txt"),  emit: trypsin_ids
    path "versions.yml",                                  emit: versions

    script:
    """
    # 1. DIAMOND BLASTp contra UniProt
    diamond blastp \\
        --query ${pep} \\
        --db ${params.uniprot_db} \\
        --outfmt 6 qseqid sseqid stitle pident length qlen slen evalue bitscore \\
        --evalue ${params.diamond_evalue} \\
        --max-target-seqs 5 \\
        --sensitive \\
        --threads ${task.cpus} \\
        --out diamond_all.tsv

    # 2. Filtrar hits com "trypsin" no stitle
    awk -F'\\t' 'tolower(\$3) ~ /trypsin/ {print \$1}' diamond_all.tsv | sort -u \
        > diamond_trypsin_ids.txt

    echo "DIAMOND trypsin IDs: \$(wc -l < diamond_trypsin_ids.txt)"

    # 3. Interseção com HMMER (os IDs HMMER devem estar disponíveis via join)
    # Este script também faz a interseção se o arquivo hmm_ids.txt existir
    python3 ${projectDir}/bin/filter_complete_trypsins.py \\
        --diamond diamond_trypsin_ids.txt \\
        --pep ${pep} \\
        --output_confident trypsins_confident.fasta \\
        --output_suggestive trypsins_suggestive.fasta \\
        --report identification_report.tsv \\
        --ids_confident trypsin_ids_confident.txt

    echo "Confident trypsins: \$(grep -c '>' trypsins_confident.fasta)"
    echo "Suggestive trypsins: \$(grep -c '>' trypsins_suggestive.fasta)"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diamond: \$(diamond --version 2>&1 | awk '{print \$3}')
        python: \$(python3 --version 2>&1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    echo ">TRINITY_DN1_trypsin_stub" > trypsins_confident.fasta
    echo "MAALGAVLLLCVLPALAARRGIPYSDGICSEVMPVKGRGKTFLVDDLCRSVAFLCGASITDVPDAMTQQIKAGKGLDEALITQNPYEGPVSEAQDLLQKLFGNASVNRFPSQARLSSQPLFIYVAGKLSSGNCATGKPIEVIDFRQKLATFHQTAAKNFFLPLGEVALQLNSTMSKVSATFGPQVSLASKLKQYAVLKPSYFNTTHSDSYTIQLTPNQFYAAASIEQGAQNIASGQNQVNQHLFPQFLDKSIRR" >> trypsins_confident.fasta
    touch trypsins_suggestive.fasta
    echo -e "id\tdiamond_hit\thmm_hit\tconfident" > identification_report.tsv
    echo "TRINITY_DN1_trypsin_stub" > trypsin_ids_confident.txt
    touch versions.yml
    """
}
