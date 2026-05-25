process INTERPROSCAN {
    tag        "${meta.id}"
    label      'process_high'
    publishDir "${params.outdir}/05_primary/interproscan", mode: params.publish_dir_mode
    conda      "${projectDir}/../envs/annotation.yml"

    errorStrategy 'retry'
    maxRetries    2

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}_interproscan.tsv"),  emit: tsv
    tuple val(meta), path("${meta.id}_interproscan.xml"),  emit: xml, optional: true
    path "versions.yml",                                    emit: versions

    script:
    """
    # Remover asteriscos de stop codon do FASTA (InterProScan não aceita)
    sed 's/\\*//g' ${fasta} > ${meta.id}_clean.fasta

    interproscan.sh \\
        -i ${meta.id}_clean.fasta \\
        -f TSV,XML \\
        -appl Pfam,TIGRFAM,Gene3D,SUPERFAMILY,PANTHER,CDD \\
        -goterms -pa \\
        -cpu ${task.cpus} \\
        -o ${meta.id}_interproscan \\
        -b ${meta.id}_interproscan

    # Renomear saídas
    mv ${meta.id}_interproscan.tsv ${meta.id}_interproscan.tsv || true
    mv ${meta.id}_interproscan.xml ${meta.id}_interproscan.xml || true

    echo "IPRScan annotations: \$(wc -l < ${meta.id}_interproscan.tsv)"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        interproscan: \$(interproscan.sh --version 2>&1 | head -1 | awk '{print \$NF}')
    END_VERSIONS
    """

    stub:
    """
    echo -e "TRINITY_DN1\tMD5\t245\tPfam\tPF00089\tTryp_SPc\t1\t245\t1e-50\tT\t01-01-2026\tIPR001254\tPeptidase S1\tGO:0004252|GO:0008233" \
        > ${meta.id}_interproscan.tsv
    touch ${meta.id}_interproscan.xml
    touch versions.yml
    """
}
