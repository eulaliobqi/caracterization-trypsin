// modules/local/busco.nf — BUSCO v5 para avaliação de completude do assembly
// Lineage: insecta_odb10

process BUSCO {
    tag        "${meta.id}"
    label      'process_medium'
    publishDir "${params.outdir}/01_qc/busco", mode: params.publish_dir_mode
    conda      "${projectDir}/../envs/assembly_qc.yml"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}_busco/"),                    emit: busco_dir
    tuple val(meta), path("${meta.id}_busco/short_summary*.txt"),  emit: summary
    path "versions.yml",                                            emit: versions

    script:
    def lineage_path = params.busco_lineage
    """
    busco \\
        --in ${fasta} \\
        --out ${meta.id}_busco \\
        --mode transcriptome \\
        --lineage_dataset ${lineage_path} \\
        --cpu ${task.cpus} \\
        --force

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        busco: \$(busco --version 2>&1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${meta.id}_busco
    echo "C:75.0%[S:70.0%,D:5.0%],F:10.0%,M:15.0%,n:1013" > ${meta.id}_busco/short_summary.specific.insecta_odb10.${meta.id}_busco.txt
    touch versions.yml
    """
}
