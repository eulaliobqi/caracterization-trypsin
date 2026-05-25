process MAFFT {
    tag        "${meta.id}"
    label      'process_medium'
    publishDir "${params.outdir}/06_phylogeny/mafft", mode: params.publish_dir_mode
    conda      "${projectDir}/../envs/phylogeny.yml"

    input:
    tuple val(meta), path(fastas)

    output:
    tuple val(meta), path("${meta.id}_aligned.fasta"), emit: alignment
    path "versions.yml",                                emit: versions

    script:
    """
    cat ${fastas} > combined.fasta

    mafft --localpair --maxiterate 1000 --thread ${task.cpus} \\
        combined.fasta > ${meta.id}_aligned.fasta

    echo "Sequências alinhadas: \$(grep -c '>' ${meta.id}_aligned.fasta)"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mafft: \$(mafft --version 2>&1)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_aligned.fasta
    touch versions.yml
    """
}
