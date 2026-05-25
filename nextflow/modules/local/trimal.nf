process TRIMAL {
    tag        "${meta.id}"
    label      'process_low'
    publishDir "${params.outdir}/06_phylogeny/trimal", mode: params.publish_dir_mode
    conda      "${projectDir}/../envs/phylogeny.yml"

    input:
    tuple val(meta), path(alignment)

    output:
    tuple val(meta), path("${meta.id}_trimmed.fasta"), emit: trimmed
    path "versions.yml",                                emit: versions

    script:
    """
    trimal -in ${alignment} -out ${meta.id}_trimmed.fasta -automated1

    echo "Colunas após trimming: \$(awk 'NR==2{print length(\$0)}' ${meta.id}_trimmed.fasta)"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trimal: \$(trimal --version 2>&1 | head -1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_trimmed.fasta
    touch versions.yml
    """
}
