process IQTREE {
    tag        "${meta.id}"
    label      'process_high'
    publishDir "${params.outdir}/06_phylogeny/iqtree", mode: params.publish_dir_mode
    conda      "${projectDir}/../envs/phylogeny.yml"

    errorStrategy 'retry'
    maxRetries    2

    input:
    tuple val(meta), path(alignment)

    output:
    tuple val(meta), path("${meta.id}.treefile"),        emit: tree
    tuple val(meta), path("${meta.id}.iqtree"),          emit: log
    tuple val(meta), path("${meta.id}.contree"),         emit: contree, optional: true
    path "versions.yml",                                  emit: versions

    script:
    """
    iqtree2 \\
        -s ${alignment} \\
        -m TEST \\            # ModelFinder para seleção automática de modelo
        -bb 1000 \\           # 1000 ultrafast bootstraps
        -alrt 1000 \\         # SH-aLRT para suporte adicional
        -T ${task.cpus} \\
        -pre ${meta.id} \\
        -redo

    echo "Árvore gerada: ${meta.id}.treefile"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iqtree2: \$(iqtree2 --version 2>&1 | head -1 | awk '{print \$4}')
    END_VERSIONS
    """

    stub:
    """
    echo "((A:0.1,B:0.2):0.3,C:0.4);" > ${meta.id}.treefile
    touch ${meta.id}.iqtree
    touch versions.yml
    """
}
