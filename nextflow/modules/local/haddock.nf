process HADDOCK {
    tag        "${meta.id}_${meta.ligand}"
    label      'process_high'
    publishDir "${params.outdir}/09_docking/haddock/${meta.id}", mode: params.publish_dir_mode
    conda      "${projectDir}/../envs/docking.yml"

    input:
    tuple val(meta), path(vina_results)

    output:
    tuple val(meta), path("${meta.id}_${meta.ligand}_complex.pdb"), emit: complex_pdb
    tuple val(meta), path("${meta.id}_${meta.ligand}_haddock_score.txt"), emit: scores
    path "versions.yml",                                              emit: versions

    script:
    """
    # HADDOCK 3 — refinamento do complexo de docking
    # Usa arquivo config gerado a partir dos resultados Vina
    python3 ${projectDir}/../scripts/prepare_haddock_config.py \\
        --vina_pdbqt ${vina_results} \\
        --output haddock.cfg

    haddock3 haddock.cfg

    # Pegar melhor complexo do cluster 1
    cp haddock_output/cluster_1/haddock_best_complex.pdb \
        ${meta.id}_${meta.ligand}_complex.pdb

    echo "HADDOCK score: \$(grep 'HADDOCK score' haddock_output/cluster_1/analysis.txt)" \
        > ${meta.id}_${meta.ligand}_haddock_score.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        haddock3: \$(haddock3 --version 2>&1 || echo "3.0")
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_${meta.ligand}_complex.pdb
    echo "HADDOCK score = -125.3" > ${meta.id}_${meta.ligand}_haddock_score.txt
    touch versions.yml
    """
}
