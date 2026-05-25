process AUTODOCK_VINA {
    tag        "${meta.id}_${meta.ligand}"
    label      'process_medium'
    publishDir "${params.outdir}/09_docking/vina/${meta.id}", mode: params.publish_dir_mode
    conda      "${projectDir}/../envs/docking.yml"

    input:
    tuple val(meta), path(receptor_pdb), path(ligand_fasta)

    output:
    tuple val(meta), path("${meta.id}_${meta.ligand}_vina.pdbqt"),   emit: complex_pdbqt
    tuple val(meta), path("${meta.id}_${meta.ligand}_vina_score.txt"), emit: results
    tuple val(meta), path("${meta.id}_${meta.ligand}_score.txt"),      emit: scores
    path "versions.yml",                                                emit: versions

    script:
    """
    # Preparar receptor (PDB → PDBQT)
    python3 ${projectDir}/../scripts/prepare_receptor.py \\
        -r ${receptor_pdb} \\
        -o receptor.pdbqt

    # Preparar ligante peptídico (FASTA → PDB → PDBQT)
    python3 ${projectDir}/../scripts/prepare_peptide_ligand.py \\
        --fasta ${ligand_fasta} \\
        --output ligand.pdbqt

    # Definir caixa de docking centrada no sítio ativo da tripsina
    # (aproximado ao redor do Ser195 conservado)
    python3 ${projectDir}/../scripts/define_binding_site.py \\
        --pdb ${receptor_pdb} \\
        --residue SER195 \\
        --output box_params.txt

    # Executar AutoDock Vina
    vina \\
        --receptor receptor.pdbqt \\
        --ligand ligand.pdbqt \\
        \$(cat box_params.txt) \\
        --exhaustiveness 32 \\
        --num_modes 9 \\
        --out ${meta.id}_${meta.ligand}_vina.pdbqt \\
        --log ${meta.id}_${meta.ligand}_vina_score.txt \\
        --cpu ${task.cpus}

    # Extrair melhor score ΔG
    grep "^   1 " ${meta.id}_${meta.ligand}_vina_score.txt | awk '{print \$2}' \
        > ${meta.id}_${meta.ligand}_score.txt
    echo "Vina ΔG: \$(cat ${meta.id}_${meta.ligand}_score.txt) kcal/mol"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vina: \$(vina --version 2>&1 | head -1)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_${meta.ligand}_vina.pdbqt
    echo "   1    -9.2" > ${meta.id}_${meta.ligand}_vina_score.txt
    echo "-9.2" > ${meta.id}_${meta.ligand}_score.txt
    touch versions.yml
    """
}
