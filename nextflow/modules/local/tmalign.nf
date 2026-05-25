process TMALIGN {
    tag        "${meta.id}"
    label      'process_low'
    publishDir "${params.outdir}/08_validation/tmalign", mode: params.publish_dir_mode
    conda      "${projectDir}/../envs/structure.yml"

    input:
    tuple val(meta), path(pdb)

    output:
    tuple val(meta), path("${meta.id}_tmalign_scores.tsv"), emit: scores
    path "versions.yml",                                     emit: versions

    script:
    """
    # Comparar com bovine trypsin (referência)
    REF_PDB="${projectDir}/../data/references/bovine_trypsin_1TGN.pdb"

    if [[ -f "\$REF_PDB" ]]; then
        TMalign ${pdb} \$REF_PDB > tmalign_output.txt
        TM_SCORE=\$(grep "TM-score=" tmalign_output.txt | head -1 | awk -F'=' '{print \$2}' | awk '{print \$1}')
        RMSD=\$(grep "RMSD=" tmalign_output.txt | head -1 | awk -F',' '{print \$1}' | awk -F'=' '{print \$2}')
    else
        TM_SCORE="NA"
        RMSD="NA"
    fi

    echo -e "id\ttm_score\trmsd_A\treference" > ${meta.id}_tmalign_scores.tsv
    echo -e "${meta.id}\t\$TM_SCORE\t\$RMSD\t1TGN_bovine_trypsin" >> ${meta.id}_tmalign_scores.tsv

    echo "TM-align vs 1TGN: TM-score=\$TM_SCORE, RMSD=\$RMSD Å"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tmalign: \$(TMalign 2>&1 | head -2 | tail -1)
    END_VERSIONS
    """

    stub:
    """
    echo -e "id\ttm_score\trmsd_A\treference" > ${meta.id}_tmalign_scores.tsv
    echo -e "${meta.id}\t0.88\t1.35\t1TGN_bovine_trypsin" >> ${meta.id}_tmalign_scores.tsv
    touch versions.yml
    """
}
