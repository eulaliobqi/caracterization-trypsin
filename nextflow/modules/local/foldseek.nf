process FOLDSEEK {
    tag        "${meta.id}"
    label      'process_medium'
    publishDir "${params.outdir}/08_validation/foldseek", mode: params.publish_dir_mode
    conda      "${projectDir}/../envs/structure.yml"

    input:
    tuple val(meta), path(pdb)

    output:
    tuple val(meta), path("${meta.id}_foldseek_results.tsv"), emit: results
    path "versions.yml",                                        emit: versions

    script:
    """
    # Busca estrutural contra PDB + AlphaFold DB
    foldseek easy-search \\
        ${pdb} \\
        /data/foldseek/pdb \\
        ${meta.id}_foldseek_results.tsv \\
        tmpfoldseek \\
        --format-output "query,target,alntmscore,rmsd,prob,bits,evalue,qlen,tlen,alnlen,qtmscore,ttmscore" \\
        --threads ${task.cpus} \\
        --exhaustive-search 0

    # Top 5 hits
    head -6 ${meta.id}_foldseek_results.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        foldseek: \$(foldseek version 2>&1)
    END_VERSIONS
    """

    stub:
    """
    echo -e "query\ttarget\talntmscore\trmsd\n${meta.id}\t1TGN_A\t0.85\t1.2" \
        > ${meta.id}_foldseek_results.tsv
    touch versions.yml
    """
}
