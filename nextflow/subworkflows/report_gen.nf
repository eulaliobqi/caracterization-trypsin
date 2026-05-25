// subworkflows/report_gen.nf — Fase 11: Geração de relatório e manuscrito
// Compila todos os resultados em figuras, tabelas e rascunho de manuscrito

nextflow.enable.dsl = 2

workflow REPORT_GEN {
    take:
    ch_qc
    ch_orfs
    ch_trypsins
    ch_complete
    ch_primary
    ch_phylogeny
    ch_structures
    ch_validation
    ch_docking
    ch_md

    main:

    // Gerar figuras e tabelas
    MAKE_FIGURES(
        ch_complete,
        ch_primary,
        ch_phylogeny,
        ch_structures,
        ch_validation,
        ch_docking,
        ch_md
    )

    // Gerar relatório HTML
    GENERATE_REPORT(
        ch_qc,
        ch_complete,
        ch_primary,
        MAKE_FIGURES.out.figures
    )

    emit:
    report  = GENERATE_REPORT.out.html
    figures = MAKE_FIGURES.out.figures
}

process MAKE_FIGURES {
    tag        "figures"
    label      'process_low'
    publishDir "${params.outdir}/11_report/figures", mode: params.publish_dir_mode
    conda      "${projectDir}/../envs/analysis.yml"

    input:
    tuple val(meta_complete),    path(complete_fasta)
    tuple val(meta_primary),     path(protparam_tsv)
    tuple val(meta_phylogeny),   path(tree_file)
    tuple val(meta_struct),      path(pdb_files)
    tuple val(meta_val),         path(molprobity_report)
    tuple val(meta_dock),        path(docking_scores)
    tuple val(meta_md),          path(md_analysis)

    output:
    path "*.svg",   emit: figures
    path "*.png",   emit: figures_png, optional: true
    path "versions.yml", emit: versions

    script:
    """
    python3 ${projectDir}/bin/make_figures.py \\
        --protparam ${protparam_tsv} \\
        --tree ${tree_file} \\
        --docking ${docking_scores} \\
        --md_analysis ${md_analysis} \\
        --outdir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | awk '{print \$2}')
        matplotlib: \$(python3 -c "import matplotlib; print(matplotlib.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch fig1_trypsin_properties.svg
    touch fig2_phylogeny.svg
    touch fig3_structures_plddt.svg
    touch fig4_docking.svg
    touch fig5_md_rmsd.svg
    touch versions.yml
    """
}

process GENERATE_REPORT {
    tag        "report"
    label      'process_low'
    publishDir "${params.outdir}/11_report", mode: params.publish_dir_mode
    conda      "${projectDir}/../envs/analysis.yml"

    input:
    tuple val(meta_qc),       path(busco_summary)
    tuple val(meta_complete), path(complete_fasta)
    tuple val(meta_primary),  path(protparam_tsv)
    path figures

    output:
    path "pipeline_report.html", emit: html
    path "manuscript_draft.md",  emit: manuscript, optional: true
    path "versions.yml",          emit: versions

    script:
    """
    python3 -c "
    import os, datetime
    n_trypsins = sum(1 for line in open('${complete_fasta}') if line.startswith('>'))
    report = f'''# Pipeline Report — Trypsin Agemmatalis Structural
    Generated: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M')}
    Complete trypsins identified: {n_trypsins}
    Assembly: ${params.input_fasta}
    '''
    open('pipeline_report.html', 'w').write(f'<html><body><pre>{report}</pre></body></html>')
    print(f'Report generated with {n_trypsins} trypsins')
    "

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    echo '<html><body>Report stub</body></html>' > pipeline_report.html
    touch manuscript_draft.md
    touch versions.yml
    """
}
