// subworkflows/completeness.nf — Fase 4: Filtro de completude
// Met inicial + ≥220 aa + tríade catalítica His-Asp-Ser + domínio PF00089 ≥80%

nextflow.enable.dsl = 2

workflow COMPLETENESS {
    take:
    ch_trypsins   // tuple val(meta), path(trypsins_confident.fasta)

    main:

    // Módulo de filtro de completude (implementado em Python)
    // Ver nextflow/bin/filter_complete_trypsins.py
    FILTER_COMPLETE(ch_trypsins)

    emit:
    complete_trypsins = FILTER_COMPLETE.out.complete  // isoformas completas
    report            = FILTER_COMPLETE.out.report
    versions          = FILTER_COMPLETE.out.versions
}

// Processo inline para filtro de completude
process FILTER_COMPLETE {
    tag        "${meta.id}"
    label      'process_low'
    publishDir "${params.outdir}/04_complete", mode: params.publish_dir_mode
    conda      "${projectDir}/../envs/annotation.yml"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("complete_trypsins.fasta"),   emit: complete
    tuple val(meta), path("completeness_report.tsv"),   emit: report
    path "versions.yml",                                emit: versions

    script:
    """
    # Filtro 1: Met inicial + comprimento ≥ ${params.min_aa_length} aa
    python3 ${projectDir}/bin/filter_complete_trypsins.py \\
        --input ${fasta} \\
        --min_len ${params.min_aa_length} \\
        --output complete_trypsins.fasta \\
        --report completeness_report.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | awk '{print \$2}')
        biopython: \$(python3 -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch complete_trypsins.fasta
    echo -e "id\tmet_start\tlength\ttriad\tpfam_coverage\tpass" > completeness_report.tsv
    touch versions.yml
    """
}
