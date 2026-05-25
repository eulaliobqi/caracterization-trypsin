process SIGNALP {
    tag        "${meta.id}"
    label      'process_medium'
    publishDir "${params.outdir}/05_primary/signalp", mode: params.publish_dir_mode
    conda      "${projectDir}/../envs/annotation.yml"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}_signalp6.tsv"),   emit: tsv
    tuple val(meta), path("${meta.id}_signalp6.gff3"),  emit: gff3, optional: true
    path "versions.yml",                                  emit: versions

    script:
    """
    # SignalP 6 (modo fast para insetos - Eukarya)
    signalp6 \\
        --fastafile ${fasta} \\
        --organism eukarya \\
        --output_dir signalp_out \\
        --format none \\
        --mode fast

    # Consolidar resultados
    mv signalp_out/output.gff3    ${meta.id}_signalp6.gff3 2>/dev/null || true
    mv signalp_out/prediction_results.txt ${meta.id}_signalp6.tsv

    echo "SignalP SP predictions: \$(grep -c 'SP(' ${meta.id}_signalp6.tsv || echo 0)"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        signalp: \$(signalp6 --version 2>&1 || echo "6.0")
    END_VERSIONS
    """

    stub:
    """
    echo -e "# SignalP-6.0 stub output" > ${meta.id}_signalp6.tsv
    echo -e "ID\tPrediction\tSP(Sec/SPI)\tOTHER\tCS Position" >> ${meta.id}_signalp6.tsv
    echo -e "TRINITY_DN1\tSP(Sec/SPI)\t0.98\t0.02\t18-19" >> ${meta.id}_signalp6.tsv
    touch versions.yml
    """
}
