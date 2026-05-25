// subworkflows/primary_char.nf — Fase 5: Caracterização primária
// ProtParam (MW, pI, GRAVY) + SignalP6 + InterProScan

nextflow.enable.dsl = 2

include { INTERPROSCAN } from '../modules/local/interproscan.nf'
include { SIGNALP      } from '../modules/local/signalp.nf'

workflow PRIMARY_CHAR {
    take:
    ch_complete   // tuple val(meta), path(complete_trypsins.fasta)

    main:

    // 1. ProtParam — propriedades físico-químicas (MW, pI, GRAVY, instabilidade)
    PROTPARAM(ch_complete)

    // 2. SignalP 6 — predição de peptídeo sinal N-terminal
    SIGNALP(ch_complete)

    // 3. InterProScan — domínios, famílias, GO terms
    INTERPROSCAN(ch_complete)

    emit:
    protparam     = PROTPARAM.out.tsv
    signalp       = SIGNALP.out.tsv
    interproscan  = INTERPROSCAN.out.tsv
    versions      = INTERPROSCAN.out.versions
}

process PROTPARAM {
    tag        "${meta.id}"
    label      'process_low'
    publishDir "${params.outdir}/05_primary/protparam", mode: params.publish_dir_mode
    conda      "${projectDir}/../envs/annotation.yml"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("protparam_results.tsv"), emit: tsv
    path "versions.yml",                            emit: versions

    script:
    """
    python3 ${projectDir}/bin/compute_protparam.py \\
        --input ${fasta} \\
        --output protparam_results.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | awk '{print \$2}')
        biopython: \$(python3 -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """

    stub:
    """
    echo -e "id\tmw_kda\tpi\tgravy\tinstability\taa_count" > protparam_results.tsv
    touch versions.yml
    """
}
