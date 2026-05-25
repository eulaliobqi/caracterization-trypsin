// subworkflows/qc_assembly.nf — Fase 1: QC do assembly Trinity
// BUSCO v5 (insecta_odb10) + CD-HIT-EST (0.95) + seqkit stats

nextflow.enable.dsl = 2

include { BUSCO       } from '../modules/local/busco.nf'
include { CD_HIT      } from '../modules/local/cd_hit.nf'

workflow QC_ASSEMBLY {
    take:
    ch_fasta   // tuple val(meta), path(fasta)

    main:

    // 1. Estatísticas do assembly bruto
    ch_fasta
        .map { meta, fasta -> [ meta, fasta ] }
        .set { ch_for_busco }

    // 2. BUSCO para avaliar completude gênica (insecta_odb10)
    BUSCO(ch_for_busco)

    // 3. CD-HIT-EST para reduzir redundância de transcritos (95%)
    CD_HIT(ch_fasta)

    emit:
    clean_fasta  = CD_HIT.out.fasta      // tuple val(meta), path(fasta)
    busco_report = BUSCO.out.summary     // tuple val(meta), path(summary.txt)
    versions     = BUSCO.out.versions
}
