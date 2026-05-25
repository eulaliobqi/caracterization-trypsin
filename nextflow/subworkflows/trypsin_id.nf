// subworkflows/trypsin_id.nf — Fase 3: Identificação de tripsinas
// DIAMOND BLASTp (UniProt) ∩ HMMER (PF00089 Tryp_SPc)

nextflow.enable.dsl = 2

include { DIAMOND_TRYPSIN } from '../modules/local/diamond.nf'
include { HMMER_TRYPSIN   } from '../modules/local/hmmer.nf'

workflow TRYPSIN_ID {
    take:
    ch_pep    // tuple val(meta), path(proteome.pep)

    main:

    // 1. DIAMOND BLASTp contra UniProt — filtrar hits com "trypsin" no título
    DIAMOND_TRYPSIN(ch_pep)

    // 2. HMMER scan com PF00089 (Tryp_SPc domain)
    HMMER_TRYPSIN(ch_pep)

    // 3. Canal com os dois outputs para interseção (via script Python no módulo)
    ch_pep
        .join(DIAMOND_TRYPSIN.out.trypsin_ids)
        .join(HMMER_TRYPSIN.out.trypsin_ids)
        .set { ch_for_intersect }

    emit:
    trypsins_confident  = DIAMOND_TRYPSIN.out.confident_fasta  // alta confiança
    trypsins_suggestive = DIAMOND_TRYPSIN.out.suggestive_fasta  // revisão manual
    report              = DIAMOND_TRYPSIN.out.report_tsv
    versions            = DIAMOND_TRYPSIN.out.versions
}
