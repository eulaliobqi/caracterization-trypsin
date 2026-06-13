// subworkflows/trypsin_id.nf — Fase 3: Identificação de tripsinas
// Pipeline: DIAMOND BLASTp → HMMER (PF00089) → INTERSECT (Python)
// Confident = presente nos DOIS métodos | Suggestive = apenas UM

nextflow.enable.dsl = 2

include { DIAMOND_TRYPSIN   } from '../modules/local/diamond.nf'
include { HMMER_TRYPSIN     } from '../modules/local/hmmer.nf'
include { INTERSECT_TRYPSINS } from '../modules/local/diamond.nf'

workflow TRYPSIN_ID {
    take:
    ch_pep    // tuple val(meta), path(proteome.pep)

    main:

    // Passo 1 — DIAMOND BLASTp: busca contra UniProt, filtra hits "trypsin"
    // Produz: diamond_trypsin_ids.txt
    DIAMOND_TRYPSIN(ch_pep)

    // Passo 2 — HMMER scan: busca pelo domínio Pfam PF00089 (Tryp_SPc)
    // Produz: trypsin_hmm_ids.txt
    HMMER_TRYPSIN(ch_pep)

    // Passo 3 — Interseção real: DIAMOND ∩ HMMER via script Python
    // confident = ambos concordam (alta especificidade)
    // suggestive = apenas um método (requer revisão manual)
    ch_pep
        .join( DIAMOND_TRYPSIN.out.trypsin_ids )
        .join( HMMER_TRYPSIN.out.trypsin_ids   )
        .set { ch_for_intersect }

    INTERSECT_TRYPSINS(ch_for_intersect)

    emit:
    trypsins_confident  = INTERSECT_TRYPSINS.out.confident_fasta   // alta confiança
    trypsins_suggestive = INTERSECT_TRYPSINS.out.suggestive_fasta   // revisão manual
    report              = INTERSECT_TRYPSINS.out.report_tsv
    versions            = DIAMOND_TRYPSIN.out.versions
}
