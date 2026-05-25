// subworkflows/orf_prediction.nf — Fase 2: Predição de ORFs com TransDecoder
// TransDecoder.LongOrfs → BLAST hint → Pfam hint → TransDecoder.Predict

nextflow.enable.dsl = 2

include { TRANSDECODER } from '../modules/local/transdecoder.nf'
include { DIAMOND      } from '../modules/local/diamond.nf'
include { HMMER        } from '../modules/local/hmmer.nf'

workflow ORF_PREDICTION {
    take:
    ch_fasta   // tuple val(meta), path(clean_fasta)

    main:

    // 1. DIAMOND BLASTp para hints de ORF (opcional mas melhora predição)
    DIAMOND(ch_fasta)

    // 2. HMMER para hints Pfam (detecta domínios conservados)
    HMMER(ch_fasta)

    // 3. TransDecoder com hints combinados
    TRANSDECODER(
        ch_fasta,
        DIAMOND.out.tsv,
        HMMER.out.domtblout
    )

    emit:
    pep      = TRANSDECODER.out.pep    // tuple val(meta), path(*.pep)
    cds      = TRANSDECODER.out.cds    // tuple val(meta), path(*.cds)
    gff3     = TRANSDECODER.out.gff3   // tuple val(meta), path(*.gff3)
    versions = TRANSDECODER.out.versions
}
