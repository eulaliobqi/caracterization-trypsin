// subworkflows/structure_pred.nf — Fase 7: Predição estrutural com AlphaFold3
// Via módulo local que chama run_alphafold3.py
// REQUER GPU (label: process_gpu)

nextflow.enable.dsl = 2

include { ALPHAFOLD3 } from '../modules/local/alphafold3.nf'

workflow STRUCTURE_PRED {
    take:
    ch_complete   // tuple val(meta), path(complete_trypsins.fasta)

    main:

    // Dividir FASTA em sequências individuais para paralelizar AF3
    ch_complete
        .splitFasta(by: 1, file: true, decompress: false)
        .map { fasta ->
            def id = fasta.name.replaceFirst(/\.fasta$/, '')
            [ [id: id, organism: 'Anticarsia_gemmatalis'], fasta ]
        }
        .set { ch_individual_seqs }

    // AlphaFold3 por sequência (paralelismo automático via Nextflow)
    ALPHAFOLD3(ch_individual_seqs)

    emit:
    pdbs       = ALPHAFOLD3.out.pdb        // tuple val(meta), path(*.pdb)
    confidence = ALPHAFOLD3.out.confidence  // tuple val(meta), path(confidence.json)
    versions   = ALPHAFOLD3.out.versions
}
