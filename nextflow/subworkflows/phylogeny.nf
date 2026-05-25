// subworkflows/phylogeny.nf — Fase 6: Análise filogenética
// MAFFT-linsi + trimAl + IQ-TREE2 (ModelFinder + 1000 UFBoot)

nextflow.enable.dsl = 2

include { MAFFT  } from '../modules/local/mafft.nf'
include { TRIMAL } from '../modules/local/trimal.nf'
include { IQTREE } from '../modules/local/iqtree.nf'

workflow PHYLOGENY {
    take:
    ch_complete    // tuple val(meta), path(complete_trypsins.fasta)
    ch_lep_refs    // path(trypsin_refs_lepidoptera.fasta) — sequências de referência

    main:

    // 1. Combinar tripsinas A. gemmatalis + referências de Lepidoptera
    ch_complete
        .combine(ch_lep_refs)
        .map { meta, query, refs ->
            [ meta, [query, refs].flatten() ]
        }
        .set { ch_combined }

    // 2. Alinhamento múltiplo MAFFT-linsi
    MAFFT(ch_combined)

    // 3. Trimagem do alinhamento (trimAl -automated1)
    TRIMAL(MAFFT.out.alignment)

    // 4. Inferência filogenética IQ-TREE2 (ModelFinder + 1000 UFBoot)
    IQTREE(TRIMAL.out.trimmed)

    emit:
    tree       = IQTREE.out.tree       // árvore .treefile
    alignment  = TRIMAL.out.trimmed    // alinhamento aparado
    versions   = IQTREE.out.versions
}
