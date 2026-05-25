// subworkflows/docking.nf — Fase 9: Docking molecular
// AutoDock Vina (screening) + HADDOCK (refinamento) + PLIP (análise de interações)

nextflow.enable.dsl = 2

include { AUTODOCK_VINA } from '../modules/local/autodock_vina.nf'
include { HADDOCK       } from '../modules/local/haddock.nf'

workflow DOCKING {
    take:
    ch_validated_pdbs   // tuple val(meta), path(validated.pdb)

    main:

    // Canal de inibidores (SKTI + BPTI-like peptides)
    ch_inhibitors = Channel.fromPath(params.inhibitors_fasta, checkIfExists: true)
        .splitFasta(by: 1, file: true)

    // Combinação: cada tripsina x cada inibidor
    ch_validated_pdbs
        .combine(ch_inhibitors)
        .map { meta, receptor, ligand ->
            def ligand_id = ligand.name.replaceFirst(/\.fasta$/, '')
            def new_meta = meta + [ligand: ligand_id]
            [ new_meta, receptor, ligand ]
        }
        .set { ch_docking_pairs }

    // 1. AutoDock Vina — screening inicial (mais rápido)
    AUTODOCK_VINA(ch_docking_pairs)

    // 2. Selecionar top 3 complexos por tripsina (ΔG mais negativo)
    AUTODOCK_VINA.out.results
        .groupTuple(by: 0)
        .flatMap { meta, results ->
            results.sort { a, b ->
                // Ordenar por score ΔG (menor = melhor)
                a.text.find(/-\d+\.\d+/)?.toFloat() <=>
                b.text.find(/-\d+\.\d+/)?.toFloat()
            }.take(3).collect { r -> [ meta, r ] }
        }
        .set { ch_top_vina }

    // 3. HADDOCK — refinamento dos top complexos
    HADDOCK(ch_top_vina)

    emit:
    top_complexes = HADDOCK.out.complex_pdb   // top complexos para MD
    vina_scores   = AUTODOCK_VINA.out.scores
    haddock_scores = HADDOCK.out.scores
    versions      = AUTODOCK_VINA.out.versions
}
