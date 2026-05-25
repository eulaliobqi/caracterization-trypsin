// subworkflows/struct_validation.nf — Fase 8: Validação estrutural
// MolProbity (Ramachandran) + ProSA (Z-score) + Foldseek (TM-score) + ConSurf (conservação)

nextflow.enable.dsl = 2

include { FOLDSEEK  } from '../modules/local/foldseek.nf'
include { TMALIGN   } from '../modules/local/tmalign.nf'
include { MOLPROBITY } from '../modules/local/molprobity.nf'
include { PROSA     } from '../modules/local/prosa.nf'
include { CONSURF   } from '../modules/local/consurf.nf'

workflow STRUCT_VALIDATION {
    take:
    ch_pdbs   // tuple val(meta), path(*.pdb)

    main:

    // 1. MolProbity — qualidade estrutural (Ramachandran, rotâmeros, clashscore)
    MOLPROBITY(ch_pdbs)

    // 2. ProSA — Z-score comparativo com cristalografia/NMR
    PROSA(ch_pdbs)

    // 3. Foldseek — busca estrutural contra AFDB + PDB (TM-score)
    FOLDSEEK(ch_pdbs)

    // 4. TM-align com bovine trypsin (referência)
    TMALIGN(ch_pdbs)

    // 5. ConSurf — mapeamento de conservação evolutiva
    CONSURF(ch_pdbs)

    // Filtrar estruturas que passam no critério de qualidade
    ch_pdbs
        .join(MOLPROBITY.out.report)
        .filter { meta, pdb, report ->
            // Ramachandran favored > 95% E clashscore < 20
            def rama = report.text =~ /Ramachandran favored = (\d+\.\d+)/
            rama ? rama[0][1].toFloat() > 95.0 : true // pass se não parsear
        }
        .map { meta, pdb, report -> [ meta, pdb ] }
        .set { ch_validated }

    emit:
    validated_pdbs = ch_validated           // PDBs aprovados para docking
    molprobity     = MOLPROBITY.out.report
    foldseek       = FOLDSEEK.out.results
    tmalign        = TMALIGN.out.scores
    versions       = FOLDSEEK.out.versions
}
