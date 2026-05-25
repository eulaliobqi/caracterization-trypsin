// subworkflows/md_simulation.nf — Fase 10: Simulações MD 100 ns com GROMACS 2025
// Minimização → NVT equil. → NPT equil. → Produção → Análises RMSD/RMSF/Rg/Hbonds
// REQUER GPU (label: process_gpu)

nextflow.enable.dsl = 2

include { GROMACS } from '../modules/local/gromacs.nf'

workflow MD_SIMULATION {
    take:
    ch_complexes   // tuple val(meta), path(complex.pdb)

    main:

    // GROMACS pipeline completo por complexo
    GROMACS(ch_complexes)

    emit:
    trajectories = GROMACS.out.xtc        // trajetórias .xtc
    analysis     = GROMACS.out.analysis   // RMSD, RMSF, Rg, H-bonds
    reports      = GROMACS.out.report
    versions     = GROMACS.out.versions
}
