#!/usr/bin/env nextflow
// main.nf — Pipeline principal para caracterização estrutural de tripsinas
// Organismo: Anticarsia gemmatalis (Lepidoptera: Erebidae)
// Autor: Eulálio Gutemberg
// Versão: 1.0.0

nextflow.enable.dsl = 2

// ── Imports de subworkflows ───────────────────────────────────────────────
include { QC_ASSEMBLY       } from './subworkflows/qc_assembly.nf'
include { ORF_PREDICTION    } from './subworkflows/orf_prediction.nf'
include { TRYPSIN_ID        } from './subworkflows/trypsin_id.nf'
include { COMPLETENESS      } from './subworkflows/completeness.nf'
include { PRIMARY_CHAR      } from './subworkflows/primary_char.nf'
include { PHYLOGENY         } from './subworkflows/phylogeny.nf'
include { STRUCTURE_PRED    } from './subworkflows/structure_pred.nf'
include { STRUCT_VALIDATION } from './subworkflows/struct_validation.nf'
include { DOCKING           } from './subworkflows/docking.nf'
include { MD_SIMULATION     } from './subworkflows/md_simulation.nf'
include { REPORT_GEN        } from './subworkflows/report_gen.nf'

// ── Help message ─────────────────────────────────────────────────────────
def helpMsg() {
    log.info """
    ╔══════════════════════════════════════════════════════════════╗
    ║    TRYPSIN-AGEMMATALIS-STRUCTURAL  v${workflow.manifest.version}           ║
    ║    Anticarsia gemmatalis trypsin characterization pipeline   ║
    ╚══════════════════════════════════════════════════════════════╝

    Usage:
      nextflow run main.nf [options]

    Required:
      --input_fasta   PATH    Trinity assembly FASTA (default: data/raw/trinity_assembly.fasta)
      --outdir        PATH    Output directory (default: results/)

    Phase control (--phase N runs phases 1 through N):
      --phase 1   QC assembly (BUSCO + CD-HIT-EST)
      --phase 2   ORF prediction (TransDecoder)
      --phase 3   Trypsin identification (DIAMOND + HMMER)
      --phase 4   Completeness filter (catalytic triad)
      --phase 5   Primary characterization (ProtParam + SignalP6)
      --phase 6   Phylogeny (MAFFT + IQ-TREE2)
      --phase 7   Structure prediction (AlphaFold3) [GPU]
      --phase 8   Structural validation (MolProbity + Foldseek)
      --phase 9   Docking (AutoDock Vina + HADDOCK)
      --phase 10  MD simulations 100ns (GROMACS) [GPU]
      --phase 11  Report and manuscript generation
      --phase null Run all phases (default)

    Profiles:
      -profile test     Minimal test (use with -stub-run)
      -profile debian   Debian server 64 cores, 256 GB RAM
      -profile gpu      Add GPU support (combine: -profile debian,gpu)

    Examples:
      # Run all phases:
      nextflow run main.nf -profile debian --input_fasta data/raw/trinity_assembly.fasta

      # Run only QC (phase 1):
      nextflow run main.nf -profile debian --phase 1 --input_fasta data/raw/trinity_assembly.fasta

      # Resume after interruption:
      nextflow run main.nf -profile debian -resume

      # Stub-run (test without executing tools):
      nextflow run main.nf -profile test -stub-run
    """.stripIndent()
}

if (params.help) {
    helpMsg()
    exit 0
}

// ── Validação de inputs ───────────────────────────────────────────────────
def validateParams() {
    // Verificar FASTA de input
    if (!file(params.input_fasta).exists()) {
        log.error "❌ Input FASTA não encontrado: ${params.input_fasta}"
        log.error "   Coloque seu assembly Trinity em data/raw/trinity_assembly.fasta"
        exit 1
    }
    // Avisos sobre referências
    if (!file(params.uniprot_db).exists()) {
        log.warn "⚠️  UniProt DIAMOND DB não encontrado: ${params.uniprot_db}"
        log.warn "   Baixe e crie com: diamond makedb --in uniprot_sprot.fasta -d uniprot_sprot"
    }
}
validateParams()

// ── Workflow principal ────────────────────────────────────────────────────
workflow {

    // Log de início
    log.info """
    ╔══════════════════════════════════════════════════════╗
    ║   Trypsin-Agemmatalis-Structural v${workflow.manifest.version}         ║
    ╚══════════════════════════════════════════════════════╝
    input_fasta : ${params.input_fasta}
    outdir      : ${params.outdir}
    phase       : ${params.phase ?: 'all (1-11)'}
    profile     : ${workflow.profile}
    """.stripIndent()

    // Canal de input: meta map + FASTA
    ch_input = Channel.fromPath(params.input_fasta, checkIfExists: true)
        .map { fasta -> [ [id: 'agemmatalis_midgut', organism: 'Anticarsia_gemmatalis'], fasta ] }

    // ── Fase 1: QC do assembly ────────────────────────────────────────────
    if (!params.phase || params.phase as int >= 1) {
        QC_ASSEMBLY(ch_input)
        ch_clean_fasta = QC_ASSEMBLY.out.clean_fasta
    } else {
        ch_clean_fasta = ch_input
    }

    // ── Fase 2: Predição de ORFs ──────────────────────────────────────────
    if (!params.phase || params.phase as int >= 2) {
        ORF_PREDICTION(ch_clean_fasta)
        ch_pep = ORF_PREDICTION.out.pep
    }

    // ── Fase 3: Identificação de tripsinas ────────────────────────────────
    if (!params.phase || params.phase as int >= 3) {
        TRYPSIN_ID(ch_pep)
        ch_trypsins_confident = TRYPSIN_ID.out.trypsins_confident
    }

    // ── Fase 4: Filtro de completude ──────────────────────────────────────
    if (!params.phase || params.phase as int >= 4) {
        COMPLETENESS(ch_trypsins_confident)
        ch_complete = COMPLETENESS.out.complete_trypsins
    }

    // ── Fase 5: Caracterização primária ───────────────────────────────────
    if (!params.phase || params.phase as int >= 5) {
        PRIMARY_CHAR(ch_complete)
    }

    // ── Fase 6: Filogenia ─────────────────────────────────────────────────
    if (!params.phase || params.phase as int >= 6) {
        ch_lep_refs = Channel.fromPath(params.lep_refs, checkIfExists: true)
        PHYLOGENY(ch_complete, ch_lep_refs)
    }

    // ── Fase 7: Predição estrutural (AlphaFold3) ──────────────────────────
    if (!params.phase || params.phase as int >= 7) {
        STRUCTURE_PRED(ch_complete)
        ch_pdbs = STRUCTURE_PRED.out.pdbs
    }

    // ── Fase 8: Validação estrutural ──────────────────────────────────────
    if (!params.phase || params.phase as int >= 8) {
        STRUCT_VALIDATION(ch_pdbs)
        ch_validated_pdbs = STRUCT_VALIDATION.out.validated_pdbs
    }

    // ── Fase 9: Docking ───────────────────────────────────────────────────
    if (!params.phase || params.phase as int >= 9) {
        DOCKING(ch_validated_pdbs)
        ch_top_complexes = DOCKING.out.top_complexes
    }

    // ── Fase 10: Simulações MD ────────────────────────────────────────────
    if (!params.phase || params.phase as int >= 10) {
        MD_SIMULATION(ch_top_complexes)
        ch_md_results = MD_SIMULATION.out.analysis
    }

    // ── Fase 11: Relatório e manuscrito ───────────────────────────────────
    if (!params.phase || params.phase as int >= 11) {
        REPORT_GEN(
            QC_ASSEMBLY.out,
            ORF_PREDICTION.out,
            TRYPSIN_ID.out,
            COMPLETENESS.out,
            PRIMARY_CHAR.out,
            PHYLOGENY.out,
            STRUCTURE_PRED.out,
            STRUCT_VALIDATION.out,
            DOCKING.out,
            MD_SIMULATION.out
        )
    }
}

// ── Handlers de eventos ───────────────────────────────────────────────────
workflow.onComplete {
    def status = workflow.success ? '✅ SUCCESS' : '❌ FAILED'
    log.info """
    ══════════════════════════════════════════════════
    Pipeline finalizado: ${status}
    Duração total:  ${workflow.duration}
    Output em:      ${params.outdir}
    CPU hours:      ${workflow.stats.computeTimeFmt ?: 'N/A'}
    ══════════════════════════════════════════════════
    """.stripIndent()

    if (!workflow.success) {
        log.error "Verifique .nextflow.log e work/ para detalhes do erro"
        log.error "Consulte LEARNINGS.md para erros conhecidos"
    }
}

workflow.onError {
    log.error "Pipeline ERRO: ${workflow.errorMessage}"
    log.error "Processo com falha: ${workflow.errorReport}"
}
