#!/usr/bin/env bash
# =============================================================================
# config.sh — Configuração central do pipeline de tripsinas
# Source este arquivo no início de cada script: source "$(dirname $0)/../config.sh"
# =============================================================================

# ── Diretórios ────────────────────────────────────────────────────────────────
PROJ_ROOT="${HOME}/gromacs/caracterization-trypsin"
DATA_DIR="${PROJ_ROOT}/data"
RAW_DIR="${DATA_DIR}/raw"
REFS_DIR="${DATA_DIR}/references"
RESULTS_DIR="${PROJ_ROOT}/results"
BIN_DIR="${PROJ_ROOT}/nextflow/bin"
LOGS_DIR="${PROJ_ROOT}/logs"

# ── Input principal ───────────────────────────────────────────────────────────
ASSEMBLY="${RAW_DIR}/trinity_assembly.fasta"

# ── Referências (bancos de dados no servidor) ─────────────────────────────────
PROTEIN_DB="/home/eulalio/databases/nr/nr.dmnd"          # NCBI NR (DIAMOND)
UNIPROT_DB="$PROTEIN_DB"                                  # alias para compatibilidade
PFAM_HMM="/home/eulalio/databases/pfam/Pfam-A.hmm"       # já pressionado (h3f/h3i/h3m/h3p)
EGGNOG_DB="/home/eulalio/databases/eggnog/eggnog.db"      # EggNOG (Fase 5)
EGGNOG_PROTEINS="/home/eulalio/databases/eggnog/eggnog_proteins.dmnd"
TAXDUMP_DIR="/home/eulalio/databases/taxdump"
BUSCO_LINEAGE="insecta_odb10"                             # baixado automaticamente pelo BUSCO
BOVINE_TRYPSIN="${REFS_DIR}/bovine_trypsin_1TGN.fasta"    # para alinhamento da tríade
TRYP_REFS="${REFS_DIR}/trypsin_refs_lepidoptera.fasta"

# ── Parâmetros de qualidade ───────────────────────────────────────────────────
CPUS=16                    # ajuste conforme seu servidor (use nproc para detectar)
MEM_GB=64                  # memória disponível em GB
MIN_ORF_AA=100             # mínimo de aa para TransDecoder.LongOrfs
DIAMOND_EVALUE="1e-10"     # e-value para DIAMOND vs UniProt
HMMER_EVALUE="1e-10"       # e-value para HMMER vs PF00089
MIN_TRYPSIN_AA=220         # comprimento mínimo de tripsina completa
CDHIT_IDENTITY=0.95        # identidade para CD-HIT-EST (Fase 1)
PFAM_COVERAGE=0.80         # cobertura mínima do domínio PF00089

# ── Ambientes conda ───────────────────────────────────────────────────────────
ENV_DISCOVERY="discovery"  # Fases 1-4
ENV_ANALYSIS="analysis"    # Fase 5 (ProtParam)
ENV_PHYLOGENY="phylogeny"  # Fase 6
ENV_STRUCTURE="structure"  # Fase 7-8
ENV_DOCKING="docking"      # Fase 9

# ── Execução em ambiente conda ───────────────────────────────────────────────
# Injeta bin/ do env no PATH — evita mamba run (gera exec -- incompatível com bash)
# e conda activate (não funciona em contexto não-interativo dentro de função).
_CONDA_BASE="$(conda info --base 2>/dev/null || echo "${HOME}/miniforge3")"

mamba_run() {
    local env="$1"; shift
    local cmd="$1"; shift
    local env_bin="${_CONDA_BASE}/envs/${env}/bin"
    if [ ! -d "$env_bin" ]; then
        echo "❌ ERRO: ambiente conda '${env}' não encontrado em ${env_bin}"
        return 1
    fi
    PATH="${env_bin}:${PATH}" "${env_bin}/${cmd}" "$@"
}
MAMBA_RUN="mamba_run"

# ── Funções utilitárias ───────────────────────────────────────────────────────
log_step() {
    echo ""
    echo "════════════════════════════════════════"
    echo "  $1"
    echo "  $(date '+%Y-%m-%d %H:%M:%S')"
    echo "════════════════════════════════════════"
}

check_file() {
    local f="$1"
    local desc="$2"
    if [ ! -f "$f" ]; then
        echo "❌ ERRO: $desc não encontrado: $f"
        return 1
    fi
    echo "✅ $desc: $f"
    return 0
}

check_env() {
    local env="$1"
    local conda_base
    conda_base=$(conda info --base 2>/dev/null || echo "${HOME}/miniforge3")
    if [ ! -d "${conda_base}/envs/${env}" ]; then
        echo "❌ ERRO: ambiente conda '${env}' não encontrado"
        echo "   Rode: mamba env create -f envs/${env}.yml --yes"
        return 1
    fi
    return 0
}

# Criar diretório de logs se não existir
mkdir -p "${LOGS_DIR}"
