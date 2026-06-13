#!/usr/bin/env bash
# =============================================================================
# FASE 1 — QC do Assembly
# Ferramentas: BUSCO v5 (insecta_odb10) + CD-HIT-EST (identidade 0.95)
#
# Input:  data/raw/trinity_assembly.fasta
# Output: results/phase1/
#   - busco_agemmatalis/         → relatório BUSCO completo
#   - assembly_nr95.fasta        → assembly sem redundância
#   - phase1_summary.txt         → resumo executivo
#
# Uso:
#   bash scripts/phase1_qc/run.sh
#   bash scripts/phase1_qc/run.sh --cpus 32  # override CPUS
# =============================================================================

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config.sh"

# ── Parse de argumentos opcionais ────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --cpus) CPUS="$2"; shift 2 ;;
        --input) ASSEMBLY="$2"; shift 2 ;;
        *) echo "Argumento desconhecido: $1"; exit 1 ;;
    esac
done

# ── Diretório de output ───────────────────────────────────────────────────────
OUT="${RESULTS_DIR}/phase1"
mkdir -p "$OUT"
LOG="${LOGS_DIR}/phase1_$(date +%Y%m%d_%H%M%S).log"

echo "Tee logs para: $LOG"
exec > >(tee -a "$LOG") 2>&1

log_step "FASE 1: QC do Assembly"
echo "Assembly:   $ASSEMBLY"
echo "Output:     $OUT"
echo "CPUs:       $CPUS"
echo ""

# ── Verificações de pré-requisitos ────────────────────────────────────────────
check_file "$ASSEMBLY" "Assembly Trinity" || exit 1
check_env "$ENV_DISCOVERY" || exit 1

# Contar sequências do assembly
N_SEQS=$(grep -c "^>" "$ASSEMBLY" || true)
echo "Sequências no assembly: $N_SEQS"
if [ "$N_SEQS" -lt 1000 ]; then
    echo "⚠️  AVISO: Assembly tem poucas sequências (${N_SEQS}). Esperado >10.000 para transcriptoma."
fi

# ═══════════════════════════════════════════════════════════════════════
# PASSO 1A: BUSCO — Completeness do transcriptoma
# ═══════════════════════════════════════════════════════════════════════
log_step "1A: BUSCO v5 (insecta_odb10, modo transcriptome)"

BUSCO_OUT="${OUT}/busco_agemmatalis"

if [ -f "${BUSCO_OUT}/run_${BUSCO_LINEAGE}/short_summary.specific.${BUSCO_LINEAGE}.busco_agemmatalis.txt" ]; then
    echo "⏩ BUSCO já executado — pulando. Delete ${BUSCO_OUT} para refazer."
else
    ${MAMBA_RUN} ${ENV_DISCOVERY} busco \
        --in "${ASSEMBLY}" \
        --out busco_agemmatalis \
        --out_path "${OUT}" \
        --mode transcriptome \
        --lineage_dataset "${BUSCO_LINEAGE}" \
        --cpu "${CPUS}" \
        --download_path "${REFS_DIR}/busco_downloads" \
        --force \
        2>&1 | tee "${LOGS_DIR}/busco_detail.log"
    echo "✅ BUSCO concluído"
fi

# Extrair resumo do BUSCO
BUSCO_SUMMARY=$(find "${BUSCO_OUT}" -name "short_summary*.txt" 2>/dev/null | head -1)
if [ -n "$BUSCO_SUMMARY" ]; then
    echo ""
    echo "── Resumo BUSCO ──"
    cat "$BUSCO_SUMMARY"
    cp "$BUSCO_SUMMARY" "${OUT}/busco_short_summary.txt"
fi

# ═══════════════════════════════════════════════════════════════════════
# PASSO 1B: CD-HIT-EST — Remoção de redundância
# ═══════════════════════════════════════════════════════════════════════
log_step "1B: CD-HIT-EST (c=${CDHIT_IDENTITY})"

CDHIT_PCT=$(echo "$CDHIT_IDENTITY * 100 / 1" | bc)
NR_FASTA="${OUT}/assembly_nr${CDHIT_PCT}.fasta"

if [ -f "$NR_FASTA" ]; then
    echo "⏩ CD-HIT-EST já executado — pulando. Delete ${NR_FASTA} para refazer."
else
    ${MAMBA_RUN} ${ENV_DISCOVERY} cd-hit-est \
        -i "${ASSEMBLY}" \
        -o "${NR_FASTA}" \
        -c "${CDHIT_IDENTITY}" \
        -n 10 \
        -M $((MEM_GB * 1000)) \
        -T "${CPUS}" \
        -d 0 \
        2>&1 | tee "${LOGS_DIR}/cdhit_detail.log"
    echo "✅ CD-HIT-EST concluído"
fi

N_NR=$(grep -c "^>" "$NR_FASTA" || true)
REDUNDANCY=$(echo "scale=1; (1 - ${N_NR}/${N_SEQS}) * 100" | bc)
echo ""
echo "Sequências após CD-HIT (${CDHIT_IDENTITY}): ${N_NR}"
echo "Redundância removida: ~${REDUNDANCY}%"

# ═══════════════════════════════════════════════════════════════════════
# RESUMO EXECUTIVO
# ═══════════════════════════════════════════════════════════════════════
log_step "Resumo da Fase 1"

cat > "${OUT}/phase1_summary.txt" << SUMMARY
FASE 1 — QC do Assembly
Executado em: $(date)
Assembly input: ${ASSEMBLY}
Sequências originais: ${N_SEQS}
Sequências após CD-HIT (${CDHIT_IDENTITY}): ${N_NR}
Redundância removida: ~${REDUNDANCY}%
BUSCO lineage: ${BUSCO_LINEAGE}
BUSCO output: ${BUSCO_SUMMARY}

PRÓXIMO PASSO:
  python scripts/phase1_qc/validate.py   # checar se PASS para avançar
SUMMARY

cat "${OUT}/phase1_summary.txt"
echo ""
echo "✅ Fase 1 concluída. Rode agora:"
echo "   python scripts/phase1_qc/validate.py"
