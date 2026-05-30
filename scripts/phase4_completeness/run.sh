#!/usr/bin/env bash
# =============================================================================
# FASE 4 — Filtro de Completude das Tripsinas
#
# Critérios: Met inicial + ≥220 aa + motivo GDSGG (Ser195)
# Opcional:  tríade His57-Asp102-Ser195 via alinhamento com tripsina bovina
#
# Input:  results/phase3/trypsins_confident.fasta
# Output: results/phase4/
#   - complete_trypsins.fasta    → tripsinas completas (para modelagem)
#   - borderline_trypsins.fasta  → Met+comprimento OK mas sem motivo GDSGG
#   - completeness_report.tsv   → detalhes por sequência
# =============================================================================

# ── Ativar ambiente conda ANTES do set -euo pipefail ─────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
_CB="$(conda info --base 2>/dev/null || echo "${HOME}/miniforge3")"
# shellcheck disable=SC1091
source "${_CB}/etc/profile.d/conda.sh" 2>/dev/null || true
source "${_CB}/etc/profile.d/mamba.sh" 2>/dev/null || true
if ! conda activate orf_prediction 2>/dev/null; then
    echo "❌ Ambiente 'orf_prediction' não encontrado."
    exit 1
fi

set -euo pipefail
source "${SCRIPT_DIR}/../config.sh"

log_step "FASE 4: Filtro de Completude"

CONFIDENT="${RESULTS_DIR}/phase3/trypsins_confident.fasta"
SUGGESTIVE="${RESULTS_DIR}/phase3/trypsins_suggestive.fasta"
OUT="${RESULTS_DIR}/phase4"
mkdir -p "$OUT"
LOG="${LOGS_DIR}/phase4_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1

check_file "$CONFIDENT" "trypsins_confident.fasta (Fase 3)" || exit 1

N_CONF=$(grep -c "^>" "$CONFIDENT")
echo "Confident da Fase 3: $N_CONF sequências"
echo ""

# ── Filtro principal ──────────────────────────────────────────────────────────
log_step "4A: Filtro Met + ≥220 aa + motivo GDSGG"

python3 "${SCRIPT_DIR}/filter.py" \
    --input    "$CONFIDENT" \
    --results4 "$OUT" \
    --min-len  "${MIN_TRYPSIN_AA}"

# ── Opcional: incluir suggestive se poucas tripsinas completas ────────────────
N_COMPLETE=$(grep -c "^>" "${OUT}/complete_trypsins.fasta" 2>/dev/null || echo 0)
if [ "$N_COMPLETE" -lt 5 ] && [ -f "$SUGGESTIVE" ]; then
    echo ""
    echo "⚠️  Apenas $N_COMPLETE completas em confident — expandindo para suggestive..."
    python3 "${SCRIPT_DIR}/filter.py" \
        --input    "$CONFIDENT" \
        --results4 "$OUT" \
        --min-len  "${MIN_TRYPSIN_AA}" \
        --also-suggestive
fi

# ── Resumo ────────────────────────────────────────────────────────────────────
log_step "Resumo da Fase 4"

N_FINAL=$(grep -c "^>" "${OUT}/complete_trypsins.fasta" 2>/dev/null || echo 0)
N_BORDER=$(grep -c "^>" "${OUT}/borderline_trypsins.fasta" 2>/dev/null || echo 0)

echo "Tripsinas completas:  $N_FINAL"
echo "Borderline:           $N_BORDER"
echo ""
echo "✅ Fase 4 concluída. Próximo passo:"
echo "   python scripts/phase4_completeness/validate.py"
