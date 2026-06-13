#!/usr/bin/env bash
# =============================================================================
# FASE 5 — Passo 2: SignalP6 (predição de peptídeo sinal)
#
# SignalP6 requer licença académica gratuita do DTU.
# Download: https://services.healthtech.dtu.dk/services/SignalP-6.0/
#
# Input:  results/phase4/complete_trypsins.fasta
# Output: results/phase5/signalp/
#   - output.gff3             → posições de clivagem
#   - output.json             → probabilidades por sequência
#   - signalp_summary.tsv    → resumo (ID, tem_sinal, score, posição_clivagem)
# =============================================================================

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config.sh"

INPUT="${RESULTS_DIR}/phase4/complete_trypsins.fasta"
OUT="${RESULTS_DIR}/phase5"
SIGNALP_OUT="${OUT}/signalp"
mkdir -p "$SIGNALP_OUT"
LOG="${LOGS_DIR}/phase5_signalp_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1

log_step "FASE 5 — Passo 2: SignalP6"
check_file "$INPUT" "Tripsinas completas (Fase 4)" || exit 1

# ── Verificar se SignalP6 está instalado ──────────────────────────────────────
SIGNALP_CMD=""

# Procurar em locais comuns
for candidate in signalp6 signalp "$(which signalp 2>/dev/null)" \
    "${HOME}/bin/signalp6" \
    "${HOME}/tools/signalp-6.0/bin/signalp"; do
    if command -v "$candidate" &>/dev/null 2>&1; then
        SIGNALP_CMD="$candidate"
        break
    fi
done

# Tentar via mamba env (se instalado no env analysis)
if [ -z "$SIGNALP_CMD" ]; then
    if mamba run -n ${ENV_ANALYSIS} which signalp6 &>/dev/null 2>&1; then
        SIGNALP_CMD="mamba run --no-capture-output -n ${ENV_ANALYSIS} signalp6"
    elif mamba run -n ${ENV_ANALYSIS} python -c "import signalp" &>/dev/null 2>&1; then
        SIGNALP_CMD="mamba run --no-capture-output -n ${ENV_ANALYSIS} python -m signalp"
    fi
fi

if [ -z "$SIGNALP_CMD" ]; then
    echo "⚠️  SignalP6 não encontrado no sistema."
    echo ""
    echo "Para instalar SignalP6 (licença acadêmica gratuita):"
    echo "  1. Registrar em: https://services.healthtech.dtu.dk/services/SignalP-6.0/"
    echo "  2. Baixar: signalp-6.0.fast.tar.gz"
    echo "  3. Instalar:"
    echo "     tar xf signalp-6.0.fast.tar.gz"
    echo "     mamba activate ${ENV_ANALYSIS}"
    echo "     pip install signalp-6.0.fast/"
    echo "     mamba deactivate"
    echo ""
    echo "⚠️  Criando resultado placeholder para continuar o pipeline..."

    # Criar arquivo placeholder para não bloquear Fase 5
    cat > "${SIGNALP_OUT}/signalp_summary.tsv" << EOF
id	has_signal_peptide	signal_score	cleavage_pos	organism_group
# SignalP6 não executado — instalar manualmente (ver instruções acima)
EOF
    echo "✅ Placeholder criado. Retome após instalar SignalP6."
    exit 0
fi

echo "SignalP6 encontrado: $SIGNALP_CMD"

# ── Executar SignalP6 ─────────────────────────────────────────────────────────
echo "Rodando SignalP6 (modo euk — eucarioto)..."

$SIGNALP_CMD \
    --fastafile "$INPUT" \
    --organism eukarya \
    --output_dir "$SIGNALP_OUT" \
    --format txt \
    --mode fast \
    2>&1

# ── Parsear resultado ─────────────────────────────────────────────────────────
echo ""
echo "Parseando resultados SignalP6..."

python3 - << PYEOF
import sys
from pathlib import Path

signalp_dir = Path("${SIGNALP_OUT}")
summary_out = signalp_dir / "signalp_summary.tsv"

# Tentar ler output no formato txt/gff3/json
results = {}

# Formato txt (SignalP6)
for txt_file in signalp_dir.glob("*.txt"):
    if "summary" in txt_file.name:
        continue
    for line in txt_file.open():
        if line.startswith("#"):
            continue
        parts = line.strip().split("\t")
        if len(parts) >= 2:
            seq_id = parts[0]
            has_signal = parts[1] not in ("OTHER", "0", "False") if len(parts) > 1 else False
            score = parts[2] if len(parts) > 2 else "NA"
            cleavage = parts[4] if len(parts) > 4 else "NA"
            results[seq_id] = {
                "has_signal": has_signal,
                "score": score,
                "cleavage": cleavage,
            }

with open(summary_out, "w") as f:
    f.write("id\thas_signal_peptide\tsignal_score\tcleavage_pos\n")
    for seq_id, data in sorted(results.items()):
        f.write(f"{seq_id}\t{data['has_signal']}\t{data['score']}\t{data['cleavage']}\n")

n_with_signal = sum(1 for d in results.values() if d["has_signal"])
print(f"Total sequências: {len(results)}")
print(f"Com peptídeo sinal: {n_with_signal} ({n_with_signal/len(results)*100:.0f}% se >0)")
PYEOF

echo ""
echo "✅ SignalP6 concluído. Resultado em: ${SIGNALP_OUT}/signalp_summary.tsv"
echo ""
echo "Próximo passo:"
echo "   bash scripts/phase5_primary_char/03_interproscan.sh"
