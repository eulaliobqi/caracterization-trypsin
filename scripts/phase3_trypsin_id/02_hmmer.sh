#!/usr/bin/env bash
# =============================================================================
# FASE 3 — Passo 2: HMMER vs Pfam PF00089 (Tryp_SPc)
#
# Busca o domínio Tryp_SPc (PF00089) em todas as ORFs.
# PF00089 é o domínio catalítico de tripsinas/serina-proteases S1.
#
# Input:  results/phase2/agemmatalis.pep
# Output: results/phase3/
#   - hmmer_tryp_spc.domtblout   → hits detalhados por domínio
#   - hmmer_trypsin_ids.txt      → IDs com domínio PF00089 (cobertura ≥ 80%)
# =============================================================================

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config.sh"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --cpus) CPUS="$2"; shift 2 ;;
        *) echo "Argumento desconhecido: $1"; exit 1 ;;
    esac
done

PEP="${RESULTS_DIR}/phase2/agemmatalis.pep"
OUT="${RESULTS_DIR}/phase3"
mkdir -p "$OUT"
LOG="${LOGS_DIR}/phase3_hmmer_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1

log_step "FASE 3 — Passo 2: HMMER vs PF00089 (Tryp_SPc)"
echo "Input:  $PEP"
echo "Pfam:   $PFAM_HMM"
echo "CPUs:   $CPUS"
echo "e-value: $HMMER_EVALUE"
echo "cobertura mínima: ${PFAM_COVERAGE}"

check_file "$PEP" "Proteínas Fase 2" || exit 1
check_env "$ENV_DISCOVERY" || exit 1

if [ ! -f "$PFAM_HMM" ] && [ ! -f "${PFAM_HMM}.h3f" ]; then
    echo "❌ ERRO: Pfam-A.hmm não encontrado: $PFAM_HMM"
    echo ""
    echo "Para baixar e preparar o banco Pfam:"
    echo "  mkdir -p ${REFS_DIR}"
    echo "  wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz -P ${REFS_DIR}"
    echo "  gunzip ${REFS_DIR}/Pfam-A.hmm.gz"
    echo "  mamba run -n ${ENV_DISCOVERY} hmmpress ${PFAM_HMM}"
    echo ""
    echo "⚠️  ALTERNATIVA: usar apenas DIAMOND (sem HMMER)"
    echo "   Neste caso, pule para: python scripts/phase3_trypsin_id/03_intersect.py --diamond-only"
    exit 1
fi

# Pressionar banco se necessário
if [ ! -f "${PFAM_HMM}.h3f" ]; then
    echo "Preparando banco Pfam (hmmpress)..."
    ${MAMBA_RUN} ${ENV_DISCOVERY} hmmpress "${PFAM_HMM}"
fi

# ── Extrair apenas o modelo PF00089 do Pfam-A.hmm ────────────────────────────
# Usar hmmsearch com o modelo específico é mais eficiente que hmmscan com Pfam inteiro
# Mas se não tiver o modelo isolado, usar hmmscan com grep

TRYP_SPC_HMM="${OUT}/Tryp_SPc_PF00089.hmm"
HMMER_DOMTBL="${OUT}/hmmer_tryp_spc.domtblout"
HMMER_IDS="${OUT}/hmmer_trypsin_ids.txt"

# Extrair modelo PF00089 do Pfam-A.hmm
if [ ! -f "$TRYP_SPC_HMM" ]; then
    echo "Extraindo modelo PF00089 (Tryp_SPc) do Pfam-A..."
    ${MAMBA_RUN} ${ENV_DISCOVERY} hmmfetch "${PFAM_HMM}" "PF00089" > "$TRYP_SPC_HMM" 2>&1 || {
        echo "⚠️  hmmfetch falhou — usando hmmscan com Pfam completo (mais lento)"
        TRYP_SPC_HMM="$PFAM_HMM"
    }
fi

# ── HMMER search ──────────────────────────────────────────────────────────────
if [ -f "$HMMER_DOMTBL" ]; then
    echo "⏩ HMMER já executado — pulando. Delete ${HMMER_DOMTBL} para refazer."
else
    if [ "$TRYP_SPC_HMM" = "$PFAM_HMM" ]; then
        echo "Rodando hmmscan (Pfam completo) — pode demorar 1-2h..."
        ${MAMBA_RUN} ${ENV_DISCOVERY} hmmscan \
            --cpu "$CPUS" \
            --domtblout "$HMMER_DOMTBL" \
            -E "$HMMER_EVALUE" \
            --noali \
            "$PFAM_HMM" \
            "$PEP" \
            > "${OUT}/hmmer_scan.log" 2>&1
        # Filtrar apenas PF00089
        grep -E "^PF00089|^#" "$HMMER_DOMTBL" > "${OUT}/hmmer_tryp_spc_filtered.domtblout"
        mv "${OUT}/hmmer_tryp_spc_filtered.domtblout" "$HMMER_DOMTBL"
    else
        echo "Rodando hmmsearch com modelo PF00089 (rápido)..."
        ${MAMBA_RUN} ${ENV_DISCOVERY} hmmsearch \
            --cpu "$CPUS" \
            --domtblout "$HMMER_DOMTBL" \
            -E "$HMMER_EVALUE" \
            --noali \
            "$TRYP_SPC_HMM" \
            "$PEP" \
            > "${OUT}/hmmer_search.log" 2>&1
    fi
    echo "✅ HMMER concluído"
fi

# ── Filtrar por cobertura do domínio ─────────────────────────────────────────
echo ""
echo "Filtrando hits por cobertura do domínio (≥ ${PFAM_COVERAGE})..."

python3 - << PYEOF
import sys
from pathlib import Path

domtbl = Path("${HMMER_DOMTBL}")
ids_out = Path("${HMMER_IDS}")
min_coverage = float("${PFAM_COVERAGE}")
evalue_thresh = float("${HMMER_EVALUE}")

hits = {}
n_total = 0
n_coverage_fail = 0

with open(domtbl) as f:
    for line in f:
        if line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 23:
            continue
        n_total += 1

        # Colunas do domtblout:
        # 0=target_name, 2=query_name(seq), 5=hmm_len,
        # 7=evalue(domain), 17=hmm_from, 18=hmm_to
        try:
            target_name = parts[0]   # modelo HMM (PF00089)
            seq_id      = parts[2]   # sequência query (nossa proteína)
            hmm_len     = int(parts[5])
            dom_evalue  = float(parts[11])
            hmm_from    = int(parts[17])
            hmm_to      = int(parts[18])
        except (ValueError, IndexError):
            continue

        if dom_evalue > evalue_thresh:
            continue

        # Cobertura: quantos % do perfil HMM foram alinhados?
        coverage = (hmm_to - hmm_from + 1) / hmm_len if hmm_len > 0 else 0

        if coverage >= min_coverage:
            hits[seq_id] = max(hits.get(seq_id, 0), coverage)
        else:
            n_coverage_fail += 1

print(f"Total hits HMMER: {n_total}")
print(f"Filtrados (cobertura < {min_coverage}): {n_coverage_fail}")
print(f"Hits com cobertura ≥ {min_coverage}: {len(hits)}")

with open(ids_out, "w") as f:
    for seq_id in sorted(hits.keys()):
        f.write(seq_id + "\n")

print(f"IDs salvos em: {ids_out}")
PYEOF

N_HMMER=$(wc -l < "$HMMER_IDS" || echo 0)
echo ""
echo "Sequências com domínio PF00089 (Tryp_SPc): $N_HMMER"

if [ "$N_HMMER" -lt 5 ]; then
    echo "⚠️  AVISO: Muito poucos hits HMMER ($N_HMMER). Verifique o Pfam-A.hmm."
elif [ "$N_HMMER" -gt 1000 ]; then
    echo "⚠️  AVISO: Muitos hits HMMER ($N_HMMER). Pode indicar E-value muito permissivo."
else
    echo "✅ Range razoável de hits HMMER ($N_HMMER)"
fi

echo ""
echo "✅ HMMER concluído. Próximo passo:"
echo "   python scripts/phase3_trypsin_id/03_intersect.py"
