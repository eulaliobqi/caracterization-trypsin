#!/usr/bin/env bash
# =============================================================================
# FASE 3 — Identificação de Tripsinas
# Passos: DIAMOND (01) → HMMER PF00089 (02) → Interseção (03)
#
# Input:  results/phase2/agemmatalis.pep
# Output: results/phase3/
#   - trypsins_confident.fasta   → DIAMOND ∩ HMMER (alta confiança)
#   - trypsins_suggestive.fasta  → apenas DIAMOND ou HMMER
#   - identification_report.tsv  → tabela completa
# =============================================================================

# ── Ativar ambiente conda ANTES do set -euo pipefail ─────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
_CB="$(conda info --base 2>/dev/null || echo "${HOME}/miniforge3")"
# shellcheck disable=SC1091
source "${_CB}/etc/profile.d/conda.sh" 2>/dev/null || true
source "${_CB}/etc/profile.d/mamba.sh" 2>/dev/null || true
if ! conda activate orf_prediction 2>/dev/null; then
    echo "❌ Ambiente 'orf_prediction' não encontrado."
    echo "   Crie com: mamba env create -f envs/orf_prediction.yml -y"
    exit 1
fi

set -euo pipefail
source "${SCRIPT_DIR}/../config.sh"

log_step "FASE 3: Identificação de Tripsinas"
echo "Passos: DIAMOND → HMMER PF00089 → Interseção"
echo ""

# Verificar input da Fase 2
PEP="${RESULTS_DIR}/phase2/agemmatalis.pep"
check_file "$PEP" "Proteínas Fase 2 (agemmatalis.pep)" || exit 1
N_PEP=$(grep -c "^>" "$PEP")
echo "Proteínas da Fase 2: $N_PEP"
echo ""

OUT="${RESULTS_DIR}/phase3"
mkdir -p "$OUT"
LOG="${LOGS_DIR}/phase3_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1

# ── PASSO 1: DIAMOND vs NR ────────────────────────────────────────────────────
log_step "3.1: DIAMOND blastp vs NR"

DIAMOND_RAW="${OUT}/diamond_all.tsv"
DIAMOND_TRYP="${OUT}/diamond_trypsin.tsv"
DIAMOND_IDS="${OUT}/diamond_trypsin_ids.txt"

if [ -f "$DIAMOND_IDS" ]; then
    echo "⏩ DIAMOND já executado — pulando."
    N_IDS=$(wc -l < "$DIAMOND_IDS")
    echo "Candidatos DIAMOND: $N_IDS"
else
    if [ ! -f "$DIAMOND_RAW" ]; then
        echo "Rodando DIAMOND blastp (15-30 min)..."
        diamond blastp \
            -q "$PEP" \
            -d "$PROTEIN_DB" \
            -o "$DIAMOND_RAW" \
            -p "$CPUS" \
            -e "$DIAMOND_EVALUE" \
            -k 5 \
            --query-cover 50 \
            --outfmt 6 qseqid sseqid pident length mismatch gapopen \
                         qstart qend sstart send evalue bitscore stitle \
            --sensitive \
            2>&1
        echo "✅ DIAMOND concluído"
    fi

    N_RAW=$(wc -l < "$DIAMOND_RAW")
    echo "Total hits DIAMOND: $N_RAW"

    echo "Filtrando tripsinas/serina-proteases..."
    grep -i -E "trypsin|serine.protease|chymotrypsin|trypsinogen|serine.endopeptidase" \
        "$DIAMOND_RAW" > "$DIAMOND_TRYP" || true

    N_TRYP=$(wc -l < "$DIAMOND_TRYP")
    echo "Hits tripsina/serina-protease: $N_TRYP"

    cut -f1 "$DIAMOND_TRYP" | sort -u > "$DIAMOND_IDS"
    N_IDS=$(wc -l < "$DIAMOND_IDS")
    echo "Sequências únicas candidatas: $N_IDS"
    echo ""
    echo "Top 10 descrições:"
    cut -f13 "$DIAMOND_TRYP" | sort | uniq -c | sort -rn | head -10
fi

# ── PASSO 2: HMMER vs PF00089 ─────────────────────────────────────────────────
log_step "3.2: HMMER hmmsearch vs PF00089 (Tryp_SPc)"

TRYP_SPC_HMM="${OUT}/Tryp_SPc_PF00089.hmm"
HMMER_DOMTBL="${OUT}/hmmer_tryp_spc.domtblout"
HMMER_IDS="${OUT}/hmmer_trypsin_ids.txt"

if [ -f "$HMMER_IDS" ]; then
    echo "⏩ HMMER já executado — pulando."
    N_HMMER=$(wc -l < "$HMMER_IDS")
    echo "Hits PF00089: $N_HMMER"
else
    # Extrair modelo PF00089 isolado (hmmsearch é muito mais rápido que hmmscan)
    if [ ! -f "$TRYP_SPC_HMM" ]; then
        echo "Extraindo modelo PF00089 do Pfam-A..."
        hmmfetch "$PFAM_HMM" "PF00089" > "$TRYP_SPC_HMM" 2>&1 || {
            echo "⚠️  hmmfetch falhou — copiando Pfam-A completo (hmmsearch mais lento)"
            TRYP_SPC_HMM="$PFAM_HMM"
        }
    fi

    if [ ! -f "$HMMER_DOMTBL" ]; then
        echo "Rodando hmmsearch vs PF00089 (rápido ~5-10 min)..."
        hmmsearch \
            --cpu "$CPUS" \
            --domtblout "$HMMER_DOMTBL" \
            -E "$HMMER_EVALUE" \
            --noali \
            "$TRYP_SPC_HMM" \
            "$PEP" \
            > "${OUT}/hmmer_search.log" 2>&1
        echo "✅ hmmsearch concluído"
    fi

    # Filtrar por cobertura do domínio
    echo "Filtrando por cobertura ≥ ${PFAM_COVERAGE}..."
    # Colunas do hmmsearch --domtblout (0-indexado):
    # p[0]  = target name  → ID da proteína (sequência buscada)
    # p[5]  = qlen         → comprimento do perfil HMM
    # p[11] = c-Evalue (domínio)
    # p[15] = hmm coord from
    # p[16] = hmm coord to
    python3 - << PYEOF
from pathlib import Path

domtbl = Path("${HMMER_DOMTBL}")
ids_out = Path("${HMMER_IDS}")
min_cov = float("${PFAM_COVERAGE}")
ev_thr  = float("${HMMER_EVALUE}")
hits = {}
n_total = n_fail = 0

with open(domtbl) as f:
    for line in f:
        if line.startswith("#"): continue
        p = line.split()
        if len(p) < 23: continue
        n_total += 1
        try:
            seq_id   = p[0]          # target = proteína (hmmsearch)
            hmm_len  = int(p[5])     # qlen = comprimento do perfil HMM
            dom_ev   = float(p[11])  # c-Evalue do domínio
            hmm_from = int(p[15])    # início no perfil HMM
            hmm_to   = int(p[16])    # fim no perfil HMM
        except (ValueError, IndexError): continue
        if dom_ev > ev_thr: continue
        cov = (hmm_to - hmm_from + 1) / hmm_len if hmm_len > 0 else 0
        if cov >= min_cov:
            hits[seq_id] = max(hits.get(seq_id, 0), cov)
        else:
            n_fail += 1

print(f"Total hits HMMER: {n_total}")
print(f"Abaixo da cobertura mínima: {n_fail}")
print(f"Aprovados (cobertura ≥ {min_cov}): {len(hits)}")
with open(ids_out, "w") as f:
    for sid in sorted(hits): f.write(sid + "\n")
PYEOF

    N_HMMER=$(wc -l < "$HMMER_IDS")
    echo "IDs com PF00089: $N_HMMER"
fi

# ── PASSO 3: INTERSEÇÃO ───────────────────────────────────────────────────────
log_step "3.3: Interseção DIAMOND ∩ HMMER"

python3 scripts/phase3_trypsin_id/03_intersect.py \
    --results3 "$OUT" \
    --pep "$PEP"

echo ""
echo "✅ Fase 3 concluída. Rode agora:"
echo "   python scripts/phase3_trypsin_id/validate.py"
