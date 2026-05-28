#!/usr/bin/env bash
# =============================================================================
# FASE 3 — Passo 1: DIAMOND blastp vs UniProt Swiss-Prot
#
# Busca todas as ORFs preditas contra UniProt para identificar
# candidatos com homologia a serina-proteases (tripsinas).
#
# Input:  results/phase2/agemmatalis.pep
# Output: results/phase3/
#   - diamond_trypsin.tsv        → hits DIAMOND filtrados para tripsinas
#   - diamond_trypsin_ids.txt    → IDs dos candidatos
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
LOG="${LOGS_DIR}/phase3_diamond_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1

log_step "FASE 3 — Passo 1: DIAMOND vs UniProt"
echo "Input:  $PEP"
echo "DB:     $UNIPROT_DB"
echo "CPUs:   $CPUS"

check_file "$PEP" "Proteínas Fase 2" || exit 1
check_env "$ENV_DISCOVERY" || exit 1

if [ ! -f "$UNIPROT_DB" ]; then
    echo "❌ ERRO: UniProt DIAMOND DB não encontrado: $UNIPROT_DB"
    echo ""
    echo "Para criar o banco DIAMOND:"
    echo "  1. Baixar Swiss-Prot:"
    echo "     wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
    echo "     gunzip uniprot_sprot.fasta.gz"
    echo "  2. Criar banco DIAMOND:"
    echo "     mamba run -n ${ENV_DISCOVERY} diamond makedb --in uniprot_sprot.fasta -d ${UNIPROT_DB%.dmnd} --threads ${CPUS}"
    exit 1
fi

# ── DIAMOND blastp ────────────────────────────────────────────────────────────
DIAMOND_RAW="${OUT}/diamond_all.tsv"

if [ -f "$DIAMOND_RAW" ]; then
    echo "⏩ DIAMOND já executado — pulando. Delete ${DIAMOND_RAW} para refazer."
else
    echo "Rodando DIAMOND blastp (pode demorar 30-60min)..."
    ${MAMBA_RUN} ${ENV_DISCOVERY} diamond blastp \
        -q "$PEP" \
        -d "$UNIPROT_DB" \
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

N_RAW=$(wc -l < "$DIAMOND_RAW" || echo 0)
echo "Total de hits DIAMOND: $N_RAW"

# ── Filtrar hits de tripsina ──────────────────────────────────────────────────
DIAMOND_TRYP="${OUT}/diamond_trypsin.tsv"
DIAMOND_IDS="${OUT}/diamond_trypsin_ids.txt"

echo ""
echo "Filtrando hits relacionados a tripsinas/serina-proteases..."

# Palavras-chave para filtrar (case insensitive)
# Inclui: trypsin, serine protease, chymotrypsin, trypsinogen
grep -i -E "trypsin|serine.protease|chymotrypsin|trypsinogen|serine.endopeptidase" \
    "$DIAMOND_RAW" > "$DIAMOND_TRYP" || true

N_TRYP=$(wc -l < "$DIAMOND_TRYP" || echo 0)
echo "Hits de tripsina/serina-protease: $N_TRYP"

# Extrair IDs únicos (coluna 1 do DIAMOND outfmt 6)
cut -f1 "$DIAMOND_TRYP" | sort -u > "$DIAMOND_IDS"
N_IDS=$(wc -l < "$DIAMOND_IDS" || echo 0)
echo "Sequências únicas candidatas (DIAMOND): $N_IDS"

# ── Validação rápida ──────────────────────────────────────────────────────────
if [ "$N_IDS" -lt 10 ]; then
    echo "⚠️  AVISO: Muito poucos candidatos DIAMOND ($N_IDS). Verifique o banco UniProt."
elif [ "$N_IDS" -gt 5000 ]; then
    echo "⚠️  AVISO: Muitos candidatos DIAMOND ($N_IDS). E-value pode estar muito permissivo."
else
    echo "✅ Range razoável de candidatos DIAMOND ($N_IDS)"
fi

echo ""
echo "── Top 10 descrições mais frequentes ──────────────────"
cut -f13 "$DIAMOND_TRYP" | sort | uniq -c | sort -rn | head -10

echo ""
echo "✅ DIAMOND concluído. Próximo passo:"
echo "   bash scripts/phase3_trypsin_id/02_hmmer.sh"
