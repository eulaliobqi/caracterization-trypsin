#!/usr/bin/env bash
# =============================================================================
# FASE 2 — Predição de ORFs
# Ferramenta: TransDecoder (+ hints DIAMOND e HMMER/Pfam opcionais)
#
# Input:  results/phase1/assembly_nr95.fasta
# Output: results/phase2/
#   - agemmatalis.pep            → todas as proteínas preditas
#   - agemmatalis.cds            → CDSs correspondentes
#   - agemmatalis.gff3           → coordenadas no transcriptoma
#   - blast_hints.outfmt6        → hits DIAMOND (se PROTEIN_DB disponível)
#   - pfam_hints.domtblout       → hits Pfam (se Pfam-A.hmm disponível)
#   - phase2_summary.txt         → resumo executivo
# =============================================================================

# ── Ativar ambiente conda ANTES do set -euo pipefail ─────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
_CB="$(conda info --base 2>/dev/null || echo "${HOME}/miniforge3")"
# shellcheck disable=SC1091
source "${_CB}/etc/profile.d/conda.sh"  2>/dev/null || true
source "${_CB}/etc/profile.d/mamba.sh"  2>/dev/null || true
if ! conda activate orf_prediction 2>/dev/null; then
    echo "❌ Ambiente 'orf_prediction' não encontrado."
    echo "   Crie com: mamba env create -f envs/orf_prediction.yml -y"
    exit 1
fi

set -euo pipefail
source "${SCRIPT_DIR}/../config.sh"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --cpus) CPUS="$2"; shift 2 ;;
        --input) INPUT_FASTA="$2"; shift 2 ;;
        *) echo "Argumento desconhecido: $1"; exit 1 ;;
    esac
done

# Input padrão: output da Fase 1
INPUT_FASTA="${INPUT_FASTA:-${RESULTS_DIR}/phase1/assembly_nr95.fasta}"
if [ ! -f "$INPUT_FASTA" ]; then
    echo "⚠️  assembly_nr95.fasta não encontrado — usando assembly original"
    INPUT_FASTA="$ASSEMBLY"
fi

OUT="${RESULTS_DIR}/phase2"
mkdir -p "$OUT"
LOG="${LOGS_DIR}/phase2_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1

log_step "FASE 2: Predição de ORFs (TransDecoder)"
echo "Input:  $INPUT_FASTA"
echo "Output: $OUT"
echo "CPUs:   $CPUS"

check_file "$INPUT_FASTA" "Assembly (input Fase 2)" || exit 1
check_env "$ENV_DISCOVERY" || exit 1

# Mudar para diretório de output (TransDecoder cria arquivos no CWD)
cd "$OUT"

# ═══════════════════════════════════════════════════════════════════════
# PASSO 2A: TransDecoder.LongOrfs — extrair ORFs longas
# ═══════════════════════════════════════════════════════════════════════
log_step "2A: TransDecoder.LongOrfs (min ${MIN_ORF_AA} aa)"

FASTA_BASE="$(basename $INPUT_FASTA)"
# Criar link simbólico se necessário
if [ ! -f "$FASTA_BASE" ]; then
    ln -sf "$INPUT_FASTA" "$FASTA_BASE"
fi

TRANSDECODER_DIR="${OUT}/${FASTA_BASE}.transdecoder_dir"
if [ -d "$TRANSDECODER_DIR" ]; then
    echo "⏩ LongOrfs já executado — pulando."
else
    TransDecoder.LongOrfs \
        -t "$FASTA_BASE" \
        -m "${MIN_ORF_AA}" \
        2>&1
    echo "✅ LongOrfs concluído"
fi

LONGEST_ORFS="${TRANSDECODER_DIR}/longest_orfs.pep"
N_ORFS=$(grep -c "^>" "$LONGEST_ORFS" || echo 0)
echo "ORFs longas identificadas: $N_ORFS"

# ═══════════════════════════════════════════════════════════════════════
# PASSO 2B: BLAST hints (opcional — requer UniProt DIAMOND DB)
# ═══════════════════════════════════════════════════════════════════════
log_step "2B: BLAST hints contra UniProt (opcional)"

BLAST_HINTS="${OUT}/blast_hints.outfmt6"
USE_BLAST_HINTS=""

if [ -f "${UNIPROT_DB}" ]; then
    echo "UniProt DIAMOND DB encontrado: $UNIPROT_DB"
    if [ -f "$BLAST_HINTS" ]; then
        echo "⏩ BLAST hints já existem — pulando."
    else
        echo "Rodando DIAMOND blastp para hints..."
        diamond blastp \
            -q "$LONGEST_ORFS" \
            -d "${UNIPROT_DB}" \
            -p "${CPUS}" \
            -e 1e-5 \
            -k 1 \
            --outfmt 6 \
            -o "$BLAST_HINTS" \
            2>&1
        echo "✅ BLAST hints concluídos"
    fi
    N_BLAST=$(wc -l < "$BLAST_HINTS" || echo 0)
    echo "BLAST hints: $N_BLAST hits"
    USE_BLAST_HINTS="--retain_blastp_hits ${BLAST_HINTS}"
else
    echo "⚠️  UniProt DB não encontrado em ${UNIPROT_DB}"
    echo "   TransDecoder rodará sem hints BLAST (predição menos precisa)"
    echo "   Para baixar: wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
fi

# ═══════════════════════════════════════════════════════════════════════
# PASSO 2C: HMMER/Pfam hints (opcional — requer Pfam-A.hmm)
# ═══════════════════════════════════════════════════════════════════════
log_step "2C: HMMER/Pfam hints (opcional)"

PFAM_HINTS="${OUT}/pfam_hints.domtblout"
USE_PFAM_HINTS=""

if [ -f "${PFAM_HMM}" ] || [ -f "${PFAM_HMM}.h3f" ]; then
    echo "Pfam HMM encontrado: $PFAM_HMM"

    # Pressionar o banco se necessário
    if [ ! -f "${PFAM_HMM}.h3f" ]; then
        echo "Preparando banco Pfam (hmmpress)..."
        hmmpress "${PFAM_HMM}"
    fi

    if [ -f "$PFAM_HINTS" ]; then
        echo "⏩ Pfam hints já existem — pulando."
    else
        echo "Rodando hmmscan contra Pfam-A (pode demorar 1-2h)..."
        hmmscan \
            --cpu "${CPUS}" \
            --domtblout "$PFAM_HINTS" \
            --noali \
            "${PFAM_HMM}" \
            "$LONGEST_ORFS" \
            > "${OUT}/pfam_scan.log" 2>&1
        echo "✅ Pfam hints concluídos"
    fi
    N_PFAM=$(grep -v "^#" "$PFAM_HINTS" | wc -l || echo 0)
    echo "Pfam hits: $N_PFAM"
    USE_PFAM_HINTS="--retain_pfam_hits ${PFAM_HINTS}"
else
    echo "⚠️  Pfam-A.hmm não encontrado em ${PFAM_HMM}"
    echo "   TransDecoder rodará sem hints Pfam"
    echo "   Para baixar: wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
fi

# ═══════════════════════════════════════════════════════════════════════
# PASSO 2D: TransDecoder.Predict
# ═══════════════════════════════════════════════════════════════════════
log_step "2D: TransDecoder.Predict"

PEP_OUT="${OUT}/agemmatalis.pep"

if [ -f "$PEP_OUT" ]; then
    echo "⏩ TransDecoder.Predict já executado — pulando."
else
    TransDecoder.Predict \
        -t "$FASTA_BASE" \
        ${USE_BLAST_HINTS} \
        ${USE_PFAM_HINTS} \
        --single_best_only \
        2>&1

    # Renomear outputs para nomes fixos
    mv "${FASTA_BASE}.transdecoder.pep" "$PEP_OUT"     2>/dev/null || true
    mv "${FASTA_BASE}.transdecoder.cds" "${OUT}/agemmatalis.cds"  2>/dev/null || true
    mv "${FASTA_BASE}.transdecoder.gff3" "${OUT}/agemmatalis.gff3" 2>/dev/null || true
    mv "${FASTA_BASE}.transdecoder.bed"  "${OUT}/agemmatalis.bed"  2>/dev/null || true
    echo "✅ TransDecoder.Predict concluído"
fi

# ═══════════════════════════════════════════════════════════════════════
# RESUMO
# ═══════════════════════════════════════════════════════════════════════
log_step "Resumo da Fase 2"

N_PEP=$(grep -c "^>" "$PEP_OUT" || echo 0)
N_COMPLETE=$(grep -c "complete" "$PEP_OUT" || echo 0)
N_5PRIME=$(grep -c "5prime_partial" "$PEP_OUT" || echo 0)
N_3PRIME=$(grep -c "3prime_partial" "$PEP_OUT" || echo 0)

cat > "${OUT}/phase2_summary.txt" << SUMMARY
FASE 2 — Predição de ORFs
Executado em: $(date)
Input: ${INPUT_FASTA}
ORFs longas (LongOrfs): ${N_ORFS}
Proteínas preditas (Predict): ${N_PEP}
  - Complete (5'+3'): ${N_COMPLETE}
  - 5' partial: ${N_5PRIME}
  - 3' partial: ${N_3PRIME}
Hints BLAST: $([ -f "$BLAST_HINTS" ] && echo "sim (${BLAST_HINTS})" || echo "não usados")
Hints Pfam: $([ -f "$PFAM_HINTS" ] && echo "sim (${PFAM_HINTS})" || echo "não usados")

PRÓXIMO PASSO:
  python scripts/phase2_orf/validate.py
SUMMARY

cat "${OUT}/phase2_summary.txt"
echo ""
echo "✅ Fase 2 concluída. Rode agora:"
echo "   python scripts/phase2_orf/validate.py"
