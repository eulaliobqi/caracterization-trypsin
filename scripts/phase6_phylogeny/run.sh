#!/usr/bin/env bash
# =============================================================================
# FASE 6 — Análise Filogenética das Tripsinas de A. gemmatalis
#
# Passos:
#   6A. Download de referências (NCBI: A. gemmatalis + Lepidoptera + outgroup)
#   6B. Renomear e combinar sequências (48 nossas + referências)
#   6C. Alinhamento múltiplo — MAFFT-linsi
#   6D. Trimagem — trimAl (automated1)
#   6E. Inferência filogenética — IQ-TREE2 (modelo auto + UFBoot 1000)
#   6F. Resumo e estatísticas da árvore
#
# Input:  results/phase5/trypsins_nr90_filtered.fasta (48 sequências)
# Output: results/phase6/
#   - refs_lepidoptera.fasta       → referências baixadas e curadas
#   - combined_for_tree.fasta      → nossas seqs + referências (renomeadas)
#   - alignment_mafft.fasta        → alinhamento completo
#   - alignment_trimmed.fasta      → alinhamento após trimAl
#   - trypsin_tree.treefile        → árvore Newick (IQ-TREE2)
#   - trypsin_tree.iqtree          → relatório IQ-TREE2 (modelo, bootstrap)
# =============================================================================

# ── Ativar ambiente conda ANTES do set -euo pipefail ─────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
_CB="$(conda info --base 2>/dev/null || echo "${HOME}/miniforge3")"
# shellcheck disable=SC1091
source "${_CB}/etc/profile.d/conda.sh" 2>/dev/null || true
source "${_CB}/etc/profile.d/mamba.sh" 2>/dev/null || true
if ! conda activate phylogeny 2>/dev/null; then
    echo "❌ Ambiente 'phylogeny' não encontrado."
    echo "   Crie com: mamba env create -f envs/phylogeny.yml -y"
    exit 1
fi

set -euo pipefail
source "${SCRIPT_DIR}/../config.sh"

# ── Configuração ──────────────────────────────────────────────────────────────
INPUT="${RESULTS_DIR}/phase5/trypsins_nr90_filtered.fasta"
OUT="${RESULTS_DIR}/phase6"
REFS_DIR_LOCAL="${OUT}/refs"
mkdir -p "$OUT" "$REFS_DIR_LOCAL"
LOG="${LOGS_DIR}/phase6_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1

log_step "FASE 6: Análise Filogenética"
check_file "$INPUT" "trypsins_nr90_filtered.fasta (Fase 5)" || exit 1

N_INPUT=$(grep -c "^>" "$INPUT" || echo 0)
echo "Sequências de A. gemmatalis: $N_INPUT"
echo ""

# ══════════════════════════════════════════════════════════════════════════════
# 6A: Download de sequências de referência
# ══════════════════════════════════════════════════════════════════════════════
log_step "6A: Download de referências (NCBI + UniProt)"

REFS_FASTA="${REFS_DIR_LOCAL}/refs_lepidoptera.fasta"

if [ -f "$REFS_FASTA" ] && [ "$(grep -c "^>" "$REFS_FASTA" || echo 0)" -gt 5 ]; then
    N_REFS=$(grep -c "^>" "$REFS_FASTA" || echo 0)
    echo "⏩ Referências já baixadas ($N_REFS seqs) — pulando."
else
    echo "Baixando referências da NCBI e UniProt..."
    > "$REFS_FASTA"  # limpar arquivo

    # ── Função de download segura ──────────────────────────────────────────
    download_ncbi() {
        local acc="$1"
        local label="$2"
        local tmp
        tmp=$(mktemp)
        if curl -sf --max-time 30 --retry 3 \
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=${acc}&rettype=fasta&retmode=text" \
            -o "$tmp" && grep -q "^>" "$tmp"; then
            # Substituir cabeçalho pelo label limpo
            python3 -c "
import sys
lines = open('$tmp').readlines()
first = True
for l in lines:
    if l.startswith('>') and first:
        print(f'>${label}')
        first = False
    else:
        print(l, end='')
" >> "$REFS_FASTA"
            echo "  ✅ $acc → $label"
        else
            echo "  ⚠️  $acc não baixado (ignorado)"
        fi
        rm -f "$tmp"
        sleep 0.35  # respeitar limite NCBI (3 req/s)
    }

    download_uniprot() {
        local acc="$1"
        local label="$2"
        local tmp
        tmp=$(mktemp)
        if curl -sf --max-time 30 --retry 3 \
            "https://rest.uniprot.org/uniprotkb/${acc}.fasta" \
            -o "$tmp" && grep -q "^>" "$tmp"; then
            python3 -c "
import sys
lines = open('$tmp').readlines()
first = True
for l in lines:
    if l.startswith('>') and first:
        print(f'>${label}')
        first = False
    else:
        print(l, end='')
" >> "$REFS_FASTA"
            echo "  ✅ $acc → $label"
        else
            echo "  ⚠️  $acc não baixado (ignorado)"
        fi
        rm -f "$tmp"
        sleep 0.35
    }

    # ── A. gemmatalis — sequências publicadas e anotadas ──────────────────
    echo ""
    echo "── Anticarsia gemmatalis (literatura) ──────────────"
    download_ncbi "AWL83213.1"  "Agem_trypsin1_AWL83213"
    download_ncbi "AWL83214.1"  "Agem_trypsin2_AWL83214"
    download_ncbi "AAK15127.1"  "Agem_tryp_Brito2001_AAK15127"
    download_ncbi "AAK15128.1"  "Agem_tryp_Brito2001_AAK15128"
    download_ncbi "AAN71630.1"  "Agem_serineP1_AAN71630"
    download_ncbi "AAN71631.1"  "Agem_serineP2_AAN71631"
    download_ncbi "AAN71632.1"  "Agem_serpin1a_AAN71632"

    # Busca dinâmica: A. gemmatalis trypsins no NCBI (até 20 resultados)
    echo ""
    echo "  Buscando sequências adicionais de A. gemmatalis no NCBI..."
    AGEM_IDS=$(curl -sf --max-time 20 \
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=protein&term=Anticarsia+gemmatalis%5BOrganism%5D+AND+trypsin%5BTitle%5D&retmax=30&retmode=json" \
        | python3 -c "
import sys,json
try:
    data = json.load(sys.stdin)
    ids = data['esearchresult']['idlist']
    print(' '.join(ids))
except: pass
" 2>/dev/null || echo "")

    if [ -n "$AGEM_IDS" ]; then
        echo "  IDs encontrados: $AGEM_IDS"
        # Baixar em lote (até 30)
        IDS_COMMA=$(echo "$AGEM_IDS" | tr ' ' ',')
        tmp_batch=$(mktemp)
        if curl -sf --max-time 60 --retry 3 \
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=${IDS_COMMA}&rettype=fasta&retmode=text" \
            -o "$tmp_batch" && grep -q "^>" "$tmp_batch"; then
            python3 - "$tmp_batch" "$REFS_FASTA" << 'PYEOF'
import sys, re
from Bio import SeqIO
batch_file, out_file = sys.argv[1], sys.argv[2]

# Ler existentes para evitar duplicata
existing = set()
try:
    for rec in SeqIO.parse(out_file, "fasta"):
        existing.add(rec.id)
except: pass

added = 0
with open(out_file, "a") as fout:
    for rec in SeqIO.parse(batch_file, "fasta"):
        if rec.id not in existing:
            # Criar label limpo
            desc = rec.description[:60]
            label = re.sub(r'\s+', '_', f"Agem_search_{rec.id}").replace("|","_")
            fout.write(f">{label}\n{str(rec.seq)}\n")
            existing.add(rec.id)
            added += 1
print(f"  {added} sequências adicionais de A. gemmatalis adicionadas")
PYEOF
        fi
        rm -f "$tmp_batch"
    fi
    sleep 0.5

    # ── Lepidoptera: Noctuidae / Erebidae (família próxima) ───────────────
    echo ""
    echo "── Noctuidae / Erebidae (família próxima) ──────────"
    download_ncbi "AAN71635.1"  "Mconf_serpinP1_AAN71635"   # Mamestra configurata
    download_ncbi "AAN71636.1"  "Mconf_serpinP2_AAN71636"
    download_ncbi "AFM28256.1"  "Hvir_chymotrypsin_AFM28256"  # Heliothis virescens
    download_ncbi "ABR88239.1"  "Hvir_chymotrypsinC9_ABR88239"
    download_ncbi "AAF71515.1"  "Aips_chymotrypsinogen_AAF71515"  # Agrotis ipsilon

    echo ""
    echo "── Spodoptera frugiperda ────────────────────────────"
    download_ncbi "AAN71637.1"  "Sfru_serineP1_AAN71637"
    download_ncbi "AAN71638.1"  "Sfru_serineP2_AAN71638"
    download_uniprot "Q9GU24"   "Sfru_trypsin_Q9GU24"

    echo ""
    echo "── Helicoverpa armigera ─────────────────────────────"
    download_ncbi "ABC96738.1"  "Harm_trypsin_ABC96738"
    download_ncbi "ABS30437.1"  "Harm_trypsin2_ABS30437"
    download_uniprot "Q9BLK7"   "Harm_trypsin_Q9BLK7"

    echo ""
    echo "── Trichoplusia ni ──────────────────────────────────"
    download_ncbi "XP_026745783.1"  "Tni_trypsin_alkC_XP026745783"
    download_ncbi "XP_026739093.1"  "Tni_trypsin_alkC2_XP026739093"
    download_ncbi "XP_026746126.1"  "Tni_serineP1_XP026746126"

    echo ""
    echo "── Manduca sexta (Sphingidae) ───────────────────────"
    download_ncbi "AAA29523.1"  "Msex_trypsin1_AAA29523"
    download_uniprot "P35035"   "Msex_trypsin_P35035"

    echo ""
    echo "── Bombyx mori (Bombycidae) ─────────────────────────"
    download_uniprot "Q9Y0B3"   "Bmor_trypsin_Q9Y0B3"
    download_ncbi "XP_004931024.1"  "Bmor_trypsin_XP004931024"

    echo ""
    echo "── Leguminivora glycinivorella (Tortricidae) ────────"
    download_ncbi "XP_047989729.1"  "Lgly_trypsinCFT1_XP047989729"

    echo ""
    echo "── Outgroup: Bos taurus trypsin (vertebrado) ───────"
    download_uniprot "P00760"  "Btau_trypsin_P00760_OUTGROUP"

    N_REFS=$(grep -c "^>" "$REFS_FASTA" || echo 0)
    echo ""
    echo "✅ Total de referências baixadas: $N_REFS sequências"
fi

# ══════════════════════════════════════════════════════════════════════════════
# 6B: Renomear e combinar sequências
# ══════════════════════════════════════════════════════════════════════════════
log_step "6B: Combinando sequências (A. gemmatalis + referências)"

COMBINED="${OUT}/combined_for_tree.fasta"

python3 - << PYEOF
import re
from Bio import SeqIO
from pathlib import Path

our_fasta  = Path("${INPUT}")
refs_fasta = Path("${REFS_FASTA}")
combined   = Path("${COMBINED}")

records = []
seen_ids = set()

# Nossas sequências — prefixo "Agem_NEW_" para distinguir das publicadas
for rec in SeqIO.parse(our_fasta, "fasta"):
    # Abreviar: TRINITY_DN1205_c0_g1_i6.p1 → Agem_NEW_DN1205_g1_i6
    short = re.sub(r'TRINITY_', '', rec.id)
    short = re.sub(r'_c\d+', '', short)
    short = re.sub(r'\.p\d+', '', short)
    label = f"Agem_NEW_{short}"[:50]
    if label not in seen_ids:
        seq_clean = str(rec.seq).replace("*","")
        records.append((label, seq_clean))
        seen_ids.add(label)

n_ours = len(records)

# Referências
for rec in SeqIO.parse(refs_fasta, "fasta"):
    label = rec.id[:50]
    if label not in seen_ids:
        seq_clean = str(rec.seq).replace("*","")
        records.append((label, seq_clean))
        seen_ids.add(label)

n_refs = len(records) - n_ours

with open(combined, "w") as f:
    for label, seq in records:
        f.write(f">{label}\n{seq}\n")

print(f"Nossas sequências: {n_ours}")
print(f"Referências:       {n_refs}")
print(f"Total combinado:   {len(records)}")
PYEOF

N_COMBINED=$(grep -c "^>" "$COMBINED" || echo 0)
echo "✅ $N_COMBINED sequências para alinhamento"

# ══════════════════════════════════════════════════════════════════════════════
# 6C: Alinhamento múltiplo — MAFFT-linsi
# ══════════════════════════════════════════════════════════════════════════════
log_step "6C: MAFFT-linsi (alinhamento máxima precisão)"

ALIGNED="${OUT}/alignment_mafft.fasta"

if [ -f "$ALIGNED" ] && [ "$(grep -c "^>" "$ALIGNED" || echo 0)" -eq "$N_COMBINED" ]; then
    echo "⏩ Alinhamento MAFFT já existe — pulando."
else
    echo "Rodando MAFFT-linsi (pode demorar 5-20 min para ~80 sequências)..."
    mafft \
        --localpair \
        --maxiterate 1000 \
        --reorder \
        --thread "$CPUS" \
        "$COMBINED" \
        > "$ALIGNED" \
        2> "${LOGS_DIR}/mafft_phase6.log"
    echo "✅ Alinhamento concluído"
fi

N_ALIGNED=$(grep -c "^>" "$ALIGNED" || echo 0)
# Comprimento do alinhamento (primeira sequência)
ALN_LEN=$(python3 -c "
from Bio import SeqIO
for r in SeqIO.parse('${ALIGNED}', 'fasta'):
    print(len(r.seq)); break
" 2>/dev/null || echo "?")
echo "Sequências alinhadas: $N_ALIGNED  |  Comprimento: ${ALN_LEN} posições"

# ══════════════════════════════════════════════════════════════════════════════
# 6D: Trimagem — trimAl
# ══════════════════════════════════════════════════════════════════════════════
log_step "6D: trimAl (remoção de colunas mal alinhadas)"

TRIMMED="${OUT}/alignment_trimmed.fasta"

if [ -f "$TRIMMED" ]; then
    echo "⏩ trimAl já executado — pulando."
else
    # -gappyout: remove apenas colunas com gap em maioria — menos agressivo que -automated1
    # Para sequências divergentes inter-espécies, -gappyout retém ~30-60% das posições
    trimal \
        -in  "$ALIGNED" \
        -out "$TRIMMED" \
        -gappyout \
        -fasta \
        2>&1

    # Verificar se reteve posições suficientes (≥15% do original)
    TRIM_CHECK=$(python3 -c "
from Bio import SeqIO
n = sum(1 for _ in open('${TRIMMED}') if _.startswith('>'))
for r in SeqIO.parse('${TRIMMED}', 'fasta'):
    print(len(r.seq)); break
" 2>/dev/null || echo "0")
    PCT_CHECK=$(python3 -c "
try:
    pct = int('${TRIM_CHECK}') / int('${ALN_LEN}') * 100
    print(f'{pct:.1f}')
except: print('0')
" 2>/dev/null || echo "0")

    if python3 -c "exit(0 if float('${PCT_CHECK}') >= 15 else 1)" 2>/dev/null; then
        echo "✅ trimAl concluído (${PCT_CHECK}% retido)"
    else
        echo "⚠️  trimAl reteve muito pouco (${PCT_CHECK}%) — usando -gt 0.5 em vez de -gappyout"
        trimal \
            -in  "$ALIGNED" \
            -out "$TRIMMED" \
            -gt  0.5 \
            -fasta \
            2>&1
        echo "✅ trimAl (gt=0.5) concluído"
    fi
fi

TRIM_LEN=$(python3 -c "
from Bio import SeqIO
for r in SeqIO.parse('${TRIMMED}', 'fasta'):
    print(len(r.seq)); break
" 2>/dev/null || echo "?")
echo "Posições após trimagem: ${TRIM_LEN}  (original: ${ALN_LEN})"
PCT_RETAINED=$(python3 -c "
try:
    pct = int('${TRIM_LEN}') / int('${ALN_LEN}') * 100
    print(f'{pct:.1f}')
except: print('?')
" 2>/dev/null || echo "?")
echo "Retidas: ${PCT_RETAINED}% das posições"

# ══════════════════════════════════════════════════════════════════════════════
# 6E: Inferência filogenética — IQ-TREE2
# ══════════════════════════════════════════════════════════════════════════════
log_step "6E: IQ-TREE2 (modelo automático + UFBoot 1000)"

TREE_PREFIX="${OUT}/trypsin_tree"

if [ -f "${TREE_PREFIX}.treefile" ]; then
    echo "⏩ IQ-TREE2 já executado — pulando."
else
    # Detectar binário: iqtree2 (instalação padrão) ou iqtree (alguns conda)
    IQTREE_BIN=""
    for bin in iqtree2 iqtree; do
        if command -v "$bin" &>/dev/null; then
            IQTREE_BIN="$bin"
            echo "IQ-TREE encontrado: $($bin --version 2>&1 | head -1)"
            break
        fi
    done
    if [ -z "$IQTREE_BIN" ]; then
        echo "❌ IQ-TREE não encontrado. Instale: mamba install -n phylogeny -c bioconda iqtree -y"
        exit 1
    fi

    echo "Rodando $IQTREE_BIN (modelo TEST + UFBoot 1000 + SH-aLRT 1000)..."
    "$IQTREE_BIN" \
        -s    "$TRIMMED" \
        --prefix "$TREE_PREFIX" \
        -m    TEST \
        -B    1000 \
        -alrt 1000 \
        -T    "$CPUS" \
        --redo \
        2>&1 | tee "${LOGS_DIR}/iqtree2_phase6.log"
    echo "✅ IQ-TREE concluído"
fi

# ══════════════════════════════════════════════════════════════════════════════
# 6F: Resumo da análise filogenética
# ══════════════════════════════════════════════════════════════════════════════
log_step "Resumo da Fase 6"

# Extrair modelo selecionado e score
BEST_MODEL=$(grep "^Best-fit model:" "${TREE_PREFIX}.iqtree" 2>/dev/null | head -1 || echo "ver ${TREE_PREFIX}.iqtree")
LNLIKE=$(grep "^Log-likelihood" "${TREE_PREFIX}.iqtree" 2>/dev/null | head -1 || echo "ver arquivo")

echo "Sequências na árvore:  $N_ALIGNED"
echo "Posições no alinhamento (após trimAl): ${TRIM_LEN}"
echo "Modelo selecionado:    $BEST_MODEL"
echo "Log-likelihood:        $LNLIKE"
echo ""

# Verificar bootstrap: quantas ramos com UFBoot ≥ 95?
N_HIGH_BS=$(python3 -c "
import re
try:
    tree = open('${TREE_PREFIX}.treefile').read()
    vals = re.findall(r'/(\d+)', tree)  # SH-aLRT/UFBoot
    ufs = [int(v) for v in vals if v.isdigit()]
    high = sum(1 for v in ufs if v >= 95)
    total = len(ufs)
    print(f'{high}/{total} ramos com UFBoot ≥ 95')
except Exception as e:
    print(f'(não calculado: {e})')
" 2>/dev/null || echo "(não calculado)")

echo "Suporte bootstrap:     $N_HIGH_BS"
echo ""
echo "── Arquivos gerados ────────────────────────────────"
echo "  ${REFS_DIR_LOCAL}/refs_lepidoptera.fasta"
echo "  ${COMBINED}"
echo "  ${ALIGNED}"
echo "  ${TRIMMED}"
echo "  ${TREE_PREFIX}.treefile   ← árvore Newick para visualização"
echo "  ${TREE_PREFIX}.iqtree     ← relatório completo IQ-TREE2"
echo ""
echo "✅ Fase 6 concluída."
echo ""
echo "Para visualizar a árvore:"
echo "  - FigTree (GUI): abrir ${TREE_PREFIX}.treefile"
echo "  - iTOL (web): upload ${TREE_PREFIX}.treefile em https://itol.embl.de/"
echo ""
echo "Próximo passo (Fase 7 — AlphaFold3):"
echo "   bash scripts/phase7_structure/run.sh   [requer GPU]"
