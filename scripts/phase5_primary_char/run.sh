#!/usr/bin/env bash
# =============================================================================
# FASE 5 — Caracterização Primária das Tripsinas de A. gemmatalis
#
# Passos:
#   5A. CD-HIT proteico (≥90% aa) — colapsar isoformas redundantes
#   5B. Inspeção de sequências longas (>500 aa) — detectar quimeras
#   5C. ProtParam — MW, pI, instabilidade, GRAVY, aliphatic index
#   5D. SignalP-6.0 — predição de peptídeo sinal (graceful fallback)
#   5E. EggNOG-mapper — anotação funcional + GO terms (usa bancos locais)
#
# Input:  results/phase4/complete_trypsins.fasta (67 sequências)
# Output: results/phase5/
#   - trypsins_nr90.fasta         → conjunto não-redundante (≥90% aa)
#   - chimera_suspects.txt        → IDs com >500 aa para verificação
#   - protparam_results.tsv       → propriedades físico-químicas completas
#   - protparam_summary.txt       → tabela para manuscrito
#   - signalp/signalp_summary.tsv → presença/posição do peptídeo sinal
#   - eggnog/eggnog_annotations.tsv → anotação funcional + COG + GO
#   - phase5_summary.tsv          → tabela integrada para publicação
# =============================================================================

# ── Ativar ambiente conda ANTES do set -euo pipefail ─────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
_CB="$(conda info --base 2>/dev/null || echo "${HOME}/miniforge3")"
# shellcheck disable=SC1091
source "${_CB}/etc/profile.d/conda.sh" 2>/dev/null || true
source "${_CB}/etc/profile.d/mamba.sh" 2>/dev/null || true
if ! conda activate annotation 2>/dev/null; then
    echo "❌ Ambiente 'annotation' não encontrado."
    echo "   Crie com: mamba env create -f envs/annotation.yml -y"
    exit 1
fi

set -euo pipefail
source "${SCRIPT_DIR}/../config.sh"

# ── Configuração ──────────────────────────────────────────────────────────────
INPUT="${RESULTS_DIR}/phase4/complete_trypsins.fasta"
INPUT_MUTANT="${RESULTS_DIR}/phase4/mutant_ser195_trypsins.fasta"
OUT="${RESULTS_DIR}/phase5"
mkdir -p "$OUT"
LOG="${LOGS_DIR}/phase5_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1

log_step "FASE 5: Caracterização Primária"
check_file "$INPUT" "complete_trypsins.fasta (Fase 4)" || exit 1

N_INPUT=$(grep -c "^>" "$INPUT" || echo 0)
echo "Input: $N_INPUT tripsinas completas"
echo ""

# ══════════════════════════════════════════════════════════════════════════════
# 5A: CD-HIT — clustering proteico ≥90% identidade
# ══════════════════════════════════════════════════════════════════════════════
log_step "5A: CD-HIT proteico (≥90% aa)"

NR_FASTA="${OUT}/trypsins_nr90.fasta"

if [ -f "$NR_FASTA" ]; then
    echo "⏩ CD-HIT já executado — pulando."
else
    cd-hit \
        -i "$INPUT" \
        -o "$NR_FASTA" \
        -c 0.90 \
        -n 5 \
        -M $((MEM_GB * 1000)) \
        -T "$CPUS" \
        -d 0 \
        2>&1
    echo "✅ CD-HIT proteico concluído"
fi

N_NR=$(grep -c "^>" "$NR_FASTA" || echo 0)
echo "Sequências originais:     $N_INPUT"
echo "Não-redundantes (≥90%):   $N_NR"
echo "Isoformas colapsadas:     $((N_INPUT - N_NR))"
echo ""

# ══════════════════════════════════════════════════════════════════════════════
# 5B: Inspeção de sequências longas (possíveis quimeras)
# ══════════════════════════════════════════════════════════════════════════════
log_step "5B: Detecção de sequências suspeitas (>500 aa)"

CHIMERA_SUSPECTS="${OUT}/chimera_suspects.txt"

python3 - << PYEOF
from Bio import SeqIO
from pathlib import Path

suspects = []
kept = []
fasta = Path("${NR_FASTA}")

for rec in SeqIO.parse(fasta, "fasta"):
    seq_len = len(str(rec.seq).replace("*",""))
    if seq_len > 500:
        suspects.append((rec.id, seq_len))
        print(f"  ⚠️  {rec.id}: {seq_len} aa — SUSPEITA DE QUIMERA")
    else:
        kept.append(rec)

with open("${CHIMERA_SUSPECTS}", "w") as f:
    f.write("id\tlength\tnote\n")
    for sid, slen in suspects:
        f.write(f"{sid}\t{slen}\tverificar_blast_web\n")

# Salvar conjunto filtrado (sem quimeras suspeitas)
nr_filtered = Path("${OUT}/trypsins_nr90_filtered.fasta")
SeqIO.write(kept, str(nr_filtered), "fasta")

print(f"\nSuspeitas de quimera: {len(suspects)}")
print(f"Para análise posterior: {len(kept)} sequências")
print(f"Arquivo filtrado: {nr_filtered}")
PYEOF

echo ""

# Usar conjunto filtrado (sem quimeras) para as análises seguintes
ANALYSIS_FASTA="${OUT}/trypsins_nr90_filtered.fasta"
N_FINAL=$(grep -c "^>" "$ANALYSIS_FASTA" || echo 0)
echo "✅ Conjunto para análise: $N_FINAL sequências"
echo ""

# ══════════════════════════════════════════════════════════════════════════════
# 5C: ProtParam — propriedades físico-químicas
# ══════════════════════════════════════════════════════════════════════════════
log_step "5C: ProtParam (MW, pI, instabilidade, GRAVY)"

PROTPARAM_TSV="${OUT}/protparam_results.tsv"

if [ -f "$PROTPARAM_TSV" ]; then
    echo "⏩ ProtParam já executado — pulando."
else
    python3 - << PYEOF
import sys, re
from pathlib import Path
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# Apenas os 20 aminoácidos padrão — ProteinAnalysis rejeita qualquer outro
STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")

def clean_seq(s):
    """Remove stop codon, gaps, e qualquer aa ambíguo (B, Z, U, J, X, O)."""
    s = s.upper().replace("*","").replace("-","").replace(" ","")
    return re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", s)

fasta = Path("${ANALYSIS_FASTA}")
out_tsv = Path("${PROTPARAM_TSV}")

if not fasta.exists():
    print(f"ERRO: arquivo não encontrado: {fasta}", file=sys.stderr)
    sys.exit(1)

rows = []
n_skipped = 0
for rec in SeqIO.parse(fasta, "fasta"):
    raw = str(rec.seq)
    seq = clean_seq(raw)
    if len(seq) < 10:
        n_skipped += 1
        continue
    try:
        pa = ProteinAnalysis(seq)
        mw  = pa.molecular_weight() / 1000
        pi  = pa.isoelectric_point()
        ii  = pa.instability_index()
        gv  = pa.gravy()
        ai  = pa.aromaticity()
        aa_comp = pa.get_amino_acids_percent()
        aliphatic = (aa_comp.get('A',0) + 2.9 * aa_comp.get('V',0) +
                     3.9 * (aa_comp.get('I',0) + aa_comp.get('L',0))) * 100
        rows.append({
            "id":               rec.id,
            "length_aa":        len(seq),
            "mw_kda":           f"{mw:.2f}",
            "pi":               f"{pi:.2f}",
            "instability_idx":  f"{ii:.2f}",
            "stable":           "yes" if ii < 40 else "no",
            "gravy":            f"{gv:.3f}",
            "aromaticity":      f"{ai:.3f}",
            "aliphatic_idx":    f"{aliphatic:.1f}",
            "mw_ok":            "yes" if 20 <= mw <= 45 else "no",
        })
    except Exception as e:
        print(f"  ⚠️  {rec.id}: erro ProtParam — {e}", file=sys.stderr)
        rows.append({"id": rec.id, "error": str(e)})

print(f"Sequências lidas: {len(rows) + n_skipped} | processadas: {len(rows)} | puladas (<10aa): {n_skipped}")
n_errors = sum(1 for r in rows if "error" in r)
if n_errors:
    print(f"  ⚠️  Com erro: {n_errors} (ver coluna 'error' no TSV)")

with open(out_tsv, "w") as f:
    cols = ["id","length_aa","mw_kda","pi","instability_idx","stable",
            "gravy","aromaticity","aliphatic_idx","mw_ok"]
    f.write("\t".join(cols) + "\n")
    for r in rows:
        if "error" in r:
            f.write(f"{r['id']}\tERROR\t\t\t\t\t\t\t\t{r['error']}\n")
        else:
            f.write("\t".join(str(r.get(c,"")) for c in cols) + "\n")

valid = [r for r in rows if "error" not in r]
if not valid:
    print("❌ ERRO: nenhuma sequência processada com sucesso pelo ProtParam.")
    print("   Verifique se o arquivo FASTA contém sequências proteicas válidas.")
    sys.exit(1)

mws = [float(r["mw_kda"]) for r in valid]
pis = [float(r["pi"])     for r in valid]
n_stable = sum(1 for r in valid if r["stable"] == "yes")
n_mw_ok  = sum(1 for r in valid if r["mw_ok"] == "yes")

print(f"\nSequências analisadas com sucesso: {len(valid)}")
print(f"MW:  {min(mws):.1f}–{max(mws):.1f} kDa  (média: {sum(mws)/len(mws):.1f})")
print(f"pI:  {min(pis):.1f}–{max(pis):.1f}       (média: {sum(pis)/len(pis):.1f})")
print(f"Estáveis (instab < 40):  {n_stable}/{len(valid)}")
print(f"MW esperada (20-45 kDa): {n_mw_ok}/{len(valid)}")
PYEOF
    echo "✅ ProtParam concluído: ${PROTPARAM_TSV}"
fi

# ══════════════════════════════════════════════════════════════════════════════
# 5D: SignalP-6.0 — predição de peptídeo sinal
# ══════════════════════════════════════════════════════════════════════════════
log_step "5D: SignalP-6.0 (peptídeo sinal)"

SIGNALP_DIR="${OUT}/signalp"
mkdir -p "$SIGNALP_DIR"
SIGNALP_SUM="${SIGNALP_DIR}/signalp_summary.tsv"

if [ -f "$SIGNALP_SUM" ]; then
    echo "⏩ SignalP já executado — pulando."
else
    # Procurar SignalP6 instalado
    if command -v signalp6 &>/dev/null 2>/dev/null; then
        echo "Rodando SignalP-6.0..."
        signalp6 \
            --fastafile "$ANALYSIS_FASTA" \
            --organism eukarya \
            --output_dir "$SIGNALP_DIR" \
            --mode fast \
            2>&1
        echo "✅ SignalP-6.0 concluído"
    else
        echo "⚠️  SignalP-6.0 não encontrado."
        echo "   Para instalar:"
        echo "   1. Licença acadêmica gratuita: https://services.healthtech.dtu.dk/services/SignalP-6.0/"
        echo "   2. pip install signalp-6.0h.fast/"
        echo ""
        echo "   Criando placeholder — retome após instalar:"
        cat > "$SIGNALP_SUM" << SPEOF
id	prediction	signal_score	cleavage_pos	note
# SignalP-6.0 não instalado — executar manualmente
SPEOF
        echo "✅ Placeholder criado em ${SIGNALP_SUM}"
    fi
fi

# ══════════════════════════════════════════════════════════════════════════════
# 5E: EggNOG-mapper — anotação funcional
# ══════════════════════════════════════════════════════════════════════════════
log_step "5E: EggNOG-mapper (COG, GO, KEGG)"

EGGNOG_DIR="${OUT}/eggnog"
mkdir -p "$EGGNOG_DIR"
EGGNOG_ANN="${EGGNOG_DIR}/eggnog_annotations.tsv"

if [ -f "${EGGNOG_DIR}/eggnog.emapper.annotations" ]; then
    echo "⏩ EggNOG-mapper já executado — pulando."
    cp "${EGGNOG_DIR}/eggnog.emapper.annotations" "$EGGNOG_ANN" 2>/dev/null || true
else
    if ! command -v emapper.py &>/dev/null 2>/dev/null; then
        echo "⚠️  emapper.py não encontrado."
        echo "   Instale: mamba install -n annotation -c bioconda eggnog-mapper -y"
        echo "   Pulando EggNOG..."
    else
        echo "Rodando EggNOG-mapper (bancos: ${EGGNOG_DB%/*})..."
        emapper.py \
            -i "$ANALYSIS_FASTA" \
            --itype proteins \
            -o eggnog \
            --output_dir "$EGGNOG_DIR" \
            --data_dir "${EGGNOG_DB%/*}" \
            --cpu "$CPUS" \
            --dmnd_db "$EGGNOG_PROTEINS" \
            --override \
            2>&1
        cp "${EGGNOG_DIR}/eggnog.emapper.annotations" "$EGGNOG_ANN" 2>/dev/null || true
        echo "✅ EggNOG-mapper concluído"
    fi
fi

# ══════════════════════════════════════════════════════════════════════════════
# TABELA INTEGRADA PARA PUBLICAÇÃO
# ══════════════════════════════════════════════════════════════════════════════
log_step "Integrando resultados — phase5_summary.tsv"

python3 - << PYEOF
import csv
from pathlib import Path

pp_file  = Path("${PROTPARAM_TSV}")
sp_file  = Path("${SIGNALP_SUM}")
eg_file  = Path("${EGGNOG_ANN}")
out_file = Path("${OUT}/phase5_summary.tsv")

# ProtParam
pp = {}
if pp_file.exists():
    for r in csv.DictReader(open(pp_file), delimiter="\t"):
        pp[r["id"]] = r

# SignalP
sp = {}
if sp_file.exists():
    for r in csv.DictReader(open(sp_file), delimiter="\t"):
        if not r["id"].startswith("#"):
            sp[r["id"]] = r.get("prediction", r.get("has_signal_peptide","NA"))

# EggNOG (colunas variáveis — pegar COG e descrição)
eg = {}
if eg_file.exists():
    for r in csv.DictReader(open(eg_file), delimiter="\t"):
        seq_id = r.get("#query", r.get("query",""))
        eg[seq_id] = {
            "cog":  r.get("COG_category","NA"),
            "desc": r.get("Description", r.get("eggNOG_OGs","NA"))[:60],
            "go":   r.get("GOs","NA")[:40],
        }

all_ids = sorted(set(list(pp.keys()) + list(sp.keys())))

with open(out_file, "w") as f:
    f.write("id\tlength_aa\tmw_kda\tpi\tinstability_idx\tstable\tgravy\t"
            "aliphatic_idx\tsignal_peptide\tcog_category\tdescription\tgo_terms\n")
    for sid in all_ids:
        p = pp.get(sid, {})
        signal = sp.get(sid, "NA")
        e = eg.get(sid, {})
        f.write("\t".join([
            sid,
            p.get("length_aa","NA"), p.get("mw_kda","NA"), p.get("pi","NA"),
            p.get("instability_idx","NA"), p.get("stable","NA"),
            p.get("gravy","NA"), p.get("aliphatic_idx","NA"),
            str(signal),
            e.get("cog","NA"), e.get("desc","NA"), e.get("go","NA"),
        ]) + "\n")

print(f"Tabela integrada: {out_file} ({len(all_ids)} sequências)")
PYEOF

# ══════════════════════════════════════════════════════════════════════════════
# RESUMO FINAL
# ══════════════════════════════════════════════════════════════════════════════
log_step "Resumo da Fase 5"

echo "Input (Fase 4):          $N_INPUT sequências"
echo "Após CD-HIT 90%:         $N_NR sequências"
echo "Para análise (sem >500aa): $N_FINAL sequências"
echo ""
echo "Outputs:"
echo "  results/phase5/trypsins_nr90_filtered.fasta  ← usar para Fase 6 (filogenia)"
echo "  results/phase5/protparam_results.tsv"
echo "  results/phase5/signalp/signalp_summary.tsv"
echo "  results/phase5/eggnog/eggnog_annotations.tsv"
echo "  results/phase5/phase5_summary.tsv  ← tabela integrada"
echo ""
echo "✅ Fase 5 concluída."
echo ""
echo "Próximo passo (Fase 6 — Filogenia):"
echo "   bash scripts/phase6_phylogeny/run.sh   [a criar]"
