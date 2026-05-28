#!/usr/bin/env bash
# =============================================================================
# FASE 5 — Passo 3: InterProScan (anotação de domínios e famílias)
#
# InterProScan integra múltiplos bancos: Pfam, PANTHER, SMART, etc.
# Confirma: domínio PF00089 (Tryp_SPc), família S01, peptídeo sinal
#
# Download: https://interproscan.org/
# Tamanho: ~12 GB
#
# Alternativa online (sem instalação local):
#   https://www.ebi.ac.uk/interpro/search/sequence/
#
# Input:  results/phase4/complete_trypsins.fasta
# Output: results/phase5/interproscan/
# =============================================================================

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config.sh"

INPUT="${RESULTS_DIR}/phase4/complete_trypsins.fasta"
OUT="${RESULTS_DIR}/phase5"
IPS_OUT="${OUT}/interproscan"
mkdir -p "$IPS_OUT"
LOG="${LOGS_DIR}/phase5_interproscan_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1

log_step "FASE 5 — Passo 3: InterProScan"
check_file "$INPUT" "Tripsinas completas (Fase 4)" || exit 1

# ── Verificar InterProScan ────────────────────────────────────────────────────
IPS_CMD=""
for candidate in interproscan interproscan.sh \
    "${HOME}/tools/interproscan/interproscan.sh" \
    "/opt/interproscan/interproscan.sh"; do
    if [ -f "$candidate" ] || command -v "$candidate" &>/dev/null 2>&1; then
        IPS_CMD="$candidate"
        break
    fi
done

if [ -z "$IPS_CMD" ]; then
    echo "⚠️  InterProScan não encontrado."
    echo ""
    echo "Para instalar InterProScan:"
    echo "  # Requer Java 11+ e ~12 GB de espaço"
    echo "  mkdir -p ~/tools && cd ~/tools"
    echo "  wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.68-100.0/interproscan-5.68-100.0-64-bit.tar.gz"
    echo "  tar xf interproscan-5.68-100.0-64-bit.tar.gz"
    echo "  IPS_CMD=~/tools/interproscan-5.68-100.0/interproscan.sh"
    echo ""
    echo "Alternativa online (para poucas sequências):"
    echo "  1. Ir para: https://www.ebi.ac.uk/interpro/search/sequence/"
    echo "  2. Colar FASTA de results/phase4/complete_trypsins.fasta"
    echo "  3. Baixar resultado como TSV"
    echo "  4. Salvar em: results/phase5/interproscan/interproscan.tsv"
    echo ""
    echo "⚠️  Criando placeholder para não bloquear pipeline..."

    # Placeholder
    echo -e "protein_accession\tmd5\tlength\tanalysis\tsignature_accession\tsignature_description\tstart\tstop\tscore\tstatus\tdate\tinterpro_accession\tinterpro_description" \
        > "${IPS_OUT}/interproscan.tsv"
    echo "# InterProScan não executado — instalar manualmente ou usar versão online" \
        >> "${IPS_OUT}/interproscan.tsv"
    exit 0
fi

echo "InterProScan encontrado: $IPS_CMD"

# ── Remover asteriscos das sequências (InterProScan não aceita stop codons) ───
CLEAN_INPUT="${IPS_OUT}/complete_trypsins_clean.fasta"
sed 's/\*//g' "$INPUT" > "$CLEAN_INPUT"

# ── Executar InterProScan ────────────────────────────────────────────────────
echo "Rodando InterProScan (pode demorar 30-120 min)..."

"$IPS_CMD" \
    -i "$CLEAN_INPUT" \
    -o "${IPS_OUT}/interproscan.tsv" \
    -f TSV \
    -appl Pfam,PANTHER,SMART,ProSiteProfiles \
    -goterms \
    -pa \
    --cpu "$CPUS" \
    -dp \
    2>&1

echo "✅ InterProScan concluído"

# ── Parsear e resumir ─────────────────────────────────────────────────────────
echo ""
echo "Resumindo resultados InterProScan..."

python3 - << PYEOF
from pathlib import Path
from collections import defaultdict, Counter

ips_tsv = Path("${IPS_OUT}/interproscan.tsv")
if not ips_tsv.exists():
    print("❌ interproscan.tsv não encontrado")
    exit(1)

# Colunas: protein_acc, md5, length, analysis, sig_acc, sig_desc, start, stop, score, status, date, ipr_acc, ipr_desc
seq_domains = defaultdict(set)
domain_counts = Counter()
pfam_present = {}
tryp_spc_count = 0

for line in ips_tsv.open():
    if line.startswith("#"):
        continue
    parts = line.strip().split("\t")
    if len(parts) < 5:
        continue
    seq_id    = parts[0]
    analysis  = parts[3]
    sig_acc   = parts[4]
    sig_desc  = parts[5] if len(parts) > 5 else ""

    if analysis in ("Pfam", "PANTHER", "SMART"):
        seq_domains[seq_id].add(f"{sig_acc} ({sig_desc[:40]})")
        domain_counts[sig_acc] += 1

    if sig_acc == "PF00089":
        pfam_present[seq_id] = True
        tryp_spc_count += 1

print(f"\nTop 10 domínios mais frequentes:")
for dom, cnt in domain_counts.most_common(10):
    print(f"  {cnt:>4}x {dom}")

print(f"\nSequências com PF00089 (Tryp_SPc): {len(pfam_present)}/{len(seq_domains)}")

# Salvar sumário por sequência
summary_path = Path("${IPS_OUT}/domain_summary.tsv")
with open(summary_path, "w") as f:
    f.write("id\tn_domains\thas_PF00089\tdomains\n")
    for seq_id, doms in sorted(seq_domains.items()):
        f.write(f"{seq_id}\t{len(doms)}\t{seq_id in pfam_present}\t{'; '.join(sorted(doms))}\n")
print(f"Sumário por sequência: {summary_path}")
PYEOF

echo ""
echo "✅ InterProScan concluído."
echo "   Resultados em: ${IPS_OUT}/"
echo ""
echo "Próximo passo:"
echo "   python scripts/phase5_primary_char/validate.py"
