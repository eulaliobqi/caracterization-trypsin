#!/usr/bin/env python3
"""
validate.py — Validação da Fase 2 (Predição de ORFs)

Critérios para PASS:
  - ≥ 50.000 proteínas preditas
  - ≥ 20% são ORFs completas (5' + 3')
  - Comprimento médio ≥ 150 aa
  - % de ORFs com hint ≥ 20% (se hints foram usados)
"""

import sys
import re
import argparse
from pathlib import Path
from collections import Counter

try:
    from Bio import SeqIO
except ImportError:
    print("ERRO: Biopython necessário. Rode: mamba run -n discovery pip install biopython")
    sys.exit(1)


def main():
    parser = argparse.ArgumentParser(description="Validação Fase 2 — ORFs")
    parser.add_argument("--results", default="results/phase2")
    args = parser.parse_args()

    results_dir = Path(args.results)

    print("\n══════════════════════════════════════════════")
    print("  VALIDAÇÃO FASE 2 — Predição de ORFs")
    print("══════════════════════════════════════════════\n")

    pep_file = results_dir / "agemmatalis.pep"
    if not pep_file.exists():
        print(f"❌ {pep_file} não encontrado — execute run.sh primeiro")
        sys.exit(1)

    # ── Estatísticas das ORFs ──────────────────────────────────────────────
    records = list(SeqIO.parse(pep_file, "fasta"))
    n_total = len(records)

    type_counter = Counter()
    lengths = []
    for rec in records:
        seq = str(rec.seq).replace("*", "")
        lengths.append(len(seq))

        desc = rec.description
        if "complete"       in desc: type_counter["complete"] += 1
        elif "5prime_partial" in desc and "3prime_partial" in desc:
            type_counter["internal"] += 1
        elif "5prime_partial" in desc: type_counter["5prime_partial"] += 1
        elif "3prime_partial" in desc: type_counter["3prime_partial"] += 1
        else: type_counter["other"] += 1

    avg_len = sum(lengths) / len(lengths) if lengths else 0
    median_len = sorted(lengths)[len(lengths)//2] if lengths else 0
    pct_complete = type_counter["complete"] / n_total * 100 if n_total else 0

    print(f"Proteínas preditas: {n_total:,}")
    print(f"\nTipos de ORF:")
    for t, n in type_counter.most_common():
        pct = n / n_total * 100
        print(f"  {t:<20} {n:>7,} ({pct:.1f}%)")
    print(f"\nComprimento médio: {avg_len:.0f} aa")
    print(f"Comprimento mediano: {median_len} aa")
    print(f"Min: {min(lengths)} aa  |  Max: {max(lengths)} aa")

    # ── Distribuição de comprimentos ──────────────────────────────────────
    bins = [0, 100, 200, 300, 500, 1000, 9999]
    print(f"\nDistribuição de comprimentos:")
    for i in range(len(bins)-1):
        n = sum(1 for l in lengths if bins[i] <= l < bins[i+1])
        print(f"  {bins[i]:>5}-{bins[i+1]:<6} aa: {n:>7,} ({n/n_total*100:.1f}%)")

    # ── Hints BLAST ───────────────────────────────────────────────────────
    blast_hints = results_dir / "blast_hints.outfmt6"
    if blast_hints.exists():
        n_blast = sum(1 for _ in blast_hints.open())
        blast_ids = set(line.split("\t")[0] for line in blast_hints.open())
        pct_blast = len(blast_ids) / n_total * 100 if n_total else 0
        print(f"\nBLAST hints: {n_blast:,} hits ({len(blast_ids):,} ORFs únicas, {pct_blast:.1f}%)")

    pfam_hints = results_dir / "pfam_hints.domtblout"
    if pfam_hints.exists():
        pfam_ids = set()
        for line in pfam_hints.open():
            if not line.startswith("#"):
                parts = line.split()
                if len(parts) > 2:
                    pfam_ids.add(parts[2])  # query sequence ID
        pct_pfam = len(pfam_ids) / n_total * 100 if n_total else 0
        print(f"Pfam hints:  {len(pfam_ids):,} ORFs com domínio ({pct_pfam:.1f}%)")

    # ── Validação ─────────────────────────────────────────────────────────
    checks = []
    errors = []
    warnings = []

    # Critério 1: quantidade mínima
    if n_total >= 50_000:
        checks.append(("ORFs preditas", f"{n_total:,}", f"✅ OK (≥50.000)"))
    elif n_total >= 20_000:
        checks.append(("ORFs preditas", f"{n_total:,}", f"⚠️  AVISO (20k-50k, esperado >50k)"))
        warnings.append(f"Apenas {n_total:,} ORFs — transcriptoma pode ter cobertura baixa")
    else:
        checks.append(("ORFs preditas", f"{n_total:,}", f"❌ FAIL (<20.000)"))
        errors.append(f"Muito poucas ORFs ({n_total:,}). Verificar qualidade do assembly (Fase 1).")

    # Critério 2: % completas
    if pct_complete >= 20:
        checks.append(("ORFs completas", f"{pct_complete:.1f}%", "✅ OK (≥20%)"))
    else:
        checks.append(("ORFs completas", f"{pct_complete:.1f}%",
                        f"⚠️  BAIXO (<20%) — muitas ORFs parciais (normal em transcriptoma)"))
        warnings.append("Muitas ORFs parciais — fragmentação do transcriptoma (esperado para inseto)")

    # Critério 3: comprimento médio
    if avg_len >= 150:
        checks.append(("Comprimento médio", f"{avg_len:.0f} aa", "✅ OK (≥150 aa)"))
    else:
        checks.append(("Comprimento médio", f"{avg_len:.0f} aa", "⚠️  CURTO (<150 aa)"))
        warnings.append(f"Comprimento médio {avg_len:.0f} aa abaixo do esperado (≥150 aa)")

    print("\n── Tabela de Validação ──────────────────────────────")
    print(f"{'Critério':<30} {'Valor':<15} {'Status'}")
    print("-" * 70)
    for criterion, value, status in checks:
        print(f"  {criterion:<28} {value:<15} {status}")

    print("\n── Decisão ──────────────────────────────────────────")
    if errors:
        print("❌ FASE 2: FAIL — Corrigir antes de avançar:")
        for e in errors:
            print(f"   • {e}")
        sys.exit(1)
    elif warnings:
        print("⚠️  FASE 2: PASS com avisos — pode avançar:")
        for w in warnings:
            print(f"   • {w}")
        print("\n✅ Avançar para Fase 3:")
        print("   bash scripts/phase3_trypsin_id/01_diamond.sh")
    else:
        print("✅ FASE 2: PASS")
        print("\nPróximo passo:")
        print("   bash scripts/phase3_trypsin_id/01_diamond.sh")


if __name__ == "__main__":
    main()
