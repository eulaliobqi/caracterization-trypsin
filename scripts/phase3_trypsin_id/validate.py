#!/usr/bin/env python3
"""
validate.py — Validação Fase 3 (Identificação de tripsinas)

Mostra tabela comparativa DIAMOND vs HMMER vs Intersecção
e decide se os números são biologicamente plausíveis.
"""

import sys
import argparse
from pathlib import Path

try:
    from Bio import SeqIO
except ImportError:
    print("ERRO: Biopython necessário")
    sys.exit(1)


def count_fasta(path):
    return sum(1 for l in open(path) if l.startswith(">")) if Path(path).exists() else 0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results3", default="results/phase3")
    args = parser.parse_args()
    d = Path(args.results3)

    print("\n══════════════════════════════════════════════")
    print("  VALIDAÇÃO FASE 3 — Identificação de Tripsinas")
    print("══════════════════════════════════════════════\n")

    # Contagens
    n_diamond   = sum(1 for _ in open(d/"diamond_trypsin_ids.txt")) if (d/"diamond_trypsin_ids.txt").exists() else 0
    n_hmmer     = sum(1 for _ in open(d/"hmmer_trypsin_ids.txt"))   if (d/"hmmer_trypsin_ids.txt").exists()   else 0
    n_confident = count_fasta(d/"trypsins_confident.fasta")
    n_suggest   = count_fasta(d/"trypsins_suggestive.fasta")

    print(f"{'Método':<35} {'Candidatos':>12}")
    print("-" * 50)
    print(f"  {'DIAMOND vs UniProt':<33} {n_diamond:>12,}")
    print(f"  {'HMMER vs PF00089 (Tryp_SPc)':<33} {n_hmmer:>12,}")
    print(f"  {'Confident (DIAMOND ∩ HMMER)':<33} {n_confident:>12,}")
    print(f"  {'Suggestive (DIAMOND | HMMER)':<33} {n_suggest:>12,}")

    # Análise da distribuição de comprimentos dos confident
    conf_fasta = d / "trypsins_confident.fasta"
    if conf_fasta.exists() and n_confident > 0:
        lengths = [len(str(r.seq).replace("*","")) for r in SeqIO.parse(conf_fasta, "fasta")]
        avg_len = sum(lengths) / len(lengths)
        print(f"\nComprimento médio (confident): {avg_len:.0f} aa")
        print(f"Range: {min(lengths)}-{max(lengths)} aa")

        # % dentro do range esperado para tripsinas
        in_range = sum(1 for l in lengths if 200 <= l <= 400)
        print(f"Dentro de 200-400 aa: {in_range}/{n_confident} ({in_range/n_confident*100:.0f}%)")

    # Top anotações DIAMOND
    diamond_tsv = d / "diamond_trypsin.tsv"
    if diamond_tsv.exists():
        print(f"\n── Top 5 descrições DIAMOND ──────────────────────")
        from collections import Counter
        descs = []
        for line in open(diamond_tsv):
            parts = line.strip().split("\t")
            if len(parts) >= 13:
                descs.append(parts[12][:60])
        for desc, cnt in Counter(descs).most_common(5):
            print(f"  ({cnt:>4}x) {desc}")

    # Relatório de identificação
    report = d / "identification_report.tsv"
    if report.exists():
        import csv
        rows = list(csv.DictReader(open(report), delimiter="\t"))
        if rows and "confident" in rows[0]:
            n_both  = sum(1 for r in rows if r.get("diamond_hit")=="True" and r.get("hmm_hit")=="True")
            n_diamo = sum(1 for r in rows if r.get("diamond_hit")=="True" and r.get("hmm_hit")=="False")
            n_hmm   = sum(1 for r in rows if r.get("diamond_hit")=="False" and r.get("hmm_hit")=="True")
            print(f"\n── Diagrama de Venn ──────────────────────────────")
            print(f"  DIAMOND only:   {n_diamo}")
            print(f"  HMMER only:     {n_hmm}")
            print(f"  DIAMOND ∩ HMMER: {n_both}  ← confident")

    # Validação final
    print(f"\n── Validação ─────────────────────────────────────")
    errors = []
    warnings = []

    if n_confident == 0:
        errors.append("Zero sequências confident — pipeline falhou")
    elif n_confident < 10:
        warnings.append(f"Apenas {n_confident} confident (esperado ≥30 para Lepidoptera midgut)")
    elif n_confident > 500:
        warnings.append(f"{n_confident} confident é alto — verificar parâmetros E-value")
    else:
        print(f"✅ {n_confident} confident — dentro do range esperado para Lepidoptera")

    if n_hmmer == 0:
        warnings.append("HMMER não encontrou hits — verificar Pfam-A.hmm ou usar --diamond-only")

    for e in errors:
        print(f"  ❌ {e}")
    for w in warnings:
        print(f"  ⚠️  {w}")

    if errors:
        sys.exit(1)

    print(f"\n✅ Fase 3 concluída.")
    print(f"\nPróximo passo:")
    print(f"   python scripts/phase4_completeness/filter.py")


if __name__ == "__main__":
    main()
