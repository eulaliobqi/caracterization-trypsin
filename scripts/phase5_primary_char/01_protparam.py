#!/usr/bin/env python3
"""
01_protparam.py — Fase 5: Caracterização físico-química (ProtParam)

Wrapper que chama nextflow/bin/compute_protparam.py e gera tabela formatada.

Input:  results/phase4/complete_trypsins.fasta
Output: results/phase5/
  - protparam_results.tsv     → tabela completa
  - protparam_summary.txt     → resumo para publicação
"""

import sys
import subprocess
import argparse
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input",    default="results/phase4/complete_trypsins.fasta")
    parser.add_argument("--results5", default="results/phase5")
    parser.add_argument("--bin-dir",  default="nextflow/bin")
    args = parser.parse_args()

    out_dir = Path(args.results5)
    out_dir.mkdir(parents=True, exist_ok=True)

    input_fasta = Path(args.input)
    if not input_fasta.exists():
        print(f"❌ {input_fasta} não encontrado — execute Fase 4 primeiro")
        sys.exit(1)

    # Localizar compute_protparam.py
    script = Path(args.bin_dir) / "compute_protparam.py"
    if not script.exists():
        for candidate in [
            Path("nextflow/bin/compute_protparam.py"),
            Path("../nextflow/bin/compute_protparam.py"),
        ]:
            if candidate.exists():
                script = candidate
                break

    if not script.exists():
        print(f"❌ compute_protparam.py não encontrado em {script}")
        sys.exit(1)

    print("\n══════════════════════════════════════════════")
    print("  FASE 5 — Passo 1: ProtParam")
    print("══════════════════════════════════════════════\n")

    output_tsv = out_dir / "protparam_results.tsv"

    result = subprocess.run([
        sys.executable, str(script),
        "--input",  str(input_fasta),
        "--output", str(output_tsv),
    ])

    if result.returncode != 0:
        print("❌ compute_protparam.py falhou")
        sys.exit(1)

    # ── Resumo formatado ──────────────────────────────────────────────────
    import csv
    rows = list(csv.DictReader(open(output_tsv), delimiter="\t"))

    if not rows:
        print("❌ Tabela protparam vazia")
        sys.exit(1)

    # Tabela para publicação (LaTeX/TSV simplificado)
    summary_path = out_dir / "protparam_summary.txt"
    with open(summary_path, "w") as f:
        f.write("# ProtParam Summary — Anticarsia gemmatalis trypsins\n")
        f.write("# Gerado por: scripts/phase5_primary_char/01_protparam.py\n\n")
        f.write(f"{'ID':<45} {'MW(kDa)':>8} {'pI':>5} {'#aa':>5} {'Instab':>7} {'GRAVY':>7} {'Notes'}\n")
        f.write("-" * 100 + "\n")
        for row in rows:
            if "ERROR" in str(row.get("mw_kda", "")):
                continue
            f.write(f"  {row['id']:<43} {row.get('mw_kda','?'):>8} {row.get('pi','?'):>5} "
                    f"{row.get('aa_count','?'):>5} {row.get('instability_index','?'):>7} "
                    f"{row.get('gravy','?'):>7}  {row.get('notes','OK')}\n")

    # Estatísticas
    def safe_float(v, default=0.0):
        try: return float(v)
        except: return default

    mws  = [safe_float(r.get("mw_kda",0))  for r in rows if "ERROR" not in str(r.get("mw_kda",""))]
    pis  = [safe_float(r.get("pi",0))       for r in rows if "ERROR" not in str(r.get("pi",""))]
    lens = [int(r.get("aa_count",0))        for r in rows if r.get("aa_count","").isdigit()]

    if mws and pis:
        print(f"\n── Estatísticas Gerais ────────────────────────────")
        print(f"  Sequências analisadas: {len(mws)}")
        print(f"  MW: {min(mws):.1f}–{max(mws):.1f} kDa  (média: {sum(mws)/len(mws):.1f})")
        print(f"  pI: {min(pis):.1f}–{max(pis):.1f}       (média: {sum(pis)/len(pis):.1f})")
        if lens:
            print(f"  aa: {min(lens)}–{max(lens)}          (média: {sum(lens)//len(lens)})")

        n_ok = sum(1 for r in rows if r.get("mw_in_range")=="True" and r.get("pi_in_range")=="True")
        print(f"\n  Dentro do range esperado (MW 20-40 kDa, pI 3-7): {n_ok}/{len(rows)}")

    print(f"\n── Arquivos gerados ────────────────────────────────")
    print(f"  {output_tsv}")
    print(f"  {summary_path}")
    print(f"\n✅ ProtParam concluído.")
    print(f"\nPróximos passos:")
    print(f"   bash scripts/phase5_primary_char/02_signalp.sh")
    print(f"   bash scripts/phase5_primary_char/03_interproscan.sh")
    print(f"   python scripts/phase5_primary_char/validate.py")


if __name__ == "__main__":
    main()
