#!/usr/bin/env python3
"""
03_intersect.py — Interseção DIAMOND ∩ HMMER para identificação de tripsinas

Estratégia de dupla validação:
  CONFIDENT  = DIAMOND hit + HMMER domínio PF00089  (alta confiança)
  SUGGESTIVE = só DIAMOND  OU só HMMER              (menor confiança)

Input:
  results/phase3/diamond_trypsin_ids.txt
  results/phase3/hmmer_trypsin_ids.txt
  results/phase2/agemmatalis.pep

Output:
  results/phase3/trypsins_confident.fasta
  results/phase3/trypsins_suggestive.fasta
  results/phase3/trypsin_ids_confident.txt
  results/phase3/identification_report.tsv
"""

import sys
import argparse
from pathlib import Path
import subprocess

def main():
    parser = argparse.ArgumentParser(description="Interseção DIAMOND + HMMER")
    parser.add_argument("--results3",    default="results/phase3")
    parser.add_argument("--pep",         default="results/phase2/agemmatalis.pep")
    parser.add_argument("--diamond-only", action="store_true",
                        help="Usar apenas DIAMOND (sem HMMER)")
    parser.add_argument("--bin-dir",     default="nextflow/bin",
                        help="Diretório com filter_complete_trypsins.py")
    args = parser.parse_args()

    out_dir    = Path(args.results3)
    pep_file   = Path(args.pep)
    bin_dir    = Path(args.bin_dir)
    script     = bin_dir / "filter_complete_trypsins.py"

    if not script.exists():
        # Tentar path relativo ao repo
        for candidate in [
            Path("nextflow/bin/filter_complete_trypsins.py"),
            Path("../nextflow/bin/filter_complete_trypsins.py"),
        ]:
            if candidate.exists():
                script = candidate
                break

    if not script.exists():
        print(f"❌ filter_complete_trypsins.py não encontrado. Esperado em: {script}")
        sys.exit(1)

    diamond_ids = out_dir / "diamond_trypsin_ids.txt"
    hmmer_ids   = out_dir / "hmmer_trypsin_ids.txt"

    print("\n══════════════════════════════════════════════")
    print("  FASE 3 — Passo 3: Interseção DIAMOND ∩ HMMER")
    print("══════════════════════════════════════════════\n")

    # Verificar inputs
    for f, desc in [(diamond_ids, "DIAMOND IDs"), (pep_file, "Proteínas Fase 2")]:
        if not f.exists():
            print(f"❌ {desc} não encontrado: {f}")
            sys.exit(1)

    if not hmmer_ids.exists() and not args.diamond_only:
        print(f"⚠️  HMMER IDs não encontrado: {hmmer_ids}")
        print("   Rodando em modo DIAMOND-only. Use --diamond-only para suprimir este aviso.")
        args.diamond_only = True

    # Contar inputs
    n_diamond = sum(1 for _ in open(diamond_ids))
    n_hmmer   = sum(1 for _ in open(hmmer_ids)) if hmmer_ids.exists() else 0

    print(f"DIAMOND candidatos: {n_diamond}")
    print(f"HMMER  candidatos:  {n_hmmer}")

    # Chamar filter_complete_trypsins.py
    cmd = [
        sys.executable, str(script),
        "--diamond",           str(diamond_ids),
        "--pep",               str(pep_file),
        "--output_confident",  str(out_dir / "trypsins_confident.fasta"),
        "--output_suggestive", str(out_dir / "trypsins_suggestive.fasta"),
        "--ids_confident",     str(out_dir / "trypsin_ids_confident.txt"),
        "--report",            str(out_dir / "identification_report.tsv"),
    ]
    if hmmer_ids.exists() and not args.diamond_only:
        cmd += ["--hmm_ids", str(hmmer_ids)]

    print(f"\nRodando: {' '.join(cmd)}\n")
    result = subprocess.run(cmd, capture_output=False)
    if result.returncode != 0:
        print("❌ filter_complete_trypsins.py falhou")
        sys.exit(1)

    # Contagens finais
    confident_fasta   = out_dir / "trypsins_confident.fasta"
    suggestive_fasta  = out_dir / "trypsins_suggestive.fasta"
    confident_ids_out = out_dir / "trypsin_ids_confident.txt"

    n_confident  = sum(1 for l in open(confident_fasta)  if l.startswith(">")) if confident_fasta.exists() else 0
    n_suggestive = sum(1 for l in open(suggestive_fasta) if l.startswith(">")) if suggestive_fasta.exists() else 0

    print(f"\n── Resultado da Interseção ──────────────────────────")
    print(f"  Confident  (DIAMOND ∩ HMMER): {n_confident}")
    print(f"  Suggestive (DIAMOND | HMMER): {n_suggestive}")
    print(f"  Total candidatos:              {n_confident + n_suggestive}")

    # Validação biológica rápida
    print(f"\n── Expectativa biológica ────────────────────────────")
    print(f"  Esperado (Lepidoptera): 30–200 confident, 5–30 completas após Fase 4")

    if n_confident == 0:
        print("❌ ZERO sequências confident! Verificar:")
        print("   1. Os DBs (UniProt, Pfam) estão corretos?")
        print("   2. As ORFs da Fase 2 foram preditas corretamente?")
        sys.exit(1)
    elif n_confident < 10:
        print(f"⚠️  Poucos candidatos confident ({n_confident}). Pode estar OK para assembly pequeno.")
    elif n_confident > 500:
        print(f"⚠️  Muitos candidatos confident ({n_confident}). Verificar E-values.")
    else:
        print(f"✅ {n_confident} candidatos confident — range esperado para Lepidoptera")

    print(f"\n✅ Interseção concluída. Próximo passo:")
    print(f"   python scripts/phase3_trypsin_id/validate.py")


if __name__ == "__main__":
    main()
