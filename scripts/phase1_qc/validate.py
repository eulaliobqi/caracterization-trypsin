#!/usr/bin/env python3
"""
validate.py — Validação da Fase 1 (QC do Assembly)

Critérios para PASS e avançar para Fase 2:
  - BUSCO completeness ≥ 70%        (mínimo aceitável para análise)
  - BUSCO completeness ≥ 80%        (ideal para publicação)
  - Redundância removida 5-40%      (CD-HIT razoável)
  - Sequências totais > 10.000      (assembly com cobertura adequada)

Uso:
  python scripts/phase1_qc/validate.py
  python scripts/phase1_qc/validate.py --results results/phase1
"""

import sys
import re
import argparse
from pathlib import Path

# Thresholds
BUSCO_PASS_MIN   = 70.0   # % completude mínima para continuar
BUSCO_IDEAL_MIN  = 80.0   # % ideal para journal Q1
MIN_SEQS         = 10_000 # sequências mínimas no assembly

def parse_busco_summary(summary_path: Path) -> dict:
    """Extrai métricas do short_summary do BUSCO."""
    text = summary_path.read_text()
    metrics = {}

    # Padrão: "C:85.2%[S:70.1%,D:15.1%],F:5.3%,M:9.5%,n:1367"
    m = re.search(r'C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%,n:(\d+)', text)
    if m:
        metrics['complete']    = float(m.group(1))
        metrics['single']      = float(m.group(2))
        metrics['duplicated']  = float(m.group(3))
        metrics['fragmented']  = float(m.group(4))
        metrics['missing']     = float(m.group(5))
        metrics['n_busco']     = int(m.group(6))
        return metrics

    # Formato alternativo (BUSCO 5.x)
    patterns = {
        'complete':   r'(\d+)\s+Complete BUSCOs',
        'single':     r'(\d+)\s+Complete and single-copy',
        'duplicated': r'(\d+)\s+Complete and duplicated',
        'fragmented': r'(\d+)\s+Fragmented',
        'missing':    r'(\d+)\s+Missing',
        'n_busco':    r'(\d+)\s+Total BUSCO groups',
    }
    for key, pat in patterns.items():
        mm = re.search(pat, text)
        if mm:
            metrics[key] = int(mm.group(1))

    # Calcular porcentagens se temos contagens
    if 'n_busco' in metrics and metrics['n_busco'] > 0:
        n = metrics['n_busco']
        for k in ['complete', 'single', 'duplicated', 'fragmented', 'missing']:
            if k in metrics:
                metrics[f'{k}_pct'] = round(metrics[k] / n * 100, 1)
        metrics['complete'] = metrics.get('complete_pct', 0)

    return metrics


def count_fasta_seqs(fasta_path: Path) -> int:
    """Conta sequências num arquivo FASTA."""
    return sum(1 for line in fasta_path.open() if line.startswith(">"))


def main():
    parser = argparse.ArgumentParser(description="Validação Fase 1 — QC Assembly")
    parser.add_argument("--results", default="results/phase1",
                        help="Diretório de resultados da Fase 1")
    args = parser.parse_args()

    results_dir = Path(args.results)
    if not results_dir.exists():
        print(f"❌ Diretório não encontrado: {results_dir}")
        print(f"   Rode primeiro: bash scripts/phase1_qc/run.sh")
        sys.exit(1)

    checks = []
    warnings = []
    errors = []

    # ── 1. BUSCO ──────────────────────────────────────────────────────
    busco_summary = list(results_dir.glob("busco_*/run_*/short_summary*.txt"))
    if not busco_summary:
        busco_summary = list(results_dir.glob("*.txt"))
        busco_summary = [f for f in busco_summary if "busco" in f.name.lower()]

    # Também checar arquivo copiado
    busco_copy = results_dir / "busco_short_summary.txt"
    if busco_copy.exists():
        busco_summary = [busco_copy] + busco_summary

    print("\n══════════════════════════════════════════════")
    print("  VALIDAÇÃO FASE 1 — QC do Assembly")
    print("══════════════════════════════════════════════\n")

    if not busco_summary:
        errors.append("BUSCO summary não encontrado — execute run.sh primeiro")
        checks.append(("BUSCO completeness", "N/A", "❌ FALTANDO"))
    else:
        metrics = parse_busco_summary(busco_summary[0])
        pct = metrics.get('complete', 0)

        print(f"BUSCO ({metrics.get('n_busco', '?')} genes insecta_odb10):")
        print(f"  Completos:     {pct:.1f}%  (single={metrics.get('single', '?')}%, dup={metrics.get('duplicated', '?')}%)")
        print(f"  Fragmentados:  {metrics.get('fragmented', '?')}%")
        print(f"  Ausentes:      {metrics.get('missing', '?')}%")

        if pct >= BUSCO_IDEAL_MIN:
            status = f"✅ PASS ({pct:.1f}% ≥ {BUSCO_IDEAL_MIN}% ideal)"
        elif pct >= BUSCO_PASS_MIN:
            status = f"⚠️  PASS-MIN ({pct:.1f}% ≥ {BUSCO_PASS_MIN}% mínimo, mas <{BUSCO_IDEAL_MIN}% ideal)"
            warnings.append(f"BUSCO {pct:.1f}% abaixo do ideal ({BUSCO_IDEAL_MIN}%). "
                           "Considere assembly de maior qualidade para publicação Q1.")
        else:
            status = f"❌ FAIL ({pct:.1f}% < {BUSCO_PASS_MIN}% mínimo)"
            errors.append(f"BUSCO completeness muito baixa: {pct:.1f}%")

        checks.append(("BUSCO completeness", f"{pct:.1f}%", status))

    # ── 2. CD-HIT ─────────────────────────────────────────────────────
    nr_fasta = next(results_dir.glob("assembly_nr*.fasta"), None)
    if nr_fasta is None:
        errors.append("assembly_nr*.fasta não encontrado — execute run.sh primeiro")
        checks.append(("CD-HIT output", "N/A", "❌ FALTANDO"))
    else:
        n_nr = count_fasta_seqs(nr_fasta)
        identity = nr_fasta.stem.replace("assembly_nr", "")
        print(f"\nCD-HIT-EST ({identity}% identidade):")
        print(f"  Sequências após clusterização: {n_nr:,}")

        if n_nr < MIN_SEQS:
            status = f"⚠️  AVISO ({n_nr:,} < {MIN_SEQS:,} esperado)"
            warnings.append(f"Apenas {n_nr:,} sequências após CD-HIT. Verifique qualidade do assembly.")
        else:
            status = f"✅ OK ({n_nr:,} sequências)"

        checks.append(("CD-HIT sequências", f"{n_nr:,}", status))

    # ── 3. Assembly original ───────────────────────────────────────────
    # Tentar encontrar via config ou path padrão
    assembly_paths = [
        Path.home() / "gromacs/caracterization-trypsin/data/raw/trinity_assembly.fasta",
        Path("data/raw/trinity_assembly.fasta"),
    ]
    for p in assembly_paths:
        if p.exists():
            n_orig = count_fasta_seqs(p)
            print(f"\nAssembly original: {n_orig:,} transcritos")
            if n_nr:
                red_pct = (1 - n_nr / n_orig) * 100
                print(f"Redundância removida: {red_pct:.1f}%")
                if 5 <= red_pct <= 40:
                    checks.append(("Redundância CD-HIT", f"{red_pct:.1f}%", "✅ OK (range 5-40%)"))
                else:
                    checks.append(("Redundância CD-HIT", f"{red_pct:.1f}%",
                                   f"⚠️  Fora do range 5-40% (verificar parâmetros)"))
                    warnings.append(f"Redundância de {red_pct:.1f}% incomum. Esperado 5-40%.")
            break

    # ── Tabela de resultados ───────────────────────────────────────────
    print("\n── Tabela de Validação ──────────────────────────────")
    print(f"{'Critério':<30} {'Valor':<15} {'Status'}")
    print("-" * 70)
    for criterion, value, status in checks:
        print(f"  {criterion:<28} {value:<15} {status}")

    # ── Decisão final ─────────────────────────────────────────────────
    print("\n── Decisão ──────────────────────────────────────────")
    if errors:
        print("❌ FASE 1: FAIL — Corrigir antes de avançar:")
        for e in errors:
            print(f"   • {e}")
        sys.exit(1)
    elif warnings:
        print("⚠️  FASE 1: PASS com avisos — pode avançar, mas considere:")
        for w in warnings:
            print(f"   • {w}")
        print("\n✅ Avançar para Fase 2:")
        print("   bash scripts/phase2_orf/run.sh")
    else:
        print("✅ FASE 1: PASS — Assembly de boa qualidade")
        print("\nPróximo passo:")
        print("   bash scripts/phase2_orf/run.sh")


if __name__ == "__main__":
    main()
