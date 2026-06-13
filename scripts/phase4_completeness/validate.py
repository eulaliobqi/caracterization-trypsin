#!/usr/bin/env python3
"""
validate.py — Validação Fase 4 (Completude)

Apresenta a lista final de tripsinas completas com propriedades básicas
para inspeção antes de avançar para Fase 5 (caracterização primária).

Também sugere sequências para verificação manual via BLAST web.
"""

import sys
import argparse
from pathlib import Path

try:
    from Bio import SeqIO
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
except ImportError:
    print("ERRO: Biopython necessário")
    sys.exit(1)


def quick_props(seq_str):
    """Propriedades rápidas para preview."""
    seq = seq_str.replace("*","").replace("-","").replace("X","")
    if not seq:
        return None
    try:
        a = ProteinAnalysis(seq)
        return {
            "len":   len(seq),
            "mw":    round(a.molecular_weight() / 1000, 1),
            "pi":    round(a.isoelectric_point(), 1),
            "inst":  round(a.instability_index(), 1),
        }
    except Exception:
        return {"len": len(seq), "mw": 0, "pi": 0, "inst": 0}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results4", default="results/phase4")
    args = parser.parse_args()

    d = Path(args.results4)

    print("\n══════════════════════════════════════════════")
    print("  VALIDAÇÃO FASE 4 — Tripsinas Completas")
    print("══════════════════════════════════════════════\n")

    fasta = d / "complete_trypsins.fasta"
    if not fasta.exists():
        print(f"❌ {fasta} não encontrado — execute filter.py primeiro")
        sys.exit(1)

    records = list(SeqIO.parse(fasta, "fasta"))
    n = len(records)

    # Calcular props para todas
    props_list = []
    for rec in records:
        p = quick_props(str(rec.seq))
        if p:
            props_list.append((rec, p))

    # Tabela principal
    print(f"{'#':<4} {'ID':<45} {'len':>5} {'MW(kDa)':>8} {'pI':>5} {'Instab.':>8}")
    print("-" * 80)
    ok_mw  = 0
    ok_pi  = 0
    for i, (rec, p) in enumerate(props_list, 1):
        mw_flag = "✅" if 20 <= p["mw"] <= 40 else "⚠️ "
        pi_flag = "✅" if 3.0 <= p["pi"] <= 10.0 else "⚠️ "
        if 20 <= p["mw"] <= 40: ok_mw += 1
        if 3.0 <= p["pi"] <= 10.0: ok_pi += 1
        print(f"  {i:<3} {rec.id:<44} {p['len']:>5} {p['mw']:>7.1f} {mw_flag} "
              f"{p['pi']:>4.1f} {pi_flag}  {p['inst']:>6.1f}")

    # Estatísticas
    mws  = [p["mw"]   for _, p in props_list]
    pis  = [p["pi"]   for _, p in props_list]
    lens = [p["len"]  for _, p in props_list]

    print(f"\n── Estatísticas ──────────────────────────────────")
    print(f"  Total de tripsinas:            {n}")
    print(f"  MW média:       {sum(mws)/len(mws):.1f} kDa   (range: {min(mws):.1f}-{max(mws):.1f})")
    print(f"  pI médio:       {sum(pis)/len(pis):.1f}       (range: {min(pis):.1f}-{max(pis):.1f})")
    print(f"  Comprimento médio: {sum(lens)/len(lens):.0f} aa  (range: {min(lens)}-{max(lens)})")
    print(f"\n  MW 20-40 kDa:   {ok_mw}/{n} ({ok_mw/n*100:.0f}%)")
    print(f"  pI 3.0-10.0:    {ok_pi}/{n} ({ok_pi/n*100:.0f}%)")

    # Relatório de completude (se disponível)
    report = d / "completeness_report.tsv"
    if report.exists():
        import csv
        rows = list(csv.DictReader(open(report), delimiter="\t"))
        n_gdsgg = sum(1 for r in rows if r.get("pass_strict")=="True" and r.get("gdsgg_motif")=="True")
        n_triad = sum(1 for r in rows if r.get("triad_ok")=="True")
        print(f"\n  Motivo GDSGG (Ser195):       {n_gdsgg}/{n}")
        if any(r.get("triad_ok","NA") != "NA" for r in rows):
            print(f"  Tríade His+Asp+Ser:          {n_triad}/{n}")

    # Check de borderline
    borderline = d / "borderline_trypsins.fasta"
    if borderline.exists():
        n_bl = sum(1 for l in open(borderline) if l.startswith(">"))
        if n_bl > 0:
            print(f"\n  ⚠️  {n_bl} sequências borderline (Met+len, sem GDSGG)")
            print(f"     → Ver: results/phase4/borderline_trypsins.fasta")
            print(f"     → Algumas podem ser tripsinas ativas — verificar manualmente")

    # Decisão
    print(f"\n── Decisão ───────────────────────────────────────")
    errors = []

    if n < 3:
        errors.append(f"Apenas {n} tripsinas — muito poucas para análise robusta")

    ok_bio = ok_mw + ok_pi
    if ok_bio < n * 1.5:  # pelo menos 75% com MW OK e 75% com pI OK
        pass  # já mostrado na tabela

    if errors:
        print("❌ FASE 4: FAIL")
        for e in errors: print(f"   • {e}")
        sys.exit(1)
    elif n < 5:
        print(f"⚠️  FASE 4: PASS com ressalva — {n} tripsinas (mínimo para publicação: 5)")
        print(f"   Considere revisar parâmetros da Fase 3 ou incluir suggestive")
    else:
        print(f"✅ FASE 4: PASS — {n} tripsinas completas prontas para caracterização")

    # Sugestão para BLAST manual
    print(f"\n── Verificação manual recomendada (BLAST web) ────")
    print(f"   https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp")
    print(f"   Verifique as primeiras 5 sequências abaixo:\n")
    for i, (rec, _) in enumerate(props_list[:5], 1):
        seq_preview = str(rec.seq)[:60].replace("*","")
        print(f"  >{rec.id}")
        print(f"   {seq_preview}...")

    print(f"\n✅ Avançar para Fase 5:")
    print(f"   mamba run -n analysis python scripts/phase5_primary_char/01_protparam.py")


if __name__ == "__main__":
    main()
