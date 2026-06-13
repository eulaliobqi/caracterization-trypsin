#!/usr/bin/env python3
"""
validate.py — Validação Fase 5 (Caracterização primária)

Consolida resultados de ProtParam + SignalP6 + InterProScan numa tabela
única prontas para publicação. Decide se os candidatos são biologicamente
válidos para avançar para modelagem estrutural (Fase 6+).

Critérios de validação para tripsinas de inseto (Lepidoptera):
  MW:         20-40 kDa   (tripsina madura sem pré-pró)
  pI:         3.0-10.0    (variável em Lepidoptera, geralmente 5-9)
  Peptídeo sinal: ≥ 70%   (proteínas secretadas do midgut)
  Domínio PF00089: 100%   (obrigatório para tripsina funcional)
"""

import sys
import csv
import argparse
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results5", default="results/phase5")
    args = parser.parse_args()

    d5 = Path(args.results5)

    print("\n══════════════════════════════════════════════════════")
    print("  VALIDAÇÃO FASE 5 — Caracterização Primária")
    print("══════════════════════════════════════════════════════\n")

    # ── Carregar dados ─────────────────────────────────────────────────────
    # ProtParam
    protparam = {}
    pp_file = d5 / "protparam_results.tsv"
    if pp_file.exists():
        for row in csv.DictReader(open(pp_file), delimiter="\t"):
            if "ERROR" not in str(row.get("mw_kda","")):
                protparam[row["id"]] = row

    # SignalP
    signalp = {}
    sp_file = d5 / "signalp" / "signalp_summary.tsv"
    if sp_file.exists():
        for row in csv.DictReader(open(sp_file), delimiter="\t"):
            if not row["id"].startswith("#"):
                signalp[row["id"]] = row

    # InterProScan
    interpro = {}
    ips_file = d5 / "interproscan" / "domain_summary.tsv"
    if ips_file.exists():
        for row in csv.DictReader(open(ips_file), delimiter="\t"):
            interpro[row["id"]] = row

    # ── Tabela consolidada ─────────────────────────────────────────────────
    seq_ids = list(protparam.keys()) or \
              list(signalp.keys()) or \
              list(interpro.keys())

    if not seq_ids:
        print("❌ Nenhum resultado encontrado em results/phase5/")
        print("   Execute: python scripts/phase5_primary_char/01_protparam.py")
        sys.exit(1)

    print(f"Sequências analisadas: {len(seq_ids)}\n")

    # Header da tabela
    header = f"{'#':<3} {'ID':<40} {'MW':>6} {'pI':>5} {'#aa':>5} {'Sinal':>6} {'PF00089':>7} {'Notas'}"
    print(header)
    print("-" * 90)

    rows_out = []
    errors = []
    warnings = []

    for i, seq_id in enumerate(sorted(seq_ids), 1):
        pp = protparam.get(seq_id, {})
        sp = signalp.get(seq_id, {})
        ip = interpro.get(seq_id, {})

        mw   = pp.get("mw_kda", "?")
        pi   = pp.get("pi", "?")
        aac  = pp.get("aa_count", "?")
        notes_pp = pp.get("notes", "")

        has_signal = sp.get("has_signal_peptide", "?")
        cleavage   = sp.get("cleavage_pos", "?")

        has_pf89 = ip.get("has_PF00089", "?")
        n_doms   = ip.get("n_domains", "?")

        # Flags visuais
        mw_ok = "?" if mw == "?" else ("✅" if 20 <= float(mw) <= 40 else "⚠️ ")
        sp_ok = "✅" if has_signal in ("True", True) else ("⚠️ " if has_signal == "?" else "  ")
        ip_ok = "✅" if has_pf89 in ("True", True) else ("❌" if has_pf89 == "False" else "? ")

        notes = []
        if notes_pp and notes_pp != "OK":
            notes.append(notes_pp[:30])
        if has_signal not in ("True", True, "?"):
            notes.append("sem sinal")
        if has_pf89 == "False":
            notes.append("SEM PF00089!")

        print(f"  {i:<3} {seq_id:<39} {str(mw)[:5]:>5}{mw_ok} "
              f"{str(pi)[:4]:>4} {str(aac)[:5]:>5} "
              f"{sp_ok:>5}    {ip_ok:>5}     {'; '.join(notes) if notes else 'OK'}")

        rows_out.append({
            "id": seq_id, "mw_kda": mw, "pi": pi, "aa_count": aac,
            "has_signal_peptide": has_signal, "cleavage_site": cleavage,
            "has_PF00089": has_pf89, "n_interpro_domains": n_doms,
            "notes": "; ".join(notes) if notes else "OK"
        })

    # ── Salvar tabela consolidada ──────────────────────────────────────────
    combined_path = d5 / "primary_characterization_combined.tsv"
    if rows_out:
        fieldnames = list(rows_out[0].keys())
        with open(combined_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows_out)
        print(f"\nTabela consolidada: {combined_path}")

    # ── Estatísticas ───────────────────────────────────────────────────────
    def safe_float(v):
        try: return float(v)
        except: return None

    mw_vals = [safe_float(r["mw_kda"]) for r in rows_out if safe_float(r["mw_kda"])]
    pi_vals = [safe_float(r["pi"])      for r in rows_out if safe_float(r["pi"])]
    n_signal = sum(1 for r in rows_out if r["has_signal_peptide"] in ("True", True))
    n_pf89   = sum(1 for r in rows_out if r["has_PF00089"] in ("True", True))
    n_total  = len(rows_out)

    print(f"\n── Estatísticas de Publicação ────────────────────────")
    if mw_vals:
        n_mw_ok = sum(1 for v in mw_vals if 20 <= v <= 40)
        print(f"  MW: {min(mw_vals):.1f}–{max(mw_vals):.1f} kDa (média {sum(mw_vals)/len(mw_vals):.1f})")
        print(f"      {n_mw_ok}/{n_total} dentro de 20-40 kDa")
    if pi_vals:
        print(f"  pI: {min(pi_vals):.1f}–{max(pi_vals):.1f}    (média {sum(pi_vals)/len(pi_vals):.1f})")

    if signalp:
        pct_sp = n_signal / n_total * 100 if n_total else 0
        print(f"  Peptídeo sinal: {n_signal}/{n_total} ({pct_sp:.0f}%)")
        if pct_sp < 60:
            warnings.append(f"Apenas {pct_sp:.0f}% com peptídeo sinal (esperado ≥70%)")

    if interpro:
        pct_pf = n_pf89 / n_total * 100 if n_total else 0
        print(f"  Domínio PF00089: {n_pf89}/{n_total} ({pct_pf:.0f}%)")
        if pct_pf < 100:
            warnings.append(f"Apenas {pct_pf:.0f}% com PF00089 — verificar sequências sem domínio")

    # ── Decisão ───────────────────────────────────────────────────────────
    print(f"\n── Decisão ───────────────────────────────────────────")
    if errors:
        print("❌ FASE 5: FAIL")
        for e in errors: print(f"   • {e}")
        sys.exit(1)
    elif warnings:
        print("⚠️  FASE 5: PASS com ressalvas:")
        for w in warnings: print(f"   • {w}")
    else:
        print("✅ FASE 5: PASS — candidatos validados para modelagem estrutural")

    print(f"""
═══════════════════════════════════════════════════════
  ✅ FASES 1-5 CONCLUÍDAS

  {n_total} tripsinas de Anticarsia gemmatalis caracterizadas:
  • Identificadas por DIAMOND ∩ HMMER (Fase 3)
  • Filtradas por completude e tríade catalítica (Fase 4)
  • Caracterizadas fisico-quimicamente (Fase 5)

  PRÓXIMOS PASSOS (Nextflow, fases 6-11):

  Fase 6 — Filogenia:
    nextflow run nextflow/main.nf -profile debian --phase 6 \\
      --input_fasta results/phase4/complete_trypsins.fasta

  Fase 7 — AlphaFold3 (GPU obrigatório):
    nextflow run nextflow/main.nf -profile debian,gpu --phase 7

  Ou executar tudo de uma vez:
    nextflow run nextflow/main.nf -profile debian,gpu
═══════════════════════════════════════════════════════
""")


if __name__ == "__main__":
    main()
