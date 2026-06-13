#!/usr/bin/env python3
"""
filter.py — Fase 4: Filtro de completude das tripsinas candidatas

Critérios de completude (OBRIGATÓRIOS para avançar para modelagem estrutural):
  1. Começa com Met (M) — indica início de proteína completo
  2. Comprimento ≥ 220 aa — protease madura + pré-proteína
  3. Motivo GDSGG — marca o Ser195 do sítio ativo (Tríade catalítica)

Critério adicional (quando bovine_trypsin_1TGN.fasta disponível):
  4. Tríade His57-Asp102-Ser195 confirmada por alinhamento múltiplo

Input:  results/phase3/trypsins_confident.fasta
Output: results/phase4/
  - complete_trypsins.fasta    → tripsinas que passam TODOS os critérios
  - borderline_trypsins.fasta  → passam Met+comprimento mas sem motivo GDSGG
  - completeness_report.tsv   → detalhes por sequência
  - triad_check.tsv           → verificação da tríade (se ref disponível)
"""

import sys
import re
import argparse
import subprocess
from pathlib import Path

try:
    from Bio import SeqIO
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
except ImportError:
    print("ERRO: Biopython necessário. Rode: mamba run -n analysis pip install biopython")
    sys.exit(1)

# Motivos ao redor de Ser195 (numeração quimotripsina) — do mais ao menos conservado
# Nível 1 — canônico:    GDS[AG]G  (Ser195 intacto)
# Nível 2 — alt:         [GA]DS[GASTVNC]G  (variação nos flancos, Ser conservado)
# Nível 3 — mutante:     GD[^SP]GG ou [GA]D[^SP][GA]G  (Ser195 substituído,
#                         contexto GDXGG preservado — enzima provavelmente inativa
#                         mas relevante para análise evolutiva)
SERINE_MOTIF  = re.compile(r"GDS[AG]G")                    # canônico
ALT_SERINE    = re.compile(r"[GA]DS[GASTVNC]G")            # alt conservado
MUTANT_SER195 = re.compile(r"GD[ACTVNILMF]GG|[GA]D[ACTVNILMF][GA]G")  # Ser195 mutado


def check_completeness(seq_str, seq_id, min_len=220):
    """Verifica critérios de completude de uma tripsina candidata."""
    seq = seq_str.replace("*", "").replace("-", "").upper()

    has_canonical = bool(SERINE_MOTIF.search(seq))
    has_alt       = bool(ALT_SERINE.search(seq))
    has_mutant    = bool(MUTANT_SER195.search(seq)) and not (has_canonical or has_alt)

    # Localizar motivo (prioridade: canônico > alt > mutante)
    m = SERINE_MOTIF.search(seq) or ALT_SERINE.search(seq) or MUTANT_SER195.search(seq)

    result = {
        "id":             seq_id,
        "length":         len(seq),
        "met_start":      seq.startswith("M"),
        "len_ok":         len(seq) >= min_len,
        "gdsgg_motif":    has_canonical,
        "alt_motif":      has_alt,
        "mutant_ser195":  has_mutant,   # Ser195 possivelmente mutado
        "motif_pos":      m.start() if m else None,
        "motif_seq":      m.group() if m else None,
        "pass_strict":    False,
        "pass_mutant":    False,
        "pass_borderline": False,
    }

    base = result["met_start"] and result["len_ok"]

    # Strict: Met + comprimento + Ser195 intacto (canônico ou alt)
    result["pass_strict"] = base and (has_canonical or has_alt)

    # Mutant: Met + comprimento + contexto GDXGG mas Ser195 substituído
    result["pass_mutant"] = base and has_mutant

    # Borderline: Met + comprimento mas sem nenhum motivo
    result["pass_borderline"] = base and not (has_canonical or has_alt or has_mutant)

    return result


def run_triad_check(fasta_input, ref_fasta, out_dir, env="analysis"):
    """Executa verificação da tríade via MAFFT + extract_catalytic_triad.py."""
    bin_dir = Path("nextflow/bin")
    triad_script = bin_dir / "extract_catalytic_triad.py"

    if not triad_script.exists():
        return None

    if not Path(ref_fasta).exists():
        print(f"  ⚠️  Referência bovina não encontrada: {ref_fasta}")
        print(f"     Baixar: wget 'https://www.uniprot.org/uniprot/P00760.fasta' -O {ref_fasta}")
        return None

    print("  Alinhando com tripsina bovina (MAFFT) para verificar tríade...")

    # Combinar ref + candidatos
    combined = out_dir / "combined_for_triad.fasta"
    with open(combined, "w") as f:
        f.write(Path(ref_fasta).read_text())
        f.write("\n")
        f.write(Path(fasta_input).read_text())

    # MAFFT
    aligned = out_dir / "aligned_for_triad.fasta"
    try:
        result = subprocess.run(
            ["mamba", "run", "--no-capture-output", "-n", "discovery",
             "mafft", "--auto", "--quiet", str(combined)],
            capture_output=True, text=True, timeout=300
        )
        if result.returncode != 0:
            print(f"  ⚠️  MAFFT falhou: {result.stderr[:200]}")
            return None
        aligned.write_text(result.stdout)
    except (subprocess.TimeoutExpired, FileNotFoundError):
        print("  ⚠️  MAFFT não disponível — pulando verificação da tríade")
        return None

    # extract_catalytic_triad.py
    triad_out = out_dir / "triad_check.tsv"
    # Tentar detectar o ID da referência
    ref_id = None
    for rec in SeqIO.parse(ref_fasta, "fasta"):
        ref_id = rec.id
        break

    if not ref_id:
        return None

    try:
        subprocess.run([
            sys.executable, str(triad_script),
            "--alignment", str(aligned),
            "--reference", ref_id,
            "--output",    str(triad_out),
        ], check=True, timeout=120)
        return triad_out
    except Exception as e:
        print(f"  ⚠️  Verificação da tríade falhou: {e}")
        return None


def main():
    parser = argparse.ArgumentParser(description="Filtro de completude de tripsinas")
    parser.add_argument("--input",     default="results/phase3/trypsins_confident.fasta")
    parser.add_argument("--also-suggestive", action="store_true",
                        help="Incluir também suggestive na análise")
    parser.add_argument("--results4",  default="results/phase4")
    parser.add_argument("--min-len",   type=int, default=220)
    parser.add_argument("--ref-trypsin", help="FASTA da tripsina bovina para check da tríade")
    args = parser.parse_args()

    out_dir = Path(args.results4)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("\n══════════════════════════════════════════════")
    print("  FASE 4 — Filtro de Completude")
    print("══════════════════════════════════════════════\n")

    input_fasta = Path(args.input)
    if not input_fasta.exists():
        print(f"❌ {input_fasta} não encontrado — execute Fase 3 primeiro")
        sys.exit(1)

    # Também processar suggestive se pedido
    inputs = [input_fasta]
    if args.also_suggestive:
        suggest = Path("results/phase3/trypsins_suggestive.fasta")
        if suggest.exists():
            inputs.append(suggest)
            print(f"Incluindo suggestive: {suggest}")

    all_strict     = []
    all_mutant     = []   # Ser195 possivelmente mutado
    all_borderline = []
    all_report     = []

    for fasta in inputs:
        source = "confident" if "confident" in fasta.name else "suggestive"
        records = list(SeqIO.parse(fasta, "fasta"))
        print(f"\nAnalisando: {fasta.name} ({len(records)} sequências)")

        for rec in records:
            res = check_completeness(str(rec.seq), rec.id, args.min_len)
            res["source"] = source

            if res["pass_strict"]:
                all_strict.append(rec)
            elif res["pass_mutant"]:
                all_mutant.append(rec)
            elif res["pass_borderline"]:
                all_borderline.append(rec)

            all_report.append(res)

    # ── Estatísticas ───────────────────────────────────────────────────────
    n_input      = len(all_report)
    n_strict     = len(all_strict)
    n_mutant     = len(all_mutant)
    n_borderline = len(all_borderline)
    n_fail       = n_input - n_strict - n_mutant - n_borderline

    print(f"\n── Resultados ────────────────────────────────────")
    print(f"  Sequências analisadas:              {n_input:>6}")
    print(f"  Completas (Ser195 intacto):         {n_strict:>6}  ← usar para modelagem")
    print(f"  Possível mutação Ser195 (GDXGG):    {n_mutant:>6}  ← análise evolutiva")
    print(f"  Borderline (sem motivo ao redor):   {n_borderline:>6}  ← investigar manualmente")
    print(f"  Rejeitadas (sem Met ou curtas):     {n_fail:>6}")

    # Por critério
    n_met = sum(1 for r in all_report if r["met_start"])
    n_len = sum(1 for r in all_report if r["len_ok"])
    n_mot = sum(1 for r in all_report if r["gdsgg_motif"] or r["alt_motif"])
    n_mut = sum(1 for r in all_report if r["mutant_ser195"])
    print(f"\n  Por critério:")
    print(f"    Met inicial:               {n_met:>4}/{n_input}")
    print(f"    ≥ {args.min_len} aa:                  {n_len:>4}/{n_input}")
    print(f"    Motivo GDSGG (canônico):   {n_mot:>4}/{n_input}")
    print(f"    Motivo GDXGG (Ser mutado): {n_mut:>4}/{n_input}")

    # Comprimento médio dos completos
    if all_strict:
        lens = [len(str(r.seq).replace("*","")) for r in all_strict]
        print(f"\n  Comprimento médio (strict): {sum(lens)/len(lens):.0f} aa")
        print(f"  Range: {min(lens)}-{max(lens)} aa")

    # ── Verificação da tríade (opcional) ───────────────────────────────────
    ref_trypsin = args.ref_trypsin
    if not ref_trypsin:
        # Tentar path padrão
        default_ref = "data/references/bovine_trypsin_1TGN.fasta"
        if Path(default_ref).exists():
            ref_trypsin = default_ref

    triad_results = {}
    if ref_trypsin and all_strict:
        strict_temp = out_dir / "_strict_temp.fasta"
        SeqIO.write(all_strict, str(strict_temp), "fasta")
        print(f"\n── Verificação da Tríade Catalítica ─────────────")
        triad_tsv = run_triad_check(strict_temp, ref_trypsin, out_dir)
        if triad_tsv and Path(triad_tsv).exists():
            import csv
            for row in csv.DictReader(open(triad_tsv), delimiter="\t"):
                triad_results[row["id"]] = row.get("triad_complete","False") == "True"
            n_triad_ok = sum(1 for v in triad_results.values() if v)
            print(f"  Tríade His+Asp+Ser confirmada: {n_triad_ok}/{len(triad_results)}")

    # ── Salvar FASTAs ──────────────────────────────────────────────────────
    # Se há dados de tríade, filtrar strict para incluir apenas com tríade
    if triad_results:
        final_trypsins = [r for r in all_strict if triad_results.get(r.id, True)]
    else:
        final_trypsins = all_strict

    SeqIO.write(final_trypsins, str(out_dir / "complete_trypsins.fasta"), "fasta")
    SeqIO.write(all_mutant,     str(out_dir / "mutant_ser195_trypsins.fasta"), "fasta")
    SeqIO.write(all_borderline, str(out_dir / "borderline_trypsins.fasta"), "fasta")

    # ── Relatório TSV ──────────────────────────────────────────────────────
    report_path = out_dir / "completeness_report.tsv"
    with open(report_path, "w") as f:
        f.write("id\tsource\tlength\tmet_start\tlen_ok\tgdsgg_motif\talt_motif\t"
                "mutant_ser195\tmotif_seq\tmotif_pos\tpass_strict\tpass_mutant\t"
                "pass_borderline\ttriad_ok\n")
        for r in all_report:
            triad = triad_results.get(r["id"], "NA")
            f.write("\t".join([
                r["id"], r.get("source","?"),
                str(r["length"]), str(r["met_start"]),
                str(r["len_ok"]), str(r["gdsgg_motif"]),
                str(r["alt_motif"]), str(r.get("mutant_ser195", False)),
                str(r.get("motif_seq","NA")), str(r.get("motif_pos","NA")),
                str(r["pass_strict"]), str(r.get("pass_mutant", False)),
                str(r["pass_borderline"]), str(triad),
            ]) + "\n")

    n_final = len(final_trypsins)
    print(f"\n── Arquivos de output ────────────────────────────")
    print(f"  complete_trypsins.fasta:       {n_final} sequências")
    print(f"  mutant_ser195_trypsins.fasta:  {n_mutant} sequências")
    print(f"  borderline_trypsins.fasta:     {n_borderline} sequências")
    print(f"  completeness_report.tsv")

    # ── Aviso biológico ────────────────────────────────────────────────────
    print(f"\n── Interpretação Biológica ───────────────────────")
    if n_final < 5:
        print(f"  ⚠️  AVISO: Apenas {n_final} tripsinas completas.")
        print(f"     Esperado: 5-30 para Lepidoptera midgut.")
        print(f"     Opções:")
        print(f"       1. Verificar --also-suggestive (python filter.py --also-suggestive)")
        print(f"       2. Reduzir --min-len para 200 (python filter.py --min-len 200)")
        print(f"       3. Verificar qualidade do assembly (Fase 1)")
    elif n_final > 50:
        print(f"  ⚠️  AVISO: {n_final} tripsinas completas (esperado 5-30).")
        print(f"     Possível redundância — verificar CD-HIT-EST na Fase 1.")
    else:
        print(f"  ✅ {n_final} tripsinas completas — excelente para Lepidoptera!")
        print(f"     Range esperado: 5-30 isoformas em midgut de Anticarsia gemmatalis")

    # IDs das 10 melhores para BLAST manual de verificação
    print(f"\n── IDs para verificação manual (BLAST web) ──────")
    for i, rec in enumerate(final_trypsins[:10]):
        seq = str(rec.seq).replace("*","")
        print(f"  {i+1:>2}. {rec.id:<50} ({len(seq)} aa)")

    print(f"\n✅ Fase 4 concluída.")
    print(f"\nPróximo passo:")
    print(f"   python scripts/phase4_completeness/validate.py")
    print(f"\nDepois de validar, Fase 5:")
    print(f"   python scripts/phase5_primary_char/01_protparam.py")


if __name__ == "__main__":
    main()
