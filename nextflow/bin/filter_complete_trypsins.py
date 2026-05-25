#!/usr/bin/env python3
"""
filter_complete_trypsins.py
Filtra tripsinas completas por critérios de completude:
1. Met inicial
2. Comprimento >= min_len
3. IDs presentes na interseção DIAMOND + HMMER (opcional)
4. Tríade catalítica detectada por alinhamento (opcional)

Uso:
    # Modo 1: Filtro simples Met + comprimento
    python3 filter_complete_trypsins.py \
        --input trypsins_confident.fasta \
        --min_len 220 \
        --output complete_trypsins.fasta \
        --report completeness_report.tsv

    # Modo 2: Interseção DIAMOND + HMMER (Fase 3)
    python3 filter_complete_trypsins.py \
        --diamond diamond_trypsin_ids.txt \
        --pep proteome.pep \
        --output_confident trypsins_confident.fasta \
        --output_suggestive trypsins_suggestive.fasta \
        --report identification_report.tsv \
        --ids_confident trypsin_ids_confident.txt
"""

import argparse
import sys
from pathlib import Path

try:
    from Bio import SeqIO
except ImportError:
    print("ERRO: Biopython não instalado. Use: pip install biopython", file=sys.stderr)
    sys.exit(1)


def parse_args():
    parser = argparse.ArgumentParser(description="Filtro de tripsinas completas")

    # Modo 1: Filtro de completude
    parser.add_argument("--input",              help="FASTA de entrada")
    parser.add_argument("--min_len",  type=int, default=220, help="Comprimento mínimo em aa")
    parser.add_argument("--output",             help="FASTA de saída (tripsinas completas)")
    parser.add_argument("--report",             help="Relatório TSV de completude")
    parser.add_argument("--triad",              help="TSV com status da tríade catalítica")
    parser.add_argument("--hmm_coverage",       help="HMMER domtblout para cobertura PF00089")
    parser.add_argument("--min_coverage", type=float, default=0.80)

    # Modo 2: Interseção DIAMOND + HMMER
    parser.add_argument("--diamond",            help="Lista de IDs DIAMOND trypsin hits")
    parser.add_argument("--hmm_ids",            help="Lista de IDs HMMER hits (opcional)")
    parser.add_argument("--pep",                help="FASTA proteoma completo (.pep)")
    parser.add_argument("--output_confident",   help="FASTA confident trypsins")
    parser.add_argument("--output_suggestive",  help="FASTA suggestive trypsins")
    parser.add_argument("--ids_confident",      help="Lista de IDs confident")

    return parser.parse_args()


def mode_intersect(args):
    """Modo 2: Interseção DIAMOND + HMMER para identificação de tripsinas."""
    # Carregar IDs DIAMOND
    diamond_ids = set()
    if args.diamond and Path(args.diamond).exists():
        with open(args.diamond) as f:
            diamond_ids = {line.strip() for line in f if line.strip()}

    # Carregar IDs HMMER (se disponível)
    hmm_ids = set()
    if args.hmm_ids and Path(args.hmm_ids).exists():
        with open(args.hmm_ids) as f:
            hmm_ids = {line.strip() for line in f if line.strip()}

    print(f"DIAMOND hits: {len(diamond_ids)}", file=sys.stderr)
    print(f"HMMER hits: {len(hmm_ids)}", file=sys.stderr)

    # Interseção e diferença
    if hmm_ids:
        confident_ids = diamond_ids & hmm_ids
        suggestive_ids = (diamond_ids | hmm_ids) - confident_ids
    else:
        # Sem HMMER: usar DIAMOND como confident
        confident_ids = diamond_ids
        suggestive_ids = set()

    print(f"Confident (intersecção): {len(confident_ids)}", file=sys.stderr)
    print(f"Suggestive: {len(suggestive_ids)}", file=sys.stderr)

    # Extrair sequências do proteoma
    all_records = {rec.id: rec for rec in SeqIO.parse(args.pep, "fasta")}
    confident_recs = [all_records[id_] for id_ in confident_ids if id_ in all_records]
    suggestive_recs = [all_records[id_] for id_ in suggestive_ids if id_ in all_records]

    # Salvar FASTAs
    SeqIO.write(confident_recs, args.output_confident, "fasta")
    SeqIO.write(suggestive_recs, args.output_suggestive, "fasta")

    # Salvar IDs confident
    if args.ids_confident:
        with open(args.ids_confident, "w") as f:
            for id_ in confident_ids:
                f.write(f"{id_}\n")

    # Relatório
    if args.report:
        with open(args.report, "w") as f:
            f.write("id\tdiamond_hit\thmm_hit\tconfident\n")
            all_ids = confident_ids | suggestive_ids
            for id_ in sorted(all_ids):
                in_diamond = id_ in diamond_ids
                in_hmm = id_ in hmm_ids if hmm_ids else "NA"
                is_confident = id_ in confident_ids
                f.write(f"{id_}\t{in_diamond}\t{in_hmm}\t{is_confident}\n")


def mode_completeness(args):
    """Modo 1: Filtro de completude (Met + comprimento + tríade)."""
    if not args.input:
        print("ERRO: --input obrigatório no modo completude", file=sys.stderr)
        sys.exit(1)

    # Carregar status da tríade (se disponível)
    triad_pass = {}
    if args.triad and Path(args.triad).exists():
        with open(args.triad) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    triad_pass[parts[0]] = parts[1].lower() in ("true", "yes", "1", "pass")

    complete_recs = []
    report_rows = []

    for rec in SeqIO.parse(args.input, "fasta"):
        seq = str(rec.seq).replace("*", "")  # remover stop codon

        # Critério 1: Met inicial
        met_ok = seq.startswith("M")

        # Critério 2: comprimento mínimo
        len_ok = len(seq) >= args.min_len

        # Critério 3: tríade catalítica
        triad_ok = triad_pass.get(rec.id, True)  # True se não há dado (não bloqueia)

        passes = met_ok and len_ok and triad_ok
        report_rows.append((rec.id, met_ok, len(seq), triad_ok, passes))

        if passes:
            complete_recs.append(rec)

    print(f"Tripsinas candidatas: {sum(1 for _ in SeqIO.parse(args.input, 'fasta'))}", file=sys.stderr)
    print(f"Tripsinas completas: {len(complete_recs)}", file=sys.stderr)

    # Salvar FASTA
    if args.output:
        SeqIO.write(complete_recs, args.output, "fasta")

    # Relatório TSV
    if args.report:
        with open(args.report, "w") as f:
            f.write("id\tmet_start\tlength_aa\ttriad_ok\tpass_completeness\n")
            for id_, met, ln, triad, passes in report_rows:
                f.write(f"{id_}\t{met}\t{ln}\t{triad}\t{passes}\n")

    # Aviso biológico
    n = len(complete_recs)
    if n < 5:
        print(f"⚠️  AVISO: Apenas {n} tripsinas completas (esperado 5-30). "
              "Considere reduzir min_len ou verificar qualidade do assembly.", file=sys.stderr)
    elif n > 50:
        print(f"⚠️  AVISO: {n} tripsinas completas (esperado 5-30). "
              "Possível redundância — verifique CD-HIT-EST na Fase 1.", file=sys.stderr)
    else:
        print(f"✅ {n} tripsinas completas — dentro do range esperado (5-30) para Lepidoptera", file=sys.stderr)


def main():
    args = parse_args()

    # Detectar modo
    if args.diamond or args.pep:
        mode_intersect(args)
    else:
        mode_completeness(args)


if __name__ == "__main__":
    main()
