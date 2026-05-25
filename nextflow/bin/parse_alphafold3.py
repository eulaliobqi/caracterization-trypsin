#!/usr/bin/env python3
"""
parse_alphafold3.py
Converte FASTA para formato JSON de input do AlphaFold3.
Também analisa outputs de confiança do AF3 após a predição.

Uso:
    # Converter FASTA para JSON (antes do AF3)
    python3 parse_alphafold3.py \
        --fasta trypsin.fasta \
        --id TRINITY_DN1_trypsin \
        --output trypsin_af3_input.json

    # Analisar output de confiança (após AF3)
    python3 parse_alphafold3.py \
        --confidence results/trypsin_confidence.json \
        --pdb results/trypsin_model_0.pdb \
        --report plddt_summary.tsv
"""

import argparse
import json
import sys
from pathlib import Path

try:
    from Bio import SeqIO
except ImportError:
    print("ERRO: Biopython necessário. pip install biopython", file=sys.stderr)
    sys.exit(1)


def fasta_to_af3_json(fasta_path: str, seq_id: str, output_path: str):
    """Converte FASTA para formato JSON do AlphaFold3."""
    records = list(SeqIO.parse(fasta_path, "fasta"))
    if not records:
        print(f"ERRO: Nenhuma sequência encontrada em {fasta_path}", file=sys.stderr)
        sys.exit(1)

    rec = records[0]
    seq_str = str(rec.seq).replace("*", "")  # remover stop codon

    af3_input = {
        "name": seq_id,
        "sequences": [
            {
                "proteinChain": {
                    "sequence": seq_str,
                    "count": 1
                }
            }
        ],
        "dialect": "alphafold3",
        "version": 1
    }

    with open(output_path, "w") as f:
        json.dump(af3_input, f, indent=2)

    print(f"JSON AF3 criado: {output_path}", file=sys.stderr)
    print(f"  ID: {seq_id}, Comprimento: {len(seq_str)} aa", file=sys.stderr)


def analyze_confidence(confidence_path: str, pdb_path: str = None, report_path: str = None):
    """Analisa arquivo de confiança do AF3 e gera sumário."""
    with open(confidence_path) as f:
        conf = json.load(f)

    # Extrair métricas
    plddts = conf.get("atom_plddts", conf.get("plddt", []))
    ptm = conf.get("ptm", conf.get("ranking_score", None))
    iptm = conf.get("iptm", None)

    if not plddts:
        print("AVISO: Nenhum dado de pLDDT encontrado no arquivo de confiança", file=sys.stderr)
        return

    plddt_mean = sum(plddts) / len(plddts)
    plddt_min  = min(plddts)
    plddt_max  = max(plddts)

    # Classificação de qualidade
    if plddt_mean >= 90:
        quality = "muito_alta"
    elif plddt_mean >= 70:
        quality = "alta"
    elif plddt_mean >= 50:
        quality = "media"
    else:
        quality = "baixa"

    # Percentagem por categoria
    n = len(plddts)
    pct_very_high = sum(1 for p in plddts if p >= 90) / n * 100
    pct_high      = sum(1 for p in plddts if 70 <= p < 90) / n * 100
    pct_medium    = sum(1 for p in plddts if 50 <= p < 70) / n * 100
    pct_low       = sum(1 for p in plddts if p < 50) / n * 100

    # Nome do modelo (a partir do path do confiança)
    model_id = Path(confidence_path).stem.replace("_confidence", "")

    print(f"\n📊 pLDDT para {model_id}:", file=sys.stderr)
    print(f"  Médio: {plddt_mean:.1f} ({quality})", file=sys.stderr)
    print(f"  Min: {plddt_min:.1f}, Max: {plddt_max:.1f}", file=sys.stderr)
    print(f"  Muito alta (≥90): {pct_very_high:.1f}%", file=sys.stderr)
    print(f"  Alta (70-90): {pct_high:.1f}%", file=sys.stderr)
    print(f"  Média (50-70): {pct_medium:.1f}%", file=sys.stderr)
    print(f"  Baixa (<50): {pct_low:.1f}%  (regiões desordenadas)", file=sys.stderr)
    if ptm:
        print(f"  pTM: {ptm:.3f}", file=sys.stderr)

    # Validação biológica
    if plddt_mean < 70:
        print(f"  ⚠️  AVISO: pLDDT médio {plddt_mean:.1f} < 70 — verificar manualmente", file=sys.stderr)
    elif plddt_mean < 80:
        print(f"  ⚠️  INFO: pLDDT médio {plddt_mean:.1f} (aceitável, mas não ideal)", file=sys.stderr)
    else:
        print(f"  ✅ pLDDT > 80 — alta confiança na estrutura", file=sys.stderr)

    # Salvar relatório
    if report_path:
        with open(report_path, "w") as f:
            f.write("id\tplddt_mean\tplddt_min\tplddt_max\tptm\t"
                    "pct_very_high\tpct_high\tpct_medium\tpct_low\tquality\n")
            f.write(f"{model_id}\t{plddt_mean:.2f}\t{plddt_min:.2f}\t{plddt_max:.2f}\t"
                    f"{ptm if ptm else 'NA'}\t"
                    f"{pct_very_high:.1f}\t{pct_high:.1f}\t{pct_medium:.1f}\t{pct_low:.1f}\t"
                    f"{quality}\n")
        print(f"Relatório salvo: {report_path}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description="AlphaFold3 helper — conversão e análise")
    # Modo input
    parser.add_argument("--fasta",   help="FASTA de entrada para conversão")
    parser.add_argument("--id",      help="ID para o JSON AF3")
    parser.add_argument("--output",  help="JSON de saída para AF3")
    # Modo análise
    parser.add_argument("--confidence", help="JSON de confiança do AF3")
    parser.add_argument("--pdb",        help="PDB do modelo (opcional)")
    parser.add_argument("--report",     help="TSV de sumário de confiança")

    args = parser.parse_args()

    if args.fasta and args.output:
        seq_id = args.id or Path(args.fasta).stem
        fasta_to_af3_json(args.fasta, seq_id, args.output)
    elif args.confidence:
        analyze_confidence(args.confidence, args.pdb, args.report)
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
