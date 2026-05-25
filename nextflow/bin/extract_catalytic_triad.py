#!/usr/bin/env python3
"""
extract_catalytic_triad.py
Detecta a tríade catalítica His57-Asp102-Ser195 por alinhamento múltiplo.
Alinha cada tripsina candidata com bovine trypsin (referência) e verifica
se as posições conservadas contêm os resíduos catalíticos corretos.

Numeração: His57, Asp102, Ser195 (quimotripsina) ≡ Tryp_SPc ativa

Uso:
    python3 extract_catalytic_triad.py \
        --alignment aligned_with_ref.fasta \
        --reference "1TGN_A" \
        --output triad_check.tsv
"""

import argparse
import sys
from collections import defaultdict

try:
    from Bio import SeqIO, AlignIO
    from Bio.Align import MultipleSeqAlignment
except ImportError:
    print("ERRO: Biopython >= 1.80 necessário. pip install biopython", file=sys.stderr)
    sys.exit(1)

# Posições da tríade catalítica na numeração da quimotripsina bovina
# (serão mapeadas via alinhamento para colunas do MSA)
CATALYTIC_RESIDUES = {
    "His57":  {"aa": "H", "pos_chymotrypsin": 57},
    "Asp102": {"aa": "D", "pos_chymotrypsin": 102},
    "Ser195": {"aa": "S", "pos_chymotrypsin": 195},
}


def find_ref_columns(alignment: MultipleSeqAlignment, ref_id: str) -> dict:
    """
    Encontra as colunas do MSA que correspondem às posições 57, 102, 195
    da sequência de referência (bovine trypsin).
    """
    ref_seq = None
    for rec in alignment:
        if rec.id == ref_id or rec.id.startswith(ref_id):
            ref_seq = str(rec.seq)
            break

    if ref_seq is None:
        raise ValueError(
            f"Sequência de referência '{ref_id}' não encontrada no alinhamento. "
            f"IDs disponíveis: {[r.id for r in alignment][:5]}"
        )

    # Mapear posição no alinhamento → posição na sequência (sem gaps)
    col_to_seqpos = {}
    seq_pos = 0
    for col_idx, aa in enumerate(ref_seq):
        if aa != "-":
            seq_pos += 1
            col_to_seqpos[seq_pos] = col_idx

    # Obter colunas das posições catalíticas
    catalytic_cols = {}
    for name, info in CATALYTIC_RESIDUES.items():
        pos = info["pos_chymotrypsin"]
        if pos in col_to_seqpos:
            catalytic_cols[name] = col_to_seqpos[pos]
            print(f"  Coluna MSA para {name}: {col_to_seqpos[pos]+1}", file=sys.stderr)
        else:
            print(f"  ⚠️  {name} (pos {pos}) não encontrada na referência", file=sys.stderr)

    return catalytic_cols


def check_triad(alignment: MultipleSeqAlignment, ref_id: str, catalytic_cols: dict) -> dict:
    """
    Para cada sequência no alinhamento (exceto referência),
    verifica se as posições catalíticas têm os aminoácidos corretos.
    """
    results = {}
    for rec in alignment:
        if rec.id == ref_id or rec.id.startswith(ref_id):
            continue

        seq = str(rec.seq)
        triad_status = {}

        for name, col_idx in catalytic_cols.items():
            expected_aa = CATALYTIC_RESIDUES[name]["aa"]
            actual_aa = seq[col_idx] if col_idx < len(seq) else "-"
            is_present = actual_aa.upper() == expected_aa
            triad_status[name] = {
                "expected": expected_aa,
                "found": actual_aa,
                "ok": is_present,
            }

        # Triade completa = todos os 3 resíduos corretos
        all_present = all(v["ok"] for v in triad_status.values())
        results[rec.id] = {
            "triad_complete": all_present,
            "details": triad_status,
        }

    return results


def main():
    parser = argparse.ArgumentParser(description="Detecção da tríade catalítica de tripsinas")
    parser.add_argument("--alignment", required=True, help="MSA FASTA (com sequência de referência)")
    parser.add_argument("--reference", required=True, help="ID da sequência de referência no MSA")
    parser.add_argument("--output", required=True, help="Arquivo TSV de saída")
    args = parser.parse_args()

    print(f"Lendo alinhamento: {args.alignment}", file=sys.stderr)
    try:
        alignment = AlignIO.read(args.alignment, "fasta")
    except Exception as e:
        print(f"ERRO ao ler alinhamento: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"Sequências no alinhamento: {len(alignment)}", file=sys.stderr)
    print(f"Comprimento do alinhamento: {alignment.get_alignment_length()} colunas", file=sys.stderr)

    # Encontrar colunas catalíticas
    print("Localizando posições catalíticas na referência...", file=sys.stderr)
    try:
        catalytic_cols = find_ref_columns(alignment, args.reference)
    except ValueError as e:
        print(f"ERRO: {e}", file=sys.stderr)
        sys.exit(1)

    if not catalytic_cols:
        print("ERRO: Nenhuma posição catalítica encontrada", file=sys.stderr)
        sys.exit(1)

    # Verificar tríade em todas as sequências
    results = check_triad(alignment, args.reference, catalytic_cols)

    # Estatísticas
    n_complete = sum(1 for v in results.values() if v["triad_complete"])
    n_total = len(results)
    print(f"\nResultados:", file=sys.stderr)
    print(f"  Tríade completa (His+Asp+Ser): {n_complete}/{n_total}", file=sys.stderr)
    print(f"  Tríade incompleta: {n_total - n_complete}", file=sys.stderr)

    # Escrever TSV
    with open(args.output, "w") as f:
        # Header
        f.write("id\ttriad_complete\t")
        for name in CATALYTIC_RESIDUES:
            f.write(f"{name}_found\t{name}_ok\t")
        f.write("\n")

        # Dados
        for seq_id, data in sorted(results.items()):
            row = [seq_id, str(data["triad_complete"])]
            for name in CATALYTIC_RESIDUES:
                if name in data["details"]:
                    d = data["details"][name]
                    row.extend([d["found"], str(d["ok"])])
                else:
                    row.extend(["NA", "False"])
            f.write("\t".join(row) + "\n")

    print(f"\nResultados salvos em: {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
