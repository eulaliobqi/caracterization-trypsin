#!/usr/bin/env python3
"""
compute_protparam.py
Calcula propriedades físico-químicas de proteínas usando BioPython ProtParam.
Equivalente ao ExPASy ProtParam para análise em batch.

Propriedades calculadas:
- Peso molecular (kDa)
- pI isoelétrico
- Índice de instabilidade
- GRAVY (Grand average of hydropathicity)
- Composição de aminoácidos
- Fórmula molecular
- Número de resíduos carregados

Uso:
    python3 compute_protparam.py \
        --input complete_trypsins.fasta \
        --output protparam_results.tsv
"""

import argparse
import sys

try:
    from Bio import SeqIO
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
except ImportError:
    print("ERRO: Biopython necessário. pip install biopython", file=sys.stderr)
    sys.exit(1)


def calculate_protparam(seq_str: str) -> dict:
    """Calcula todas as propriedades físico-químicas de uma sequência proteica."""
    # Remover stop codon e gaps
    seq_clean = seq_str.replace("*", "").replace("-", "").replace("X", "")

    if not seq_clean:
        return None

    try:
        analysis = ProteinAnalysis(seq_clean)

        # Propriedades básicas
        mw_kda = analysis.molecular_weight() / 1000
        pi = analysis.isoelectric_point()
        instability = analysis.instability_index()
        gravy = analysis.gravy()
        aa_count = len(seq_clean)

        # Composição
        aa_comp = analysis.get_amino_acids_percent()

        # Resíduos carregados
        neg_charged = aa_comp.get("D", 0) + aa_comp.get("E", 0)
        pos_charged = aa_comp.get("K", 0) + aa_comp.get("R", 0)
        neg_n = round(aa_comp.get("D", 0) * aa_count) + round(aa_comp.get("E", 0) * aa_count)
        pos_n = round(aa_comp.get("K", 0) * aa_count) + round(aa_comp.get("R", 0) * aa_count)

        # Classificação
        stability = "estável" if instability < 40 else "instável"
        hydrophobicity = "hidrofóbica" if gravy > 0 else "hidrofílica"

        return {
            "mw_kda": round(mw_kda, 2),
            "pi": round(pi, 2),
            "instability_index": round(instability, 2),
            "stability": stability,
            "gravy": round(gravy, 4),
            "hydrophobicity": hydrophobicity,
            "aa_count": aa_count,
            "neg_charged_residues": neg_n,
            "pos_charged_residues": pos_n,
            "pct_hydrophobic": round(
                (aa_comp.get("A", 0) + aa_comp.get("I", 0) + aa_comp.get("L", 0) +
                 aa_comp.get("M", 0) + aa_comp.get("F", 0) + aa_comp.get("V", 0) +
                 aa_comp.get("W", 0) + aa_comp.get("P", 0)) * 100, 1
            ),
            "pct_ala": round(aa_comp.get("A", 0) * 100, 1),
            "pct_cys": round(aa_comp.get("C", 0) * 100, 1),
        }
    except Exception as e:
        print(f"  AVISO: Erro ao calcular ProtParam: {e}", file=sys.stderr)
        return None


def main():
    parser = argparse.ArgumentParser(description="Cálculo de propriedades físico-químicas (ProtParam)")
    parser.add_argument("--input", required=True, help="FASTA de proteínas")
    parser.add_argument("--output", required=True, help="TSV de saída")
    args = parser.parse_args()

    records = list(SeqIO.parse(args.input, "fasta"))
    print(f"Sequências para análise: {len(records)}", file=sys.stderr)

    # Avaliar range esperado para tripsinas de Lepidoptera
    expected_mw_min = 20.0
    expected_mw_max = 40.0
    expected_pi_min = 3.0
    expected_pi_max = 7.0

    with open(args.output, "w") as f:
        # Header
        f.write("\t".join([
            "id", "mw_kda", "pi", "instability_index", "stability",
            "gravy", "hydrophobicity", "aa_count",
            "neg_charged", "pos_charged", "pct_hydrophobic",
            "mw_in_range", "pi_in_range", "notes"
        ]) + "\n")

        n_ok = 0
        for rec in records:
            props = calculate_protparam(str(rec.seq))
            if props is None:
                f.write(f"{rec.id}\tERROR\n")
                continue

            # Validação biológica
            mw_ok = expected_mw_min <= props["mw_kda"] <= expected_mw_max
            pi_ok = expected_pi_min <= props["pi"] <= expected_pi_max
            notes = []
            if not mw_ok:
                notes.append(f"MW fora do range ({props['mw_kda']:.1f} kDa)")
            if not pi_ok:
                notes.append(f"pI fora do range ({props['pi']:.1f})")
            if mw_ok and pi_ok:
                n_ok += 1

            f.write("\t".join([
                rec.id,
                str(props["mw_kda"]),
                str(props["pi"]),
                str(props["instability_index"]),
                props["stability"],
                str(props["gravy"]),
                props["hydrophobicity"],
                str(props["aa_count"]),
                str(props["neg_charged_residues"]),
                str(props["pos_charged_residues"]),
                str(props["pct_hydrophobic"]),
                str(mw_ok),
                str(pi_ok),
                "; ".join(notes) if notes else "OK"
            ]) + "\n")

    print(f"✅ ProtParam concluído: {n_ok}/{len(records)} dentro do range esperado para tripsinas", file=sys.stderr)
    print(f"   Range MW: {expected_mw_min}-{expected_mw_max} kDa", file=sys.stderr)
    print(f"   Range pI: {expected_pi_min}-{expected_pi_max}", file=sys.stderr)


if __name__ == "__main__":
    main()
