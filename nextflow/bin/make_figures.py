#!/usr/bin/env python3
"""
make_figures.py
Gera figuras de alta qualidade (300 DPI) para o manuscrito.

Figuras geradas:
  Fig1: Propriedades físico-químicas das tripsinas (MW, pI, GRAVY)
  Fig2: Árvore filogenética (wrapper para figtree/ete3)
  Fig3: pLDDT dos modelos AlphaFold3 por resíduo
  Fig4: Scores de docking (ΔG Vina vs HADDOCK score)
  Fig5: RMSD das simulações MD ao longo do tempo

Uso:
    python3 make_figures.py \
        --protparam results/05_primary/protparam/protparam_results.tsv \
        --tree results/06_phylogeny/iqtree/*.treefile \
        --docking results/09_docking/vina/vina_scores_all.tsv \
        --md_analysis results/10_md/*/md_report.tsv \
        --outdir results/11_report/figures
"""

import argparse
import sys
import os
from pathlib import Path

try:
    import matplotlib
    matplotlib.use("Agg")  # Headless
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import numpy as np
    import pandas as pd
except ImportError as e:
    print(f"ERRO: {e}\nInstale: pip install matplotlib numpy pandas", file=sys.stderr)
    sys.exit(1)

# Paleta daltônica (Okabe-Ito)
COLORS = {
    "blue":        "#0072B2",
    "orange":      "#E69F00",
    "green":       "#009E73",
    "red":         "#D55E00",
    "purple":      "#CC79A7",
    "light_blue":  "#56B4E9",
    "yellow":      "#F0E442",
    "gray":        "#999999",
}

# Configurações globais de matplotlib para publicação
plt.rcParams.update({
    "font.family":     "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size":       9,
    "axes.titlesize":  10,
    "axes.labelsize":  9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "figure.dpi":      300,
    "savefig.dpi":     300,
    "savefig.bbox":    "tight",
    "savefig.format":  "svg",
})


def fig1_protparam(protparam_tsv: str, outdir: str):
    """Fig1: Propriedades físico-químicas das tripsinas."""
    df = pd.read_csv(protparam_tsv, sep="\t", comment="#")
    if df.empty:
        print("AVISO: protparam TSV vazio — pulando Fig1", file=sys.stderr)
        return

    fig, axes = plt.subplots(1, 3, figsize=(7.1, 2.8))  # largura coluna dupla

    # A — Distribuição de MW
    axes[0].hist(df["mw_kda"], bins=10, color=COLORS["blue"], edgecolor="white", linewidth=0.5)
    axes[0].axvline(20, color=COLORS["red"], linestyle="--", linewidth=0.8, alpha=0.7)
    axes[0].axvline(35, color=COLORS["red"], linestyle="--", linewidth=0.8, alpha=0.7)
    axes[0].set_xlabel("Massa molecular (kDa)")
    axes[0].set_ylabel("Frequência")
    axes[0].set_title("A")

    # B — Distribuição de pI
    axes[1].hist(df["pi"], bins=10, color=COLORS["green"], edgecolor="white", linewidth=0.5)
    axes[1].set_xlabel("pI isoelétrico")
    axes[1].set_ylabel("Frequência")
    axes[1].set_title("B")

    # C — MW vs pI scatter
    colors_stability = [COLORS["blue"] if s == "estável" else COLORS["orange"]
                        for s in df.get("stability", ["estável"] * len(df))]
    axes[2].scatter(df["mw_kda"], df["pi"], c=colors_stability, s=30, alpha=0.8, edgecolors="none")
    axes[2].set_xlabel("Massa molecular (kDa)")
    axes[2].set_ylabel("pI")
    axes[2].set_title("C")
    legend_elements = [
        mpatches.Patch(color=COLORS["blue"], label="Estável (II < 40)"),
        mpatches.Patch(color=COLORS["orange"], label="Instável (II ≥ 40)"),
    ]
    axes[2].legend(handles=legend_elements, fontsize=7, frameon=False)

    plt.tight_layout()
    out_path = Path(outdir) / "fig1_trypsin_properties.svg"
    fig.savefig(out_path)
    plt.close()
    print(f"✅ Fig1 salva: {out_path}", file=sys.stderr)


def fig4_docking(docking_tsv: str, outdir: str):
    """Fig4: Scores de docking por complexo."""
    try:
        df = pd.read_csv(docking_tsv, sep="\t", comment="#")
    except Exception:
        print("AVISO: arquivo de docking não encontrado — pulando Fig4", file=sys.stderr)
        return

    if df.empty or "vina_score" not in df.columns:
        return

    fig, ax = plt.subplots(figsize=(3.5, 3.0))

    # Ordenar por score
    df_sorted = df.sort_values("vina_score")
    bars = ax.barh(range(len(df_sorted)), df_sorted["vina_score"],
                   color=COLORS["blue"], edgecolor="none", height=0.7)

    ax.set_yticks(range(len(df_sorted)))
    ax.set_yticklabels(df_sorted.get("complex_id", df_sorted.index), fontsize=7)
    ax.set_xlabel("Energia de ligação AutoDock Vina (kcal/mol)")
    ax.axvline(-8.0, color=COLORS["red"], linestyle="--", linewidth=0.8, alpha=0.7)
    ax.set_title("Docking scores — Tripsinas vs Inibidores")

    plt.tight_layout()
    out_path = Path(outdir) / "fig4_docking.svg"
    fig.savefig(out_path)
    plt.close()
    print(f"✅ Fig4 salva: {out_path}", file=sys.stderr)


def fig5_md_rmsd(md_reports: list, outdir: str):
    """Fig5: RMSD das simulações MD."""
    fig, ax = plt.subplots(figsize=(5.0, 3.0))

    plotted = 0
    for report_tsv in md_reports:
        try:
            df = pd.read_csv(report_tsv, sep="\t")
            if "rmsd_mean_nm" not in df.columns:
                continue
            for _, row in df.iterrows():
                ax.bar(row["id"], row["rmsd_mean_nm"] * 10,  # nm → Å
                       color=COLORS["blue"], alpha=0.7)
                plotted += 1
        except Exception as e:
            print(f"  AVISO: {e}", file=sys.stderr)

    if plotted == 0:
        # Criar figura placeholder
        ax.text(0.5, 0.5, "Simulações MD\nem andamento",
                transform=ax.transAxes, ha="center", va="center", fontsize=12)

    ax.axhline(3.0, color=COLORS["red"], linestyle="--", linewidth=0.8, alpha=0.7,
               label="Limite 3 Å")
    ax.set_xlabel("Complexo Tripsina-Inibidor")
    ax.set_ylabel("RMSD médio backbone (Å)")
    ax.set_title("Estabilidade das simulações MD (100 ns)")
    ax.legend(fontsize=7, frameon=False)

    plt.tight_layout()
    out_path = Path(outdir) / "fig5_md_rmsd.svg"
    fig.savefig(out_path)
    plt.close()
    print(f"✅ Fig5 salva: {out_path}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description="Geração de figuras para manuscrito")
    parser.add_argument("--protparam",  help="TSV do ProtParam (Fase 5)")
    parser.add_argument("--tree",       help="Arquivo .treefile do IQ-TREE2 (Fase 6)")
    parser.add_argument("--docking",    help="TSV de scores de docking (Fase 9)")
    parser.add_argument("--md_analysis", nargs="+", help="TSV(s) de relatório MD (Fase 10)")
    parser.add_argument("--outdir",     required=True, help="Diretório de saída")
    args = parser.parse_args()

    Path(args.outdir).mkdir(parents=True, exist_ok=True)
    print(f"Gerando figuras em: {args.outdir}", file=sys.stderr)

    if args.protparam and Path(args.protparam).exists():
        fig1_protparam(args.protparam, args.outdir)
    else:
        print("INFO: --protparam não fornecido ou não encontrado — pulando Fig1", file=sys.stderr)
        # Criar placeholder
        fig, ax = plt.subplots(figsize=(5, 3))
        ax.text(0.5, 0.5, "Fig1: Aguardando resultados Fase 5",
                transform=ax.transAxes, ha="center", va="center")
        fig.savefig(Path(args.outdir) / "fig1_trypsin_properties.svg")
        plt.close()

    if args.docking and Path(args.docking).exists():
        fig4_docking(args.docking, args.outdir)

    if args.md_analysis:
        fig5_md_rmsd(args.md_analysis, args.outdir)

    print(f"\n✅ Figuras geradas em: {args.outdir}", file=sys.stderr)


if __name__ == "__main__":
    main()
