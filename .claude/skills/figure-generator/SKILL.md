---
name: figure-generator
description: Gera figuras científicas (300 DPI, vetorial) para manuscrito usando Matplotlib, PyMOL e Seaborn. Trigger: "gerar figuras", "fazer gráfico", "visualização", "figura para manuscrito", "PyMOL", na Fase 11.
---

# Figure Generator (Fase 11)

## Quando usar
- Após todas as fases analíticas concluídas (Fases 1-10)
- Usuário pede "gerar figuras", "fazer visualização", "figuras para o artigo"
- Chamado pelo `publication-formatter` automaticamente

## Figuras obrigatórias para manuscrito (IF > 5)

### Figura 1 — Visão geral do pipeline
- Diagrama de fluxo das 11 fases
- Ferramenta: `graphviz` ou `matplotlib.patches`

### Figura 2 — Identificação e completude das tripsinas
- Painel A: Diagrama de Venn (DIAMOND ∩ HMMER)
- Painel B: Gráfico de barras — contagem confident/suggestive
- Painel C: Distribuição de comprimento (violin plot) vs. range esperado (220-320 aa)

### Figura 3 — Árvore filogenética
- Saída IQ-TREE2 com suporte bootstrap
- Colorir por espécie (Lepidoptera) com `ETE3` ou `iTOL`
- Escala de tempo (se ultrametric tree com r8s)

### Figura 4 — Estruturas 3D (PyMOL)
- Representação ribbon colorida por pLDDT (azul=alto, vermelho=baixo)
- Sítio ativo com stick e superfície de bolso S1
- Sobreposição com tripsina bovina (1TGN) — RMSD no label

### Figura 5 — Caracterização primária
- Heatmap: peso molecular × pI × instabilidade × GRAVY
- Gráfico: distribuição de pLDDT dos modelos AF3

### Figura 6 — Docking
- Melhor pose de cada inibidor (SKTI, BPTI-TGPCK, BPTI-AVIMK)
- Superfície do sítio S1 com ligante em stick
- Tabela de scores Vina/HADDOCK no painel

### Figura 7 — Dinâmica Molecular
- RMSD ao longo do tempo (100 ns × 3 réplicas)
- RMSF por resíduo (destacar loop de ativação e sítio ativo)
- Raio de giro (Rg) — indicador de estabilidade global

## Procedimento

### Setup (matplotlib + seaborn)
```python
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Estilo padrão para publicação
plt.rcParams.update({
    'figure.dpi': 300,
    'font.family': 'Arial',
    'font.size': 10,
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'legend.fontsize': 9,
    'savefig.bbox': 'tight',
    'savefig.format': 'svg',   # vetorial para edição
})
```

### Gerar Figura 4 — PyMOL script
```python
# pymol_figure4.py
from pymol import cmd

cmd.load("tripsina_best.pdb", "tripsina")
cmd.load("1TGN.pdb", "bovine")

# Colorir por B-factor (pLDDT no AF3)
cmd.spectrum("b", "red_white_blue", "tripsina")

# Sítio ativo
cmd.select("catalytic_triad", "resi 57+102+195 and tripsina")
cmd.show("sticks", "catalytic_triad")
cmd.color("yellow", "catalytic_triad")

# Sobrepor com tripsina bovina
cmd.align("tripsina", "bovine")

# Exportar
cmd.ray(2400, 2400)
cmd.png("figure4_structures.png", dpi=300)
cmd.save("figure4_structures.svg")
```

### Gerar Figura 7 — MD plots
```python
import MDAnalysis as mda
import matplotlib.pyplot as plt

# RMSD por réplica
fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
colors = ['#2196F3', '#F44336', '#4CAF50']

for i, (rep, color) in enumerate(zip(['rep1','rep2','rep3'], colors)):
    data = np.loadtxt(f"{rep}_rmsd.xvg", comments=['@','#'])
    time_ns = data[:,0] / 1000  # ps → ns
    axes[0].plot(time_ns, data[:,1]*10, color=color, alpha=0.8, label=f"Réplica {i+1}")

axes[0].set_ylabel("RMSD (Å)")
axes[0].legend()
axes[0].set_title("Estabilidade do backbone (Cα RMSD)")
plt.tight_layout()
plt.savefig("figure7_md.svg", dpi=300)
```

## Outputs garantidos
- `manuscript/figures/figure1_pipeline.svg`
- `manuscript/figures/figure2_identification.svg`
- `manuscript/figures/figure3_phylogeny.svg`
- `manuscript/figures/figure4_structures.png` (300 DPI)
- `manuscript/figures/figure5_primary_char.svg`
- `manuscript/figures/figure6_docking.png` (300 DPI)
- `manuscript/figures/figure7_md.svg`

## Padrões das revistas alvo
| Revista | Max figures | Formato | Max size |
|---|---|---|---|
| Insect Biochem Mol Biol | 8 | TIFF/EPS/PDF | 10 MB |
| Int J Biol Macromol | 10 | PDF/SVG/EPS | 15 MB |
| Comput Struct Biotechnol J | 10 | SVG/PDF/PNG 300 DPI | Sem limite |

## Erros comuns e correções
| Erro | Causa | Fix |
|---|---|---|
| Figura borrada na impressão | PNG com DPI < 300 | Sempre use SVG ou PNG 300 DPI |
| PyMOL fecha sem salvar | GUI mode ativo | Use `pymol -c` (modo linha de comando) |
| Fonte não incorporada no PDF | Arial não embutida | Use `matplotlib.backends.backend_pdf` |
| Cores diferem no print | RGB vs CMYK | Converta para CMYK no Inkscape antes de submeter |
