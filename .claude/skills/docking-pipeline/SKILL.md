---
name: docking-pipeline
description: Executa docking molecular de tripsinas com inibidores SKTI/BPTI via AutoDock Vina + HADDOCK + análise PLIP. Trigger: "docking", "interação tripsina-inibidor", "AutoDock", "HADDOCK", "SKTI", "BPTI", após validação estrutural (Fase 9).
---

# Docking Pipeline (Fase 9)

## Quando usar
- Após `STRUCT_VALIDATION` produzir PDBs validados
- Usuário pede "docking", "afinidade de inibidor", "sítio de ligação"
- Para testar hipóteses de resistência a inibidores de soja (SKTI, BPTI)

## Inputs esperados
- `tripsina_validated.pdb` — estrutura AF3 validada (pLDDT > 80, Ramachandran > 95%)
- `skti.pdb` — Soybean Kunitz Trypsin Inhibitor (PDB: 1AVU cadeia I)
- `bpti_tgpck.pdb` / `bpti_avimk.pdb` — peptídeos BPTI-RCL (preparados de novo)

## Inibidores alvo do projeto
| Inibidor | PDB source | Tipo | Relevância |
|---|---|---|---|
| SKTI (Kunitz) | PDB 1AVU | Proteína completa (~18 kDa) | Principal inibidor da soja |
| BPTI-RCL (TGPCK) | Sintético | Peptídeo loop ativo | Variante resistente |
| BPTI-RCL (AVIMK) | Sintético | Peptídeo loop ativo | Variante alternativa |

## Procedimento

### Etapa 1 — Preparação das estruturas (AutoDockTools)
```bash
# Preparar receptor (tripsina)
prepare_receptor \
    -r tripsina.pdb \
    -o tripsina.pdbqt \
    -A hydrogens \
    -U nphs_lps_waters_deleteAltB

# Preparar ligante (inibidor pequeno)
# Para proteína grande como SKTI: usar HADDOCK (passo 2)
prepare_ligand \
    -l bpti_peptide.pdb \
    -o bpti_peptide.pdbqt
```

### Etapa 2 — Definir caixa de docking (sítio ativo da tripsina)
```python
# Centro do sítio ativo: resíduos Ser195, His57, Asp102
# (numeração quimotripsina — usar extract_catalytic_triad.py para obter coords)
center_x, center_y, center_z = get_catalytic_center("tripsina.pdb")

# Caixa: 25x25x25 Å ao redor do sítio catalítico
box_size = 25  # Å
```

### Etapa 3 — AutoDock Vina (peptídeos BPTI-RCL)
```bash
vina \
    --receptor tripsina.pdbqt \
    --ligand bpti_peptide.pdbqt \
    --center_x ${center_x} \
    --center_y ${center_y} \
    --center_z ${center_z} \
    --size_x 25 --size_y 25 --size_z 25 \
    --exhaustiveness 32 \
    --num_modes 10 \
    --energy_range 3 \
    --out vina_poses.pdbqt \
    --log vina_results.log
```

### Etapa 4 — HADDOCK (docking proteína-proteína para SKTI)
```bash
# HADDOCK requer definição de resíduos de interface (AIRs)
# Interface da tripsina: loop de ativação + sítio S1 (Asp189 bolso especificidade)
# Interface do SKTI: loop reativo (P1-Arg63, P2-Thr64)
haddock3 haddock.cfg \
    --receptor tripsina.pdb \
    --ligand skti.pdb \
    --ambig_fname AIRs_trypsin_skti.tbl \
    --run_dir haddock_run
```

### Etapa 5 — Análise PLIP (interações proteína-ligante)
```bash
# PLIP: detecta H-bonds, hidrofóbicos, sal bridges, Pi-stacking
plip -f melhor_complexo.pdb \
     -o plip_report \
     --xml --pymol
```

## Outputs garantidos
- `results/09_docking/<meta.id>/vina_poses.pdbqt` — poses AutoDock Vina
- `results/09_docking/<meta.id>/haddock_run/` — run HADDOCK (SKTI)
- `results/09_docking/<meta.id>/plip_report.xml` — interações moleculares
- `results/09_docking/<meta.id>/docking_summary.tsv` — scores consolidados

## Critérios de qualidade para publicação
| Métrica | Threshold |
|---|---|
| Vina score (melhor pose) | ≤ -7 kcal/mol |
| HADDOCK score | ≤ -100 (unidades HADDOCK) |
| RMSD de cluster | ≤ 2.0 Å (poses convergentes) |
| Número de H-bonds (PLIP) | ≥ 3 na interface |
| Resíduo P1 do inibidor | Arg/Lys em bolso S1 (Asp189) |

## Erros comuns e correções
| Erro | Causa | Fix |
|---|---|---|
| Vina score 0 / sem poses | Caixa não cobre sítio ativo | Redefina centro usando coords do Ser195 |
| HADDOCK não converge | AIRs mal definidos | Verifique distâncias < 5 Å entre átomos de interface |
| PLIP: "No interactions" | Complexo não está em contato | Use melhor pose Vina/HADDOCK como input |
| `pdbqt` corrompido | Ligante com valência errada | Use `obabel` para corrigir: `obabel -ipdb -opdbqt` |

## Referência
- Trott & Olson 2010 *J Comp Chem* 31:455 (AutoDock Vina)
- Dominguez et al. 2003 *JACS* 125:1731 (HADDOCK)
- Salentin et al. 2015 *Nucleic Acids Res* 43:W443 (PLIP)
