---
name: foldseek-search
description: Busca estrutural ultra-rápida com Foldseek contra AFDB + PDB por TM-score. Trigger: "busca estrutural", "Foldseek", "TM-score", "estruturas homólogas", após AlphaFold3 produzir PDBs (Fase 8).
---

# Foldseek Search (Fase 8)

## Quando usar
- Após `STRUCTURE_PRED` produzir PDBs/mmCIF das tripsinas
- Usuário pede "busca estrutural", "encontrar homólogos estruturais", "TM-score"
- Validação de qualidade dos modelos AF3 (Fase 8)

## Inputs esperados
- `model.pdb` (ou `.cif`) — estrutura predita pelo AlphaFold3
- `afdb` — base AFDB (AlphaFold Database, comprimida pelo Foldseek)
- `pdb100` — base PDB (todas as estruturas depositadas)

## Procedimento

### Download das bases (uma vez no servidor)
```bash
# Base PDB (atualizada semanalmente)
foldseek databases PDB pdb_db /tmp --threads 16

# Base AFDB (AlphaFold Database)
foldseek databases Alphafold/UniProt50 afdb_db /tmp --threads 16

# Opção mais compacta: apenas SwissProt no AFDB
foldseek databases Alphafold/Swiss-Prot afdb_sprot_db /tmp --threads 16
```

### Busca estrutural (cada modelo de tripsina)
```bash
# Busca contra PDB + AFDB em série
foldseek easy-search \
    ${meta.id}_model_0.pdb \
    pdb_db \
    foldseek_pdb.tsv \
    /tmp/foldseek_tmp \
    --format-output "query,target,alntmscore,qtmscore,ttmscore,evalue,bits,rmsd,qlen,tlen,qcov,tcov" \
    --exhaustive-search 0 \
    --threads ${task.cpus} \
    -e 0.001

foldseek easy-search \
    ${meta.id}_model_0.pdb \
    afdb_db \
    foldseek_afdb.tsv \
    /tmp/foldseek_tmp \
    --format-output "query,target,alntmscore,qtmscore,ttmscore,evalue,bits,rmsd,qlen,tlen,qcov,tcov" \
    --exhaustive-search 0 \
    --threads ${task.cpus} \
    -e 0.001
```

### Filtrar hits de alta qualidade
```python
import pandas as pd

cols = ["query","target","alntmscore","qtmscore","ttmscore",
        "evalue","bits","rmsd","qlen","tlen","qcov","tcov"]

df = pd.read_csv("foldseek_pdb.tsv", sep="\t", names=cols)

# Threshold: TM-score > 0.5 (estruturalmente similar)
# TM-score > 0.7 = mesma topologia (alta confiança)
high_conf = df[df["alntmscore"] > 0.5].sort_values("alntmscore", ascending=False)
high_conf.to_csv("foldseek_filtered.tsv", sep="\t", index=False)
```

## Outputs garantidos
- `results/08_validation/foldseek/<meta.id>_pdb_hits.tsv`   — hits no PDB
- `results/08_validation/foldseek/<meta.id>_afdb_hits.tsv`  — hits no AFDB
- `results/08_validation/foldseek/<meta.id>_filtered.tsv`   — hits TM-score > 0.5

## Interpretação dos resultados
| TM-score | Interpretação |
|---|---|
| > 0.9 | Estruturas quasi-idênticas (mesmo fold, alta identidade) |
| 0.7–0.9 | Mesma topologia (estruturalmente equivalentes) |
| 0.5–0.7 | Mesmo fold, divergência moderada (esperado para tripsinas de insetos) |
| < 0.5 | Diferentes folds (descartar ou investigar) |

**Para tripsinas de *A. gemmatalis*:** esperar TM-score > 0.7 contra tripsina bovina (1TGN) e outros insetos Lepidoptera.

## Qualidade mínima para publicação (Q1 journals)
- ≥ 1 hit PDB com TM-score > 0.7 para cada modelo
- Melhor hit esperado: tripsinas de insetos (Manduca, Heliothis, Spodoptera)
- RMSD < 2.5 Å para região do sítio ativo

## Erros comuns e correções
| Erro | Causa | Fix |
|---|---|---|
| `Database not found` | Path da DB errado | Verifique `${params.foldseek_db_pdb}` no nextflow.config |
| `0 hits` | Arquivo PDB corrompido | Verifique se AlphaFold3 gerou PDB válido com `check_structure.py` |
| Timeout | Busca exhaustiva em DB grande | Use `--exhaustive-search 0` (aproximado mas 1000x mais rápido) |
| TM-score todos < 0.5 | Modelo ruim (pLDDT < 60) | Filtre modelos por pLDDT > 70 antes de buscar |

## Referência
- van Kempen et al. 2024 *Nature* 625:832-839 (Foldseek)
- Versão recomendada: Foldseek ≥ 9.427df8a
