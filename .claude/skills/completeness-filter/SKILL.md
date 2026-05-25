---
name: completeness-filter
description: Filtra tripsinas completas por tríade catalítica (His57-Asp102-Ser195) + Met inicial + ≥220 aa + cobertura PF00089. Trigger: "filtrar tripsinas completas", "tríade catalítica", "sequências completas".
---

# Completeness Filter

## Quando usar
- Após `trypsin-identifier` gerar `trypsins_confident.fasta`
- Para selecionar apenas isoformas funcionais completas

## Critérios (todos obrigatórios — AND lógico)
1. **Met inicial:** ORF começa com Metionina
2. **Comprimento ≥ 220 aa:** cobre pré-proteína com propeptídeo
3. **Tríade catalítica:** His, Asp, Ser nas posições conservadas do alinhamento
4. **Domínio Pfam PF00089:** cobertura ≥ 80%, E-value ≤ 1e-10

## Procedimento

### 1. Filtro por Met inicial e comprimento
```python
from Bio import SeqIO
passing = []
for rec in SeqIO.parse("results/03_trypsin_ids/trypsins_confident.fasta", "fasta"):
    if rec.seq[0] == 'M' and len(rec.seq) >= 220:
        passing.append(rec)
SeqIO.write(passing, "results/04_complete/step1_met_len.fasta", "fasta")
```

### 2. Detecção da tríade catalítica por alinhamento
```bash
# Alinhar contra referência bovine trypsin (1TGN_A)
mafft --add results/04_complete/step1_met_len.fasta \
      --reorder \
      data/references/bovine_trypsin_1TGN.fasta \
      > results/04_complete/aligned_with_ref.fasta

# Script Python: verificar His57, Asp102, Ser195 na coluna do alinhamento
python3 nextflow/bin/extract_catalytic_triad.py \
    --alignment results/04_complete/aligned_with_ref.fasta \
    --reference "1TGN_A" \
    --output results/04_complete/triad_check.tsv
```

### 3. Filtro final
```python
# Manter apenas sequências que passam TODOS os critérios
# Ver nextflow/bin/filter_complete_trypsins.py
python3 nextflow/bin/filter_complete_trypsins.py \
    --input results/04_complete/step1_met_len.fasta \
    --triad results/04_complete/triad_check.tsv \
    --hmm_coverage results/03_trypsin_ids/trypsin_hmm.out \
    --min_coverage 0.80 \
    --output results/04_complete/complete_trypsins.fasta
```

## Outputs
- `results/04_complete/complete_trypsins.fasta` → candidatos finais para publicação
- `results/04_complete/completeness_report.tsv` → critérios por sequência

## Expectativa biológica
- *A. gemmatalis*: 5–30 tripsinas completas esperadas
- Se < 5: revisar parâmetros do filtro (talvez cobertura mínima muito alta)
- Se > 50: possível redundância — revisar CD-HIT-EST da Fase 1
