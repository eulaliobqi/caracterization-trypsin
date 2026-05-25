---
name: alphafold3-runner
description: Executa AlphaFold3 via nf-core/proteinfold para predição estrutural das tripsinas completas. Trigger: "predizer estrutura", "AlphaFold3", "modelagem 3D". Requer GPU.
---

# AlphaFold3 Runner

## Quando usar
- Fase 7: após confirmar lista final de tripsinas completas
- Requer GPU (RTX 5070 Ti ou similar CUDA 12.4+)

## Inputs
- `results/04_complete/complete_trypsins.fasta`
- AlphaFold3 model weights (baixar separadamente do Google)
- Bancos de sequências: UniRef90, UniClust30, PDB, BFD

## Preparação dos inputs

### Converter FASTA para JSON (formato AF3)
```python
# nextflow/bin/parse_alphafold3.py
import json
from Bio import SeqIO

for rec in SeqIO.parse("complete_trypsins.fasta", "fasta"):
    af3_input = {
        "name": rec.id,
        "sequences": [{"proteinChain": {"sequence": str(rec.seq), "count": 1}}]
    }
    with open(f"af3_inputs/{rec.id}.json", "w") as f:
        json.dump(af3_input, f, indent=2)
```

## Execução via módulo Nextflow
```bash
nextflow run nextflow/main.nf \
    -profile debian,gpu \
    --phase 7 \
    --outdir results \
    --input_fasta data/raw/trinity_assembly.fasta \
    -resume
```

## Parâmetros recomendados para tripsinas
```
--num_recycles 10          # maior precisão (padrão é 3)
--max_template_date 2024-01-01
--use_gpu_relax true
```

## Avaliação dos modelos
```python
# Critérios de qualidade (ver literature-validator)
# pLDDT > 80 = alta confiança
# pLDDT 60-80 = confiança moderada (loop de ativação normal)
# pLDDT < 60 = regiões desordenadas (propeptídeo N-terminal)

import json
conf = json.load(open("confidence.json"))
plddt_mean = sum(conf['atom_plddts']) / len(conf['atom_plddts'])
print(f"pLDDT médio: {plddt_mean:.1f}")
```

## Outputs
- `results/07_structures/alphafold3/{id}/{id}_model_0.pdb`
- `results/07_structures/alphafold3/{id}/{id}_confidence.json`
- pLDDT global e por resíduo

## Erros comuns
| Erro | Causa | Fix |
|---|---|---|
| CUDA OOM | GPU VRAM insuficiente | Reduzir max_recycling_iters=5 |
| `weights not found` | Caminho errado | Checar $ALPHAFOLD3_MODELS |
| Timeout após 24h | Muitas sequências | Processar em lotes de 10 |
