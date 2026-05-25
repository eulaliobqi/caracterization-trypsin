---
name: trypsin-identifier
description: Identifica tripsinas em proteoma predito via dupla validação DIAMOND + HMMER (Pfam Tryp_SPc PF00089). Trigger: "identificar tripsinas", "DIAMOND HMMER", após TransDecoder produzir .pep.
---

# Trypsin Identifier (DIAMOND + HMMER)

## Quando usar
- Após TransDecoder produzir `*.pep` (Fase 2 concluída)
- Para qualquer screen de serina-proteases classe S1
- Quando usuário pede "identificar tripsinas" ou "buscar serino-proteases"

## Inputs
- `results/02_orfs/transdecoder.pep` (TransDecoder output)
- `data/references/uniprot_sprot.dmnd` (DIAMOND DB — criar com `diamond makedb`)
- `data/references/Pfam-A.hmm` (versão ≥ 36.0, com PF00089)

## Procedimento

### 1. DIAMOND BLASTp (mais rápido que BLAST, padrão 2025)
```bash
diamond blastp \
  --query results/02_orfs/transdecoder.pep \
  --db data/references/uniprot_sprot.dmnd \
  --outfmt 6 qseqid sseqid stitle pident length qlen slen evalue bitscore \
  --evalue 1e-10 \
  --max-target-seqs 5 \
  --sensitive \
  --threads ${SLURM_CPUS_PER_TASK:-8} \
  --out results/03_trypsin_ids/diamond_hits.tsv
```

### 2. HMMER scan com PF00089 (Tryp_SPc)
```bash
# Extrair apenas o perfil de tripsina
hmmfetch data/references/Pfam-A.hmm Tryp_SPc > data/references/trypsin.hmm
hmmpress data/references/trypsin.hmm

# Buscar
hmmsearch \
  --domtblout results/03_trypsin_ids/trypsin_hmm.out \
  --cpu ${SLURM_CPUS_PER_TASK:-8} \
  -E 1e-10 \
  data/references/trypsin.hmm \
  results/02_orfs/transdecoder.pep
```

### 3. Interseção robusta (Python)
```python
import pandas as pd

# Parse DIAMOND — filtrar por keyword "trypsin" no stitle
diamond = pd.read_csv("results/03_trypsin_ids/diamond_hits.tsv",
                      sep='\t', header=None,
                      names=['qseqid','sseqid','stitle','pident','length',
                             'qlen','slen','evalue','bitscore'])
diamond_trypsin = set(
    diamond[diamond['stitle'].str.lower().str.contains('trypsin|serine protease')]['qseqid']
)

# Parse HMMER
hmm_ids = set()
with open("results/03_trypsin_ids/trypsin_hmm.out") as f:
    for line in f:
        if not line.startswith('#'):
            hmm_ids.add(line.split()[0])

# Interseção = alta confiança
confident_ids = diamond_trypsin & hmm_ids
suggestive_ids = (diamond_trypsin | hmm_ids) - confident_ids

print(f"Confident trypsins: {len(confident_ids)}")
print(f"Suggestive: {len(suggestive_ids)}")
```

### 4. Extrair FASTAs
```bash
seqkit grep -n -f confident_ids.txt results/02_orfs/transdecoder.pep \
  > results/03_trypsin_ids/trypsins_confident.fasta
seqkit grep -n -f suggestive_ids.txt results/02_orfs/transdecoder.pep \
  > results/03_trypsin_ids/trypsins_suggestive.fasta
```

## Outputs
- `results/03_trypsin_ids/trypsins_confident.fasta` → interseção DIAMOND+HMMER
- `results/03_trypsin_ids/trypsins_suggestive.fasta` → para revisão manual
- `results/03_trypsin_ids/identification_report.tsv` → tabela completa

## Erros comuns
| Erro | Causa | Fix |
|---|---|---|
| 0 hits HMMER | Pfam DB desatualizado ou PF00089 ausente | Baixe Pfam ≥ 36.0 |
| Muitos hits "trypsin-like" | E-value frouxo | Use 1e-15 ao invés de 1e-10 |
| DIAMOND db error | Versão diamond ≠ db | Recrie .dmnd com mesma versão |
| Muitos suggestive | Keywords muito amplas | Restringir a "trypsin" apenas |
