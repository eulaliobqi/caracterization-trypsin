---
name: trinity-parser
description: Lê, valida e extrai estatísticas de assembly FASTA do Trinity. Trigger: "verificar assembly", "estatísticas FASTA", "quantas sequências".
---

# Trinity Assembly Parser

## Quando usar
- Ao receber um novo assembly Trinity para validar
- Antes da Fase 1 (QC)
- Para gerar estatísticas básicas do input

## Inputs esperados
- `data/raw/trinity_assembly.fasta` (formato Trinity TRINITY_DN*_c*_g*_i*)

## Procedimento

### 1. Validar formato
```bash
# Verificar se é FASTA válido
grep -c "^>" data/raw/trinity_assembly.fasta

# Verificar nomes Trinity
head -20 data/raw/trinity_assembly.fasta | grep "^>"

# Checar sequências com caracteres inválidos
grep -v "^>" data/raw/trinity_assembly.fasta | grep -iP "[^ACGTNRYSWKMBDHV]" | head -5
```

### 2. Estatísticas básicas (seqkit)
```bash
seqkit stats -a data/raw/trinity_assembly.fasta > results/01_qc/assembly_stats.txt
cat results/01_qc/assembly_stats.txt
```

### 3. Distribuição de comprimentos
```python
from Bio import SeqIO
import statistics

seqs = list(SeqIO.parse("data/raw/trinity_assembly.fasta", "fasta"))
lengths = [len(s) for s in seqs]
print(f"N sequências: {len(seqs)}")
print(f"Comprimento médio: {statistics.mean(lengths):.0f} nt")
print(f"N50: {compute_N50(lengths)}")
print(f"Max: {max(lengths)} nt")
print(f"Min: {min(lengths)} nt")
```

## Outputs garantidos
- `results/01_qc/assembly_stats.txt`
- Relatório de N total de transcritos, N50, comprimento médio

## Erros comuns
| Erro | Causa | Fix |
|---|---|---|
| `seqkit: command not found` | Env não ativada | `conda activate assembly_qc` |
| Sequências com `\r` | Arquivo Windows | `sed -i 's/\r//' assembly.fasta` |
| Headers duplicados | Problema no assembly | `seqkit rmdup -n` |
