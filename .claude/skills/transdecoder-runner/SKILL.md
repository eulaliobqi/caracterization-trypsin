---
name: transdecoder-runner
description: Executa predição de ORFs com TransDecoder em 4 passos (LongOrfs → BLAST hint → Pfam hint → Predict). Trigger: "predizer ORFs", "TransDecoder", "proteoma predito", após QC assembly da Fase 1.
---

# TransDecoder Runner (Fase 2)

## Quando usar
- Após `QC_ASSEMBLY` produzir `clean_fasta`
- Usuário pede "predizer ORFs", "gerar proteoma", "executar TransDecoder"
- Antes de qualquer identificação de tripsinas (Fase 3 depende do .pep)

## Inputs esperados
- `clean_fasta` — FASTA do assembly pós-CD-HIT-EST (Fase 1)
- `uniprot_sprot.dmnd` — DIAMOND DB (para hints BLAST)
- `Pfam-A.hmm` — Base Pfam completa (para hints de domínio)

## Procedimento (4 passos obrigatórios)

### Passo 1 — Extrair ORFs longas
```bash
TransDecoder.LongOrfs \
    -t clean_fasta \
    -m 100            # mínimo 100 aa para ORFs candidatas
```

### Passo 2 — BLAST hint (melhora predição com evidência de homologia)
```bash
diamond blastp \
    --query longest_orfs.pep \
    --db uniprot_sprot.dmnd \
    --outfmt 6 qseqid sseqid evalue \
    --evalue 1e-5 --max-target-seqs 1 --sensitive \
    --threads ${task.cpus} \
    --out blastp.outfmt6
```

### Passo 3 — Pfam hint (detecta domínios conservados como âncora)
```bash
hmmscan \
    --domtblout pfam.domtblout \
    --cpu ${task.cpus} \
    -E 1e-10 \
    Pfam-A.hmm \
    longest_orfs.pep > /dev/null
```

### Passo 4 — Predição final com hints
```bash
TransDecoder.Predict \
    -t clean_fasta \
    --retain_blastp_hits blastp.outfmt6 \
    --retain_pfam_hits pfam.domtblout \
    --single_best_only    # uma ORF por transcrito (menos redundância)
```

### Passo 5 — Converter coordenadas para genoma (opcional mas bom para GFF3)
```bash
TransDecoder.Util/gtf_to_alignment_gff3.pl trinity.gtf > trinity.gff3
TransDecoder.Util/cdna_alignment_orf_to_genome_orf.pl \
    clean_fasta.transdecoder.gff3 \
    trinity.gff3 \
    clean_fasta > final.gff3
```

## Outputs garantidos
- `results/02_orfs/<meta.id>.transdecoder.pep`  — proteínas preditas (input Fase 3)
- `results/02_orfs/<meta.id>.transdecoder.cds`  — CDS correspondentes
- `results/02_orfs/<meta.id>.transdecoder.gff3` — coordenadas genômicas

## Parâmetros críticos para *A. gemmatalis*
| Parâmetro | Valor | Justificativa |
|---|---|---|
| `-m` (min aa) | 100 | Conservador; filtro real de ≥220 aa é na Fase 4 |
| `--single_best_only` | ativo | Reduz redundância de isoformas |
| `--retain_blastp_hits` | obrigatório | Melhora recall de tripsinas divergentes |
| `--retain_pfam_hits` | obrigatório | Ancora ORFs em domínios Pfam conhecidos |

## Erros comuns e correções
| Erro | Causa | Fix |
|---|---|---|
| `No such file: longest_orfs.pep` | LongOrfs não rodou antes de Predict | Garanta ordem dos passos |
| `0 ORFs predicted` | FASTA com headers problemáticos | Use `seqkit replace` para limpar headers |
| Muitas ORFs fragmentadas | `-m` muito baixo | Aumente para 100-150 |
| BLAST hint vazio | DB errada ou FASTA nucleotídeo passado ao DIAMOND | Confirme que `.pep` (proteína) vai ao DIAMOND |
| `hmmscan: command not found` | HMMER não instalado | `mamba install -n orf_prediction bioconda::hmmer` |

## Referência
- Haas et al. 2013 *Nature Protocols* 8:1494-512 (TransDecoder original)
- Versão recomendada: TransDecoder ≥ 5.7.0
