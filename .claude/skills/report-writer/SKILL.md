---
name: report-writer
description: Gera seções Methods e Results do manuscrito científico com linguagem passiva, estatísticas e citações. Trigger: "escrever manuscrito", "seção de métodos", "seção de resultados", "redigir artigo", na Fase 11.
---

# Report Writer (Fase 11)

## Quando usar
- Após conclusão de todas as fases analíticas (Fases 1-10)
- Usuário pede "escrever manuscrito", "redigir métodos", "redigir resultados"
- Chamado pelo `publication-formatter` para gerar texto base

## Revistas alvo
| Revista | IF | Limite palavras | Foco |
|---|---|---|---|
| Insect Biochem Mol Biol | ~4 | 8.000 words | Bioquímica de insetos |
| Int J Biol Macromol | ~8 | 10.000 words | Caracterização de proteínas |
| Comput Struct Biotechnol J | ~6 | Sem limite | Bioinformática estrutural |

## Template de seção Methods

### 2. Materials and Methods

#### 2.1 Transcriptome Assembly and Quality Control
```
The transcriptome of *Anticarsia gemmatalis* larval midgut was assembled using
Trinity v2.15.1 [CITE], yielding 42,372 transcripts. Assembly completeness was
evaluated with BUSCO v5.7.1 [CITE] against the *Insecta* lineage (insecta_odb10,
n=1,013 genes). Transcript redundancy was reduced using CD-HIT-EST v4.8.1 [CITE]
with a 95% identity threshold (-c 0.95 -aS 0.9).
```

#### 2.2 ORF Prediction
```
Open reading frames (ORFs) were predicted using TransDecoder v5.7.0 [CITE].
Candidate ORFs (≥100 aa) were retained for initial evaluation. Prediction
accuracy was improved using two independent evidence sources: (i) homology
searches against UniProt/Swiss-Prot (release YYYY-MM) with DIAMOND v2.1 [CITE]
(E-value ≤ 1×10⁻⁵), and (ii) domain-based evidence from HMMER v3.4 [CITE]
against Pfam-A (v36.0) [CITE].
```

#### 2.3 Trypsin Identification
```
Putative trypsins were identified by dual-evidence consensus. Protein sequences
were searched against UniProt/Swiss-Prot with DIAMOND BLASTp (E-value ≤ 1×10⁻¹⁰,
--sensitive mode) and annotated as trypsin-containing when the best hit included
"trypsin" in the description. Concurrently, HMMER 3.4 was used to detect the
serine protease trypsin domain (Pfam: PF00089/Tryp_SPc, E-value ≤ 1×10⁻¹⁰)
using hmmfetch/hmmsearch. Sequences detected by *both* methods were classified as
"confident trypsins"; those detected by only one method were classified as
"suggestive" and subjected to manual inspection.
```

#### 2.4 Completeness Filtering (Criteria)
```
Confident trypsin sequences were subjected to four completeness criteria:
(i) presence of an N-terminal methionine (start codon verified);
(ii) minimum length of 220 amino acids (consistent with full-length trypsin
zymogen in Lepidoptera [CITE]);
(iii) conservation of the catalytic triad (His⁵⁷–Asp¹⁰²–Ser¹⁹⁵, chymotrypsin
numbering) verified by pairwise alignment with bovine trypsin (UniProt: P00760)
using Biopython v1.83 [CITE]; and
(iv) Pfam PF00089 domain coverage ≥ 80% at E-value ≤ 1×10⁻¹⁰.
Sequences meeting all four criteria were designated "complete trypsins."
```

## Template de seção Results

### 3. Results

#### 3.1 Assembly QC
```
De novo transcriptome assembly yielded 42,372 transcripts (median length: XXX nt;
N50: XXX nt). BUSCO analysis revealed XX% completeness (C:XX%[S:XX%,D:XX%],
F:XX%, M:XX%, n=1,013), indicating [high/moderate] representation of conserved
insect genes. Following CD-HIT-EST clustering at 95% identity, the non-redundant
assembly comprised XX,XXX transcripts, representing a XX% reduction in
redundancy.
```

#### 3.2 Trypsin Identification
```
TransDecoder predicted XX,XXX ORFs from the non-redundant assembly, of which
XX,XXX were ≥ 220 amino acids. DIAMOND BLASTp identified XX sequences with
trypsin-matching annotations, while HMMER detected XX sequences carrying the
Tryp_SPc domain (PF00089). The intersection of both approaches yielded XX
confident trypsin candidates; an additional XX sequences detected by a single
method were retained as "suggestive" for manual review (Figure 2A–B).
Completeness filtering reduced the confident set to N complete trypsin isoforms
(N ranges 5–30 in Lepidoptera [CITE]), all possessing a canonical N-terminal
methionine, a complete catalytic triad, and PF00089 coverage > 80%.
```

## Procedimento de geração

### Compilar estatísticas dos resultados
```python
# collect_stats.py — gera tabela unificada de todos resultados
import json, pandas as pd
from pathlib import Path

stats = {}

# Fase 1
busco = open("results/01_qc/busco/short_summary*.txt").read()
stats['busco_complete'] = float(re.search(r'C:(\d+\.\d+)', busco).group(1))

# Fase 3
df = pd.read_csv("results/03_trypsin_ids/identification_report.tsv", sep='\t')
stats['confident'] = df['confident'].sum()
stats['suggestive'] = (~df['confident']).sum()

# Fase 4
df4 = pd.read_csv("results/04_complete/completeness_report.tsv", sep='\t')
stats['complete_trypsins'] = df4['pass_completeness'].sum()

json.dump(stats, open("manuscript/stats.json", "w"), indent=2)
```

## Outputs garantidos
- `manuscript/manuscript.md` — rascunho completo (Methods + Results)
- `manuscript/tables/table1_trypsins.tsv` — tabela principal
- `manuscript/tables/table2_structures.tsv` — tabela de estruturas
- `manuscript/stats.json` — todas estatísticas consolidadas

## Erros comuns e correções
| Erro | Causa | Fix |
|---|---|---|
| Estatísticas inconsistentes | Arquivos de resultados desatualizados | Re-rodar `collect_stats.py` |
| Citações faltando | Referências não mapeadas | Consultar `LEARNINGS.md` para papers chave |
| Métodos sem versões de software | Esqueceu versões | Extrair de `versions.yml` de cada processo |
| Resultados sem IC/p-values | Análises estatísticas não feitas | Adicionar teste de Kruskal-Wallis para RMSD entre réplicas MD |
