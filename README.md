# Caracterização Estrutural de Tripsinas de *Anticarsia gemmatalis*

<p align="center">
  <img src="https://img.shields.io/badge/Nextflow-%E2%89%A524.04-brightgreen?logo=nextflow" />
  <img src="https://img.shields.io/badge/Conda%2FMamba-supported-blue?logo=anaconda" />
  <img src="https://img.shields.io/badge/AlphaFold3-v3.0-orange" />
  <img src="https://img.shields.io/badge/GROMACS-2025.1-red" />
  <img src="https://img.shields.io/badge/License-MIT-yellow" />
</p>

> **Pipeline Nextflow DSL2** para identificação, caracterização primária, predição estrutural (AlphaFold3) e simulação de dinâmica molecular de tripsinas do intestino médio larval de *Anticarsia gemmatalis* (lagarta-da-soja, Lepidoptera: Erebidae).

---

## Índice

1. [O que este pipeline faz](#1-o-que-este-pipeline-faz)
2. [Contexto biológico](#2-contexto-biológico)
3. [Requisitos de sistema](#3-requisitos-de-sistema)
4. [Instalação passo a passo](#4-instalação-passo-a-passo)
   - 4.1 [Instalar Nextflow](#41-instalar-nextflow)
   - 4.2 [Instalar Mamba](#42-instalar-mamba)
   - 4.3 [Clonar o repositório](#43-clonar-o-repositório)
   - 4.4 [Criar ambientes Mamba](#44-criar-ambientes-mamba)
   - 4.5 [Baixar bancos de dados de referência](#45-baixar-bancos-de-dados-de-referência)
5. [Preparar seu dado de entrada](#5-preparar-seu-dado-de-entrada)
6. [Executar o pipeline](#6-executar-o-pipeline)
   - 6.1 [Pipeline completo](#61-pipeline-completo)
   - 6.2 [Fase por fase](#62-fase-por-fase)
   - 6.3 [Com GPU (AlphaFold3 e GROMACS)](#63-com-gpu-alphafold3-e-gromacs)
   - 6.4 [Retomar execução interrompida](#64-retomar-execução-interrompida)
   - 6.5 [Teste rápido sem ferramentas instaladas](#65-teste-rápido-sem-ferramentas-instaladas)
7. [Resultados esperados](#7-resultados-esperados)
8. [Estrutura de diretórios](#8-estrutura-de-diretórios)
9. [Referência de parâmetros](#9-referência-de-parâmetros)
10. [Solução de problemas](#10-solução-de-problemas)
11. [Citação](#11-citação)

---

## 1. O que este pipeline faz

```
Assembly Trinity (42.372 transcritos de midgut larval)
         │
         ▼ Fase 1
   ┌─────────────┐
   │  QC Assembly │  BUSCO v5 (completude gênica) + CD-HIT-EST (remove redundância)
   └─────────────┘
         │
         ▼ Fase 2
   ┌──────────────────┐
   │  Predição de ORFs │  TransDecoder + hints BLAST + hints Pfam
   └──────────────────┘
         │
         ▼ Fase 3
   ┌──────────────────────────────────┐
   │  Identificação de Tripsinas      │  DIAMOND BLASTp ∩ HMMER PF00089
   │  (dupla validação)               │  Confident = ambos | Suggestive = um só
   └──────────────────────────────────┘
         │
         ▼ Fase 4
   ┌───────────────────────────────────────────────────────────┐
   │  Filtro de Completude                                     │
   │  Met inicial + ≥220 aa + tríade His57-Asp102-Ser195 + PF00089 ≥80%
   └───────────────────────────────────────────────────────────┘
         │
         ▼ Fase 5
   ┌────────────────────────────────────────────┐
   │  Caracterização Primária                   │
   │  ProtParam (MW, pI) + SignalP6 + InterProScan
   └────────────────────────────────────────────┘
         │
         ▼ Fase 6
   ┌──────────────────────────────┐
   │  Filogenia                   │  MAFFT-linsi + trimAl + IQ-TREE2
   └──────────────────────────────┘
         │
         ▼ Fase 7  [GPU obrigatório]
   ┌────────────────────────────────────────────┐
   │  Predição Estrutural 3D                    │
   │  AlphaFold3 via nf-core/proteinfold        │
   └────────────────────────────────────────────┘
         │
         ▼ Fase 8
   ┌──────────────────────────────────────────────────┐
   │  Validação Estrutural                            │
   │  MolProbity + ProSA + Foldseek + ConSurf         │
   └──────────────────────────────────────────────────┘
         │
         ▼ Fase 9
   ┌──────────────────────────────────────────────────────────┐
   │  Docking Molecular                                       │
   │  AutoDock Vina (peptídeos) + HADDOCK (SKTI) + PLIP       │
   └──────────────────────────────────────────────────────────┘
         │
         ▼ Fase 10  [GPU recomendado]
   ┌──────────────────────────────────────────────────────────┐
   │  Simulações de Dinâmica Molecular (100 ns × 3 réplicas)  │
   │  GROMACS 2025.1 — complexos tripsina-inibidor            │
   └──────────────────────────────────────────────────────────┘
         │
         ▼ Fase 11
   ┌──────────────────────────┐
   │  Relatório e Manuscrito  │  Figuras + Tabelas + Texto
   └──────────────────────────┘
```

---

## 2. Contexto biológico

| Item | Detalhe |
|---|---|
| Organismo | *Anticarsia gemmatalis* Hübner 1818 — lagarta-da-soja |
| Tecido | Intestino médio larval (midgut, pH ~10–11) |
| Alvo | Tripsinas — serina-proteases família S1A (Pfam PF00089) |
| Função biológica | Digestão de proteínas de reserva da soja (*Glycine max*) |
| Tríade catalítica | His⁵⁷ – Asp¹⁰² – Ser¹⁹⁵ (numeração quimotripsina) |
| Tamanho esperado | 20–35 kDa (220–320 aa como pré-zimogênio) |
| Inibidores testados | SKTI (Kunitz) + peptídeos BPTI-like (TGPCK, AVIMK) |
| Isoformas esperadas | 5–30 por espécie (Lepidoptera) |

**Relevância:** A resistência de *A. gemmatalis* a inibidores de tripsina da soja (SKTI) é um mecanismo central de adaptação ao hospedeiro. Caracterizar estruturalmente as isoformas e seus complexos com inibidores é essencial para o desenvolvimento de estratégias de controle.

---

## 3. Requisitos de sistema

### Hardware mínimo (fases 1–6, CPU-only)
| Componente | Mínimo | Recomendado |
|---|---|---|
| CPU | 8 cores | 32–64 cores |
| RAM | 32 GB | 128–256 GB |
| Disco | 500 GB SSD | 2 TB SSD |
| Sistema operacional | Debian 11+ / Ubuntu 20.04+ | Debian 12 |

### Hardware para fases GPU (Fases 7 e 10)
| Componente | Mínimo | Recomendado |
|---|---|---|
| GPU | NVIDIA RTX 3080 (10 GB VRAM) | NVIDIA A100 (40–80 GB) |
| CUDA | 12.0 | 12.4 |
| Driver NVIDIA | 525+ | 550+ |

### Software pré-instalado no sistema
```bash
# Verificar versões
java --version          # Java 17+ (para Nextflow)
python3 --version       # Python 3.10+ (para scripts auxiliares)
git --version           # Git 2.30+
nvidia-smi              # Para GPU (opcional para fases 1-6)
```

---

## 4. Instalação passo a passo

### 4.1 Instalar Nextflow

```bash
# Método 1: Script oficial (recomendado)
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/

# Verificar instalação
nextflow -version
# Saída esperada: nextflow version 24.10.x build XXXX
```

```bash
# Método 2: Via SDKMAN (alternativa)
sdk install java 17.0.9-ms
sdk install nextflow
```

> ⚠️ **Nextflow requer Java 17 ou superior.** Se `java --version` retornar versão antiga:
> ```bash
> sudo apt update && sudo apt install -y openjdk-17-jre
> ```

---

### 4.2 Instalar Mamba

Mamba é um gerenciador de ambientes Conda mais rápido. É obrigatório para instalar todas as ferramentas bioinformáticas.

```bash
# Instalar Miniforge (inclui Mamba)
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -b -p $HOME/miniforge3
source $HOME/miniforge3/bin/activate
conda init bash && source ~/.bashrc

# Verificar
mamba --version
# Saída esperada: mamba X.Y.Z
```

---

### 4.3 Clonar o repositório

```bash
# Clonar
git clone https://github.com/eulaliobqi/caracterization-trypsin.git
cd caracterization-trypsin

# Ver estrutura
ls -la
```

Você deve ver:
```
├── CLAUDE.md
├── LEARNINGS.md
├── README.md
├── envs/           ← ambientes Mamba por fase
├── nextflow/       ← pipeline Nextflow DSL2
│   ├── main.nf
│   ├── modules/
│   └── subworkflows/
├── data/
│   └── raw/        ← coloque seu assembly aqui
└── results/        ← saídas do pipeline
```

---

### 4.4 Criar ambientes Mamba

Cada fase do pipeline usa um ambiente isolado. Crie todos de uma vez:

```bash
cd caracterization-trypsin

# Criar todos os 8 ambientes (pode demorar 20-40 minutos no total)
for env in envs/*.yml; do
    echo "=== Criando ambiente: $env ==="
    mamba env create -f "$env" --yes
done

# Verificar ambientes criados
mamba env list
```

Saída esperada:
```
# conda environments:
assembly_qc          /home/user/miniforge3/envs/assembly_qc
orf_prediction       /home/user/miniforge3/envs/orf_prediction
annotation           /home/user/miniforge3/envs/annotation
phylogeny            /home/user/miniforge3/envs/phylogeny
structure            /home/user/miniforge3/envs/structure
docking              /home/user/miniforge3/envs/docking
md                   /home/user/miniforge3/envs/md
analysis             /home/user/miniforge3/envs/analysis
```

> 💡 **Dica:** Se a criação de um ambiente falhar, tente individualmente:
> ```bash
> mamba env create -f envs/assembly_qc.yml --yes --no-deps --override-channels
> ```

---

### 4.5 Baixar bancos de dados de referência

Os bancos de dados **não estão incluídos** no repositório (são muito grandes). Execute os comandos abaixo **uma única vez** no servidor.

```bash
cd data/references/

# ── 1. UniProt Swiss-Prot (para DIAMOND) ─────────────────────────────────────
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz

# Criar banco DIAMOND (~5 min, requer 10 GB)
mamba run -n annotation diamond makedb \
    --in uniprot_sprot.fasta \
    -d uniprot_sprot \
    --threads 16
# Resultado: uniprot_sprot.dmnd (~3 GB)

# ── 2. Pfam-A (para HMMER — Fases 2 e 3) ────────────────────────────────────
wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz

# Preparar banco HMMER (~10 min)
mamba run -n annotation hmmpress Pfam-A.hmm
# Cria: Pfam-A.hmm.h3f, .h3i, .h3m, .h3p

# ── 3. BUSCO — lineagem insecta_odb10 (para Fase 1) ─────────────────────────
mamba run -n assembly_qc busco --download insecta_odb10 \
    --download_path data/references/
# Resultado: data/references/busco_downloads/lineages/insecta_odb10/

# ── 4. Referências filogenéticas (Lepidoptera trypsins) ──────────────────────
# Baixe manualmente do UniProt: busca por "trypsin Lepidoptera reviewed:yes"
# Salve como: data/references/trypsin_refs_lepidoptera.fasta
# Mínimo recomendado: 30 sequências de Lepidoptera (Spodoptera, Helicoverpa, Manduca...)

echo "=== Bancos de dados prontos ==="
ls -lh data/references/
```

> 📦 **Espaço necessário em disco:**
> | Banco | Tamanho |
> |---|---|
> | UniProt Swiss-Prot (FASTA) | ~280 MB |
> | UniProt Swiss-Prot (DIAMOND .dmnd) | ~3 GB |
> | Pfam-A.hmm | ~8 GB |
> | insecta_odb10 | ~100 MB |
> | **Total** | **~12 GB** |

---

## 5. Preparar seu dado de entrada

O pipeline espera um arquivo FASTA de assembly Trinity:

```bash
# Copie seu assembly para o local correto
cp /caminho/para/seu/assembly.fasta data/raw/trinity_assembly.fasta

# Verificar o arquivo
grep -c "^>" data/raw/trinity_assembly.fasta
# Saída esperada: número de transcritos (ex.: 42372)

# Verificar se é formato Trinity (headers esperados)
head -3 data/raw/trinity_assembly.fasta
# Saída esperada:
# >TRINITY_DN10000_c0_g1_i1 len=500 path=[...]
# ATCGATCGATCG...
```

> ⚠️ **Formato obrigatório:** Headers no padrão Trinity (`TRINITY_DN*_c*_g*_i*`). Se seu assembly veio de outro montador, os módulos continuam funcionando, mas o pipeline foi otimizado para Trinity.

---

## 6. Executar o pipeline

### 6.1 Pipeline completo

```bash
cd caracterization-trypsin

# Executar todas as 11 fases (exceto GPU) no servidor Debian
nextflow run nextflow/main.nf \
    -profile debian \
    --input_fasta data/raw/trinity_assembly.fasta \
    --outdir results \
    -with-report results/pipeline_info/report.html \
    -with-timeline results/pipeline_info/timeline.html \
    -with-dag results/pipeline_info/dag.svg

# Monitorar progresso em tempo real (em outro terminal)
tail -f .nextflow.log
```

---

### 6.2 Fase por fase

Execute apenas até a fase desejada usando `--phase N`. Isso é útil para inspecionar resultados intermediários antes de prosseguir.

```bash
# ── Fase 1: QC do Assembly ──────────────────────────────────────────────────
nextflow run nextflow/main.nf \
    -profile debian \
    --phase 1 \
    --input_fasta data/raw/trinity_assembly.fasta \
    --outdir results
# Saída: results/01_qc/busco/ + results/01_qc/cd_hit/
# Tempo estimado: 1-2 horas

# ── Fase 2: Predição de ORFs ─────────────────────────────────────────────────
nextflow run nextflow/main.nf \
    -profile debian \
    --phase 2 \
    --input_fasta data/raw/trinity_assembly.fasta \
    --outdir results \
    -resume    # ← sempre use -resume para não repetir fases anteriores
# Saída: results/02_orfs/*.transdecoder.pep
# Tempo estimado: 3-6 horas

# ── Fase 3: Identificar Tripsinas (DIAMOND ∩ HMMER) ─────────────────────────
nextflow run nextflow/main.nf \
    -profile debian \
    --phase 3 \
    --input_fasta data/raw/trinity_assembly.fasta \
    --outdir results \
    -resume
# Saída: results/03_trypsin_ids/trypsins_confident.fasta
# Tempo estimado: 30 minutos

# ── Fase 4: Filtro de Completude (tríade catalítica) ────────────────────────
nextflow run nextflow/main.nf \
    -profile debian \
    --phase 4 \
    --input_fasta data/raw/trinity_assembly.fasta \
    --outdir results \
    -resume
# Saída: results/04_complete/complete_trypsins.fasta
# Tempo estimado: 10 minutos

# ── Fase 5: Caracterização primária ─────────────────────────────────────────
nextflow run nextflow/main.nf \
    -profile debian \
    --phase 5 \
    --input_fasta data/raw/trinity_assembly.fasta \
    --outdir results \
    -resume
# Saída: results/05_primary/ (tabela ProtParam, SignalP, InterPro)
# Tempo estimado: 1 hora

# ── Fase 6: Filogenia ────────────────────────────────────────────────────────
nextflow run nextflow/main.nf \
    -profile debian \
    --phase 6 \
    --input_fasta data/raw/trinity_assembly.fasta \
    --outdir results \
    -resume
# Saída: results/06_phylogeny/iqtree2.treefile
# Tempo estimado: 2-4 horas
```

---

### 6.3 Com GPU (AlphaFold3 e GROMACS)

As fases 7 e 10 requerem GPU. Combine o perfil `gpu` com `debian`:

```bash
# ── Fase 7: AlphaFold3 (predição estrutural) — REQUER GPU ───────────────────
nextflow run nextflow/main.nf \
    -profile debian,gpu \
    --phase 7 \
    --input_fasta data/raw/trinity_assembly.fasta \
    --outdir results \
    -resume
# Saída: results/07_structures/alphafold3/*.pdb
# Tempo estimado: 1-4 horas por tripsina (dependendo da GPU)

# ── Fase 10: GROMACS MD 100 ns — REQUER GPU ──────────────────────────────────
nextflow run nextflow/main.nf \
    -profile debian,gpu \
    --phase 10 \
    --input_fasta data/raw/trinity_assembly.fasta \
    --outdir results \
    -resume
# Saída: results/10_md/<tripsina>_<inibidor>/md_100ns_analysis/
# Tempo estimado: 24-72 horas por complexo (3 réplicas)
```

> 🖥️ **Verificar disponibilidade de GPU:**
> ```bash
> nvidia-smi
> # Ou via Nextflow:
> nextflow run nextflow/main.nf -profile debian,gpu --phase 7 -stub-run
> ```

---

### 6.4 Retomar execução interrompida

O Nextflow armazena o estado de cada processo. Para retomar de onde parou:

```bash
# Simples: adicione -resume
nextflow run nextflow/main.nf \
    -profile debian \
    --input_fasta data/raw/trinity_assembly.fasta \
    --outdir results \
    -resume

# Ver histórico de execuções
nextflow log

# Ver processos de uma execução específica
nextflow log <run-name> -f hash,name,status,duration
```

---

### 6.5 Teste rápido sem ferramentas instaladas

Use o modo `stub-run` para verificar se o pipeline parseia corretamente sem executar nenhuma ferramenta:

```bash
# Testar sintaxe e fluxo sem executar nada de verdade
nextflow run nextflow/main.nf \
    -profile test \
    -stub-run

# Saída esperada (sucesso):
# executor > local (N)
# [XX/XXXXXX] process > QC_ASSEMBLY:BUSCO  [100%] 1 of 1 ✔
# [XX/XXXXXX] process > QC_ASSEMBLY:CD_HIT [100%] 1 of 1 ✔
# ...
# Pipeline finished: SUCCESS
```

---

## 7. Resultados esperados

Após execução completa, você terá:

| Fase | Diretório | Arquivos chave | O que esperar |
|---|---|---|---|
| 1 | `results/01_qc/` | `short_summary*.txt` | BUSCO C: > 70% |
| 2 | `results/02_orfs/` | `*.transdecoder.pep` | 10.000–50.000 proteínas |
| 3 | `results/03_trypsin_ids/` | `trypsins_confident.fasta` | 10–100 candidatos |
| 4 | `results/04_complete/` | `complete_trypsins.fasta` | **5–30 tripsinas completas** |
| 5 | `results/05_primary/` | `protparam_summary.tsv` | MW 20–35 kDa, pI 5–9 |
| 6 | `results/06_phylogeny/` | `*.treefile` | Agrupamento por isoforma |
| 7 | `results/07_structures/` | `*_model_0.pdb` | pLDDT > 80 (sítio ativo) |
| 8 | `results/08_validation/` | `foldseek_*.tsv` | TM-score > 0.7 vs. PDB |
| 9 | `results/09_docking/` | `docking_summary.tsv` | Vina ≤ −7 kcal/mol |
| 10 | `results/10_md/` | `rmsd_analysis.tsv` | RMSD < 3 Å em 100 ns |
| 11 | `manuscript/` | `manuscript.md` + figuras | Pronto para submissão |

### Critérios de publicação em Q1 (IF > 5)

```
✓ ≥ 8 isoformas completas identificadas
✓ pLDDT médio dos modelos AF3 > 80
✓ Ramachandran favored > 95%
✓ MD simulations: ≥ 100 ns × ≥ 3 réplicas por complexo
✓ Inibidores testados: ≥ 3 (SKTI + 2 peptídeos BPTI-RCL)
✓ Filogenia: ≥ 30 sequências de Lepidoptera
✓ Foldseek: TM-score > 0.7 vs. referências PDB
✓ Reprodutibilidade: pipeline público no GitHub + DOI Zenodo
```

---

## 8. Estrutura de diretórios

```
caracterization-trypsin/
│
├── CLAUDE.md                     # Contexto persistente (para Claude Code)
├── LEARNINGS.md                  # Decisões, erros e insights documentados
├── ISSUES.md                     # Bugs conhecidos e workarounds
├── README.md                     # Este arquivo
│
├── .claude/
│   ├── agents/                   # 8 subagentes especializados
│   │   ├── bioinformatics-architect.md
│   │   ├── nextflow-developer.md
│   │   ├── mamba-env-manager.md
│   │   ├── code-reviewer.md
│   │   ├── debugger.md
│   │   ├── literature-validator.md
│   │   ├── publication-formatter.md
│   │   └── token-optimizer.md
│   │
│   ├── skills/                   # 10 habilidades reutilizáveis
│   │   ├── trinity-parser/
│   │   ├── transdecoder-runner/
│   │   ├── trypsin-identifier/
│   │   ├── completeness-filter/
│   │   ├── alphafold3-runner/
│   │   ├── foldseek-search/
│   │   ├── docking-pipeline/
│   │   ├── gromacs-md/
│   │   ├── figure-generator/
│   │   └── report-writer/
│   │
│   ├── hooks/                    # Validações automáticas
│   │   ├── pre-commit.sh         # lint + segredos + YAML
│   │   ├── pre-push.sh           # stub-run obrigatório
│   │   └── post-execute.sh       # log de fase concluída
│   │
│   └── commands/                 # Slash commands customizados
│       ├── run-phase.md          # /run-phase N
│       ├── validate.md           # /validate
│       ├── sync-server.md        # /sync-server
│       └── compact-context.md    # /compact-context
│
├── nextflow/
│   ├── main.nf                   # Pipeline principal (controle de fase)
│   ├── nextflow.config           # Configuração geral + parâmetros
│   ├── conf/
│   │   ├── base.config           # Defaults de recursos
│   │   ├── debian.config         # Perfil servidor (64 cores, 256 GB)
│   │   ├── gpu.config            # Perfil GPU (A100)
│   │   └── test.config           # Mini-dataset para CI
│   │
│   ├── modules/local/            # 19 módulos Nextflow DSL2
│   │   ├── busco.nf              # Fase 1
│   │   ├── cd_hit.nf             # Fase 1
│   │   ├── transdecoder.nf       # Fase 2
│   │   ├── diamond.nf            # Fases 2, 3 (incl. INTERSECT_TRYPSINS)
│   │   ├── hmmer.nf              # Fases 2, 3
│   │   ├── interproscan.nf       # Fase 5
│   │   ├── signalp.nf            # Fase 5
│   │   ├── mafft.nf              # Fase 6
│   │   ├── trimal.nf             # Fase 6
│   │   ├── iqtree.nf             # Fase 6
│   │   ├── alphafold3.nf         # Fase 7
│   │   ├── foldseek.nf           # Fase 8
│   │   ├── tmalign.nf            # Fase 8
│   │   ├── consurf.nf            # Fase 8
│   │   ├── molprobity.nf         # Fase 8
│   │   ├── prosa.nf              # Fase 8
│   │   ├── autodock_vina.nf      # Fase 9
│   │   ├── haddock.nf            # Fase 9
│   │   └── gromacs.nf            # Fase 10
│   │
│   ├── subworkflows/             # 11 subworkflows
│   └── bin/                      # Scripts Python auxiliares
│       ├── filter_complete_trypsins.py   # Filtro de completude + interseção
│       ├── extract_catalytic_triad.py    # Verifica His/Asp/Ser por alinhamento
│       ├── compute_protparam.py          # Massa molecular, pI, GRAVY
│       ├── parse_alphafold3.py           # Extrai pLDDT, PTM dos JSONs AF3
│       └── make_figures.py               # Gera figuras para manuscrito
│
├── envs/                         # 8 ambientes Mamba (.yml)
│   ├── assembly_qc.yml           # BUSCO, CD-HIT, seqkit
│   ├── orf_prediction.yml        # TransDecoder, DIAMOND, HMMER
│   ├── annotation.yml            # DIAMOND, HMMER, Biopython
│   ├── phylogeny.yml             # MAFFT, trimAl, IQ-TREE2
│   ├── structure.yml             # AlphaFold3, Foldseek, PyMOL, ConSurf
│   ├── docking.yml               # AutoDock Vina, HADDOCK, PLIP
│   ├── md.yml                    # GROMACS 2025.1, MDAnalysis, MDTraj
│   └── analysis.yml              # Python, matplotlib, seaborn, pandas
│
├── data/
│   ├── raw/
│   │   └── trinity_assembly.fasta   # ← Seu assembly aqui (não versionado)
│   ├── references/                  # ← Bancos de dados aqui (não versionados)
│   │   ├── uniprot_sprot.dmnd
│   │   ├── Pfam-A.hmm
│   │   ├── insecta_odb10/
│   │   ├── trypsin_refs_lepidoptera.fasta
│   │   └── mdp/                     # Parâmetros GROMACS
│   └── intermediate/
│
├── results/                      # Saídas do pipeline (não versionadas)
│   ├── 01_qc/ ... 11_report/
│   └── pipeline_info/            # Timeline, report, DAG do Nextflow
│
├── manuscript/                   # Manuscrito científico
│   ├── manuscript.md
│   ├── figures/
│   ├── tables/
│   └── supplementary/
│
└── tests/
    └── data/
        └── test_assembly.fasta   # Mini-assembly para -stub-run e CI
```

---

## 9. Referência de parâmetros

Todos os parâmetros podem ser sobrescritos na linha de comando com `--nome valor`.

| Parâmetro | Padrão | Descrição |
|---|---|---|
| `--input_fasta` | `data/raw/trinity_assembly.fasta` | Assembly Trinity de entrada |
| `--outdir` | `results` | Diretório de saída |
| `--phase` | `null` (todas) | Executar até fase N (1–11) |
| `--diamond_evalue` | `1e-10` | E-value mínimo para DIAMOND |
| `--hmmer_evalue` | `1e-10` | E-value mínimo para HMMER |
| `--min_aa_length` | `220` | Comprimento mínimo das tripsinas (aa) |
| `--blast_coverage` | `0.80` | Cobertura mínima do domínio PF00089 |
| `--blast_identity` | `30` | Identidade mínima DIAMOND (%) |
| `--max_cpus` | `32` | Máximo de CPUs por processo |
| `--max_memory` | `128.GB` | Máximo de RAM por processo |
| `--max_time` | `72.h` | Tempo máximo por processo |
| `--use_gpu` | `false` | Ativar suporte a GPU |
| `--uniprot_db` | `data/references/uniprot_sprot.dmnd` | Banco DIAMOND |
| `--pfam_db` | `data/references/Pfam-A.hmm` | Banco Pfam |
| `--busco_lineage` | `data/references/insecta_odb10` | Lineagem BUSCO |
| `--lep_refs` | `data/references/trypsin_refs_lepidoptera.fasta` | Refs filogenia |

**Exemplo — E-value mais restrito e CPU máximo aumentado:**
```bash
nextflow run nextflow/main.nf \
    -profile debian \
    --diamond_evalue 1e-15 \
    --hmmer_evalue 1e-15 \
    --max_cpus 64 \
    --max_memory 256.GB \
    -resume
```

---

## 10. Solução de problemas

### ❌ `Error: Command not found: mamba`
```bash
# Ativar Miniforge
source $HOME/miniforge3/bin/activate
conda init bash && source ~/.bashrc
```

### ❌ `No such file: data/references/uniprot_sprot.dmnd`
```bash
# O banco DIAMOND não foi criado
cd data/references/
diamond makedb --in uniprot_sprot.fasta -d uniprot_sprot --threads 16
```

### ❌ `BUSCO error: offline mode`
```bash
# BUSCO não encontrou a lineagem localmente
# Especifique o caminho completo
nextflow run nextflow/main.nf \
    --busco_lineage /caminho/absoluto/para/insecta_odb10 \
    -profile debian --phase 1
```

### ❌ `0 confident trypsins found`
Possíveis causas:
1. **E-value muito restrito:** tente `--diamond_evalue 1e-5 --hmmer_evalue 1e-5`
2. **Assembly de má qualidade:** verifique BUSCO C < 50%
3. **Pfam desatualizado:** confirme versão ≥ 36.0 com `head Pfam-A.hmm | grep "# HMMER"`
4. **TransDecoder sem ORFs longas:** verifique `results/02_orfs/` — zero ORFs = assembly fragmentado

### ❌ `AlphaFold3: CUDA out of memory`
```bash
# Reduzir número de recycles ou processar uma tripsina por vez
nextflow run nextflow/main.nf \
    -profile debian,gpu \
    --phase 7 \
    -resume \
    --max_memory 32.GB  # reduzir se GPU com < 40 GB VRAM
```

### ❌ `Nextflow hung / process never ends`
```bash
# Ver logs detalhados
tail -200 .nextflow.log | grep -E "ERROR|WARN|FAILED"

# Cancelar e retomar
Ctrl+C
nextflow run nextflow/main.nf [mesmos parâmetros] -resume
```

### ℹ️ Ver relatório detalhado da execução
```bash
# Abrir no navegador
firefox results/pipeline_info/report.html   # relatório HTML
firefox results/pipeline_info/timeline.html  # timeline de processos
firefox results/pipeline_info/dag.svg        # grafo do pipeline
```

---

## 11. Citação

Se você usar este pipeline em sua pesquisa, por favor cite:

```bibtex
@misc{gutemberg2026trypsins,
  author       = {Gutemberg, Eulálio},
  title        = {Structural characterization of digestive trypsins from
                  \textit{Anticarsia gemmatalis} (Lepidoptera: Erebidae)},
  year         = {2026},
  publisher    = {GitHub},
  journal      = {GitHub repository},
  howpublished = {\url{https://github.com/eulaliobqi/caracterization-trypsin}},
  note         = {DOI: 10.5281/zenodo.XXXXXXX (to be assigned)}
}
```

### Ferramentas que devem ser citadas em seu manuscrito

| Ferramenta | Referência |
|---|---|
| Nextflow | Di Tommaso et al. 2017 *Nat Biotechnol* 35:316 |
| Trinity | Grabherr et al. 2011 *Nat Biotechnol* 29:644 |
| BUSCO v5 | Manni et al. 2021 *Mol Biol Evol* 38:4647 |
| CD-HIT | Fu et al. 2012 *Bioinformatics* 28:3150 |
| TransDecoder | Haas et al. 2013 *Nat Protoc* 8:1494 |
| DIAMOND | Buchfink et al. 2021 *Nat Methods* 18:366 |
| HMMER 3.4 | Eddy 2011 *PLOS Comput Biol* 7:e1002195 |
| Pfam | Mistry et al. 2021 *Nucleic Acids Res* 49:D412 |
| MAFFT | Katoh & Standley 2013 *Mol Biol Evol* 30:772 |
| trimAl | Capella-Gutiérrez et al. 2009 *Bioinformatics* 25:1972 |
| IQ-TREE2 | Minh et al. 2020 *Mol Biol Evol* 37:1530 |
| AlphaFold3 | Abramson et al. 2024 *Nature* 630:493 |
| Foldseek | van Kempen et al. 2024 *Nature* 625:832 |
| GROMACS 2025 | Abraham et al. 2015 *SoftwareX* 1-2:19 |
| AutoDock Vina | Trott & Olson 2010 *J Comp Chem* 31:455 |
| HADDOCK | Dominguez et al. 2003 *JACS* 125:1731 |
| PLIP | Salentin et al. 2015 *Nucleic Acids Res* 43:W443 |

---

## Licença

MIT License — livre para uso acadêmico e comercial com atribuição.

---

<p align="center">
  Desenvolvido com ❤️ para caracterização de tripsinas de <em>Anticarsia gemmatalis</em><br>
  Contato: <a href="mailto:eulalio.santos@ufv.br">eulalio.santos@ufv.br</a>
</p>
