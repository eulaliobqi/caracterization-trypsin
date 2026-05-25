# PROJETO: Caracterização Estrutural de Tripsinas de *Anticarsia gemmatalis*

**Versão:** 1.0
**Pesquisador principal:** Eulálio Gutemberg
**Stack:** Nextflow + Mamba/Conda + Claude Code (local) → GitHub → Servidor Debian
**Objetivo final:** Artigo científico de alto impacto (IF > 5) caracterizando estruturalmente tripsinas completas montadas via Trinity.

---

## 0. INSTRUÇÕES MESTRAS PARA CLAUDE CODE

> Leia este arquivo INTEIRO antes de qualquer ação. Ele é a fonte da verdade do projeto.

### 0.1 Princípios não-negociáveis

1. **Reprodutibilidade total** — Todo passo executável via `nextflow run main.nf` em qualquer máquina Debian com Mamba.
2. **Versionamento atômico** — Toda mudança vai para o GitHub via commits semânticos (`feat:`, `fix:`, `docs:`, `refactor:`, `perf:`).
3. **Economia de tokens** — NUNCA repita conteúdo de arquivos no contexto; use `@filename` para referenciar.
4. **Aprendizado persistente** — Toda decisão/erro/insight vai para `LEARNINGS.md`.
5. **Quality gates obrigatórios** — Nenhuma fase avança sem aprovação do agente revisor.
6. **Validação biológica** — Toda saída técnica passa pelo agente `literature-validator` antes do commit.

### 0.2 Modo de operação

```
LOCAL (Claude Code) ──[ git push ]──> GitHub ──[ webhook/pull ]──> Servidor Debian
        │                                                                │
        └─── desenvolvimento + revisão                    execução pesada (GPU/CPU)
```

- **Local:** desenvolvimento de código Nextflow, scripts, ambientes Mamba, documentação.
- **GitHub:** repositório `trypsin-agemmatalis-structural`, branch `main` (estável), `dev` (trabalho ativo).
- **Servidor Debian:** executa `nextflow run` com GPU para AlphaFold3 e GROMACS.

### 0.3 Antes de qualquer comando

```
[STOP] Antes de executar qualquer comando, verifique:
  1. O CLAUDE.md raiz está carregado?
  2. O LEARNINGS.md foi consultado para erros conhecidos?
  3. O agente apropriado foi invocado?
  4. A mudança é commitável atomicamente?
```

---

## 1. ESTRUTURA DE DIRETÓRIOS A CRIAR

```
trypsin-agemmatalis-structural/
├── CLAUDE.md                      # Contexto persistente para Claude Code
├── LEARNINGS.md                   # Aprendizados, erros e decisões
├── ISSUES.md                      # Bugs e workarounds documentados
├── README.md                      # Documentação do projeto
├── projeto_caracteriz_trypsin.md  # Este arquivo (fonte da verdade)
│
├── .claude/
│   ├── agents/                    # Subagentes especializados
│   │   ├── bioinformatics-architect.md
│   │   ├── nextflow-developer.md
│   │   ├── mamba-env-manager.md
│   │   ├── code-reviewer.md
│   │   ├── debugger.md
│   │   ├── literature-validator.md
│   │   ├── publication-formatter.md
│   │   └── token-optimizer.md
│   │
│   ├── skills/                    # Habilidades reutilizáveis
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
│   ├── hooks/                     # Hooks de validação automática
│   │   ├── pre-commit.sh
│   │   ├── pre-push.sh
│   │   └── post-execute.sh
│   │
│   └── commands/                  # Comandos slash customizados
│       ├── /run-phase.md
│       ├── /validate.md
│       ├── /sync-server.md
│       └── /compact-context.md
│
├── nextflow/
│   ├── main.nf                    # Pipeline principal
│   ├── nextflow.config            # Configs (local, debian, gpu)
│   ├── conf/
│   │   ├── base.config
│   │   ├── debian.config          # Perfil servidor
│   │   ├── gpu.config             # Para AF3 e GROMACS
│   │   └── test.config
│   │
│   ├── modules/                   # Módulos atômicos
│   │   ├── local/
│   │   │   ├── busco.nf
│   │   │   ├── cd_hit.nf
│   │   │   ├── transdecoder.nf
│   │   │   ├── diamond.nf
│   │   │   ├── hmmer.nf
│   │   │   ├── interproscan.nf
│   │   │   ├── signalp.nf
│   │   │   ├── mafft.nf
│   │   │   ├── trimal.nf
│   │   │   ├── iqtree.nf
│   │   │   ├── alphafold3.nf
│   │   │   ├── foldseek.nf
│   │   │   ├── tmalign.nf
│   │   │   ├── consurf.nf
│   │   │   ├── autodock_vina.nf
│   │   │   ├── haddock.nf
│   │   │   ├── gromacs.nf
│   │   │   ├── molprobity.nf
│   │   │   └── prosa.nf
│   │   └── nf-core/
│   │       └── proteinfold/       # Submódulo nf-core
│   │
│   ├── subworkflows/
│   │   ├── qc_assembly.nf         # Fase 1
│   │   ├── orf_prediction.nf      # Fase 2
│   │   ├── trypsin_id.nf          # Fase 3
│   │   ├── completeness.nf        # Fase 4
│   │   ├── primary_char.nf        # Fase 5
│   │   ├── phylogeny.nf           # Fase 6
│   │   ├── structure_pred.nf      # Fase 7
│   │   ├── struct_validation.nf   # Fase 8
│   │   ├── docking.nf             # Fase 9
│   │   ├── md_simulation.nf       # Fase 10
│   │   └── report_gen.nf          # Fase 11
│   │
│   └── bin/                       # Scripts auxiliares
│       ├── filter_complete_trypsins.py
│       ├── extract_catalytic_triad.py
│       ├── parse_alphafold3.py
│       ├── compute_protparam.py
│       └── make_figures.py
│
├── envs/                          # Ambientes Mamba
│   ├── assembly_qc.yml
│   ├── orf_prediction.yml
│   ├── annotation.yml
│   ├── phylogeny.yml
│   ├── structure.yml
│   ├── docking.yml
│   ├── md.yml
│   └── analysis.yml
│
├── data/
│   ├── raw/
│   │   └── trinity_assembly.fasta  # Input do usuário
│   ├── references/
│   │   ├── uniprot_sprot.fasta
│   │   ├── Pfam-A.hmm
│   │   ├── insecta_odb10/
│   │   └── trypsin_refs_lepidoptera.fasta
│   └── intermediate/               # Outputs de cada fase
│
├── results/
│   ├── 01_qc/
│   ├── 02_orfs/
│   ├── 03_trypsin_ids/
│   ├── 04_complete/
│   ├── 05_primary/
│   ├── 06_phylogeny/
│   ├── 07_structures/
│   ├── 08_validation/
│   ├── 09_docking/
│   ├── 10_md/
│   └── 11_report/
│
├── manuscript/
│   ├── manuscript.md
│   ├── figures/
│   ├── tables/
│   └── supplementary/
│
└── tests/
    ├── unit/
    ├── integration/
    └── data/                       # Mini dataset para testes
```

---

## 2. AGENTES (SUBAGENTS) — DEFINIÇÃO COMPLETA

> Cada agente vive em `.claude/agents/<nome>.md`. Use `Task` tool com `subagent_type: "<nome>"`.

### 2.1 `bioinformatics-architect`

```yaml
---
name: bioinformatics-architect
description: Arquiteto chefe do pipeline. Use PROATIVAMENTE no início de cada fase e quando houver decisão de design. Avalia o "porquê" antes do "como".
tools: Read, Grep, Glob, WebSearch, WebFetch
---

Você é um arquiteto bioinformático sênior com 15+ anos em pipelines de proteômica
estrutural e RNA-Seq. Sua única missão: garantir que cada fase do pipeline esteja
cientificamente justificada e use o método STATE-OF-THE-ART (2025-2026).

Antes de aprovar qualquer fase, verifique:
1. O método é o mais atual? (cite paper de 2024-2026)
2. Há alternativa melhor? (compare ≥2 opções)
3. Os parâmetros estão calibrados para Lepidoptera/insetos?
4. A saída é interpretável biologicamente?

Saída obrigatória: bloco JSON com {decisao, justificativa, referencias, parametros, alternativas_consideradas}.
NÃO escreva código. Apenas decida e justifique.
```

### 2.2 `nextflow-developer`

```yaml
---
name: nextflow-developer
description: Escreve módulos Nextflow DSL2 limpos. Use SEMPRE que precisar criar/editar arquivos .nf. NUNCA escreva Nextflow fora deste agente.
tools: Read, Write, Edit, Bash, Grep
---

Você é especialista em Nextflow DSL2 e nf-core. Regras inegociáveis:

1. SEMPRE use DSL2 (`nextflow.enable.dsl=2`).
2. Cada processo tem: `tag`, `publishDir`, `conda` (env file), `cpus`, `memory`, `time`.
3. Inputs/outputs com `tuple val(meta), path(file)`.
4. Use `Channel.fromPath` com `.checkIfExists()`.
5. `errorStrategy 'retry'`, `maxRetries 2` em processos pesados.
6. Documentação inline (// comentários explicando o porquê).
7. Toda variável de config vem de `params.<x>`, nunca hardcoded.

Antes de finalizar, RODE: `nextflow run -stub-run -profile test`.
Se falhar, ITERE até passar. NUNCA entregue código não testado.
```

### 2.3 `mamba-env-manager`

```yaml
---
name: mamba-env-manager
description: Cria e mantém ambientes Conda/Mamba (.yml). Use quando ferramenta nova for adicionada ao pipeline.
tools: Read, Write, Edit, Bash
---

Especialista em resolução de dependências bioinformáticas. Regras:

1. UM ambiente por subworkflow (não 1 mega-env).
2. Sempre fixe versões: `bioconda::trinity=2.15.1` (NUNCA `latest`).
3. Canais na ordem: `conda-forge`, `bioconda`, `defaults`.
4. Teste o env: `mamba env create -f env.yml --dry-run`.
5. Documente conflitos conhecidos em LEARNINGS.md.
6. Para ferramentas GPU (AlphaFold3, GROMACS), separe em env próprio.

Saída: arquivo .yml + comando de teste + tempo estimado de criação.
```

### 2.4 `code-reviewer`

```yaml
---
name: code-reviewer
description: Revisor crítico de TODO código antes de commit. Use PROATIVAMENTE após qualquer escrita de código. BLOQUEIA push se falhar.
tools: Read, Grep, Glob, Bash
---

Revisor sênior. Sua função é REJEITAR código ruim. Checklist obrigatório:

[ ] Sintaxe válida (rode linter: nf-core lint, pylint, shellcheck)
[ ] Sem secrets/credenciais hardcoded
[ ] Sem caminhos absolutos (use ${baseDir} ou ${projectDir})
[ ] Tratamento de erros explícito
[ ] Logging adequado (stdout + .log)
[ ] Idempotência (rodar 2x produz o mesmo resultado)
[ ] Limites de recurso definidos (CPU, RAM, tempo)
[ ] Compatível com Debian 12+ (kernel 6.x)

Se 1 item falhar: bloqueie com mensagem específica. NÃO arrume você mesmo —
delegue de volta ao desenvolvedor com instruções precisas.
```

### 2.5 `debugger`

```yaml
---
name: debugger
description: Acionado QUANDO erro ocorre. Não use para escrever código novo. Apenas diagnóstico + correção mínima.
tools: Read, Bash, Grep, Edit, WebSearch
---

Debugger metódico. Processo OBRIGATÓRIO:

1. Reproduza o erro (mostre comando + saída).
2. Isole a causa raiz (NÃO o sintoma).
3. Consulte LEARNINGS.md — esse erro já apareceu?
4. Se SIM: aplique solução conhecida e adicione contador.
5. Se NÃO: hipóteses ranqueadas (3 no máximo), teste 1 por 1.
6. Após corrigir, ADICIONE ao LEARNINGS.md em formato:
   ```
   ## ERRO-NNN [data]
   - Sintoma: ...
   - Causa raiz: ...
   - Solução: ...
   - Como evitar: ...
   ```

NUNCA aplique "fix" sem entender a causa. NUNCA suprima warnings.
```

### 2.6 `literature-validator`

```yaml
---
name: literature-validator
description: Valida cientificamente cada saída biológica. Acione após cada fase para confirmar plausibilidade.
tools: Read, WebSearch, WebFetch, Grep
---

Validador científico. Para cada resultado biológico, responda:

1. O número de tripsinas encontradas é plausível? (Lepidoptera: 5-30 isoformas típicas)
2. Os pesos moleculares estão no range esperado? (20-35 kDa)
3. A tríade catalítica está presente nos modelos? (His57, Asp102, Ser195)
4. O pI está em range de proteínas de midgut alcalino?
5. Os resultados batem com papers de A. gemmatalis ou parentes (Spodoptera, Helicoverpa)?

Cite SEMPRE 2-3 papers de 2023-2026 como referência.
Se algo destoar, FLAGUE com severidade (info/warning/critical).
```

### 2.7 `publication-formatter`

```yaml
---
name: publication-formatter
description: Prepara figuras, tabelas e seções do manuscrito em formato pronto para submissão.
tools: Read, Write, Edit, Bash
---

Especialista em padrões editoriais (Nature, Cell, IJBM, IBMB). Regras:

1. Figuras: 300 DPI, formato vetorial (SVG/PDF), Arial 8-12pt, escala em barras.
2. Tabelas: CSV → LaTeX ou Markdown estendido.
3. Métodos: passivo, passo a passo, com versões de softwares + DOIs.
4. Resultados: sentenças declarativas com estatística (n, p, IC).
5. Sempre cite o pipeline GitHub no Data Availability.

Saídas: arquivos prontos em /manuscript/.
```

### 2.8 `token-optimizer`

```yaml
---
name: token-optimizer
description: Acione QUANDO contexto exceder 60% da janela. Faz compactação cirúrgica preservando estado essencial.
tools: Read, Write, Edit
---

Otimizador de contexto. Sua tarefa: reduzir contexto SEM perder informação crítica.

Estratégia:
1. Identifique blocos repetidos → substitua por referência (@file:linha).
2. Resuma outputs de comandos antigos em 1 linha.
3. Mova decisões para LEARNINGS.md e referencie.
4. Preserve INTEGRAL: CLAUDE.md, fase atual, último erro, próximos 3 passos.
5. NUNCA descarte: parâmetros experimentais, paths de dados, IDs de amostras.

Saída: relatório do que foi compactado + economia em tokens estimada.
```

---

## 3. SKILLS (HABILIDADES) — REUTILIZÁVEIS

> Cada skill vive em `.claude/skills/<nome>/SKILL.md`. São invocadas automaticamente quando o contexto bate.

### 3.1 Template padrão de SKILL.md

```markdown
---
name: <nome-da-skill>
description: <quando-usar, com palavras-chave de trigger>
---

# <Nome>

## Quando usar
- Trigger 1
- Trigger 2

## Inputs esperados
- arquivo X (formato Y)

## Procedimento
1. Validar input
2. Executar comando
3. Validar output
4. Atualizar LEARNINGS.md

## Outputs garantidos
- arquivo Z em /results/...

## Erros comuns e correções
| Erro | Causa | Fix |
|---|---|---|
| ... | ... | ... |
```

### 3.2 Skills a criar (lista completa)

| Skill | Função | Fase |
|---|---|---|
| `trinity-parser` | Lê e valida FASTA do Trinity, extrai estatísticas | 1 |
| `transdecoder-runner` | ORF prediction com suporte BLAST+Pfam | 2 |
| `trypsin-identifier` | DIAMOND + HMMER (Pfam PF00089) com interseção | 3 |
| `completeness-filter` | Tríade catalítica + Met inicial + ≥220aa + cobertura | 4 |
| `alphafold3-runner` | Wrap nf-core/proteinfold modo AF3 | 7 |
| `foldseek-search` | Busca estrutural contra AFDB + PDB | 8 |
| `docking-pipeline` | AutoDock Vina + HADDOCK + PLIP analysis | 9 |
| `gromacs-md` | 100ns MD em complex tripsina-inibidor com análise RMSD/RMSF/SASA | 10 |
| `figure-generator` | Matplotlib/PyMOL/seaborn → figuras 300dpi | 11 |
| `report-writer` | Geração de seções Methods/Results | 11 |

### 3.3 Exemplo completo: `.claude/skills/trypsin-identifier/SKILL.md`

```markdown
---
name: trypsin-identifier
description: Identifica tripsinas em proteoma predito via dupla validação DIAMOND + HMMER (Pfam Tryp_SPc PF00089). Use quando houver arquivo .pep do TransDecoder.
---

# Trypsin Identifier (DIAMOND + HMMER)

## Quando usar
- Após TransDecoder produzir `*.pep`
- Quando usuário pede "identificar tripsinas"
- Para qualquer screen de serino-proteases

## Inputs
- `proteome.pep` (TransDecoder output)
- `uniprot_sprot.dmnd` (DIAMOND DB)
- `Pfam-A.hmm` (com PF00089 incluído)

## Procedimento

1. **DIAMOND BLASTp** (mais rápido que BLAST tradicional, padrão 2025):
   ```bash
   diamond blastp \
     --query proteome.pep \
     --db uniprot_sprot.dmnd \
     --outfmt 6 qseqid sseqid stitle pident length qlen slen evalue bitscore \
     --evalue 1e-10 --max-target-seqs 5 \
     --threads ${task.cpus} \
     --out diamond_hits.tsv
   ```

2. **HMMER scan com PF00089**:
   ```bash
   hmmfetch Pfam-A.hmm Tryp_SPc > trypsin.hmm
   hmmsearch --domtblout trypsin_hmm.out --cpu ${task.cpus} \
             -E 1e-10 trypsin.hmm proteome.pep
   ```

3. **Interseção robusta** (Python):
   ```python
   diamond_ids = set(filtered_diamond_hits_with_trypsin_keyword)
   hmm_ids = set(hmm_hits)
   confident_trypsins = diamond_ids & hmm_ids   # AMBOS métodos concordam
   suggestive = (diamond_ids | hmm_ids) - confident_trypsins  # divergentes
   ```

4. **Salvar dois conjuntos**:
   - `trypsins_confident.fasta` (interseção)
   - `trypsins_suggestive.fasta` (revisão manual)

## Outputs
- `results/03_trypsin_ids/trypsins_confident.fasta`
- `results/03_trypsin_ids/trypsins_suggestive.fasta`
- `results/03_trypsin_ids/identification_report.tsv`

## Erros comuns
| Erro | Causa | Fix |
|---|---|---|
| 0 hits HMMER | Pfam DB desatualizado | Baixe versão ≥36.0 |
| Muitos hits "trypsin-like" | E-value frouxo | Use 1e-15 ao invés de 1e-10 |
| Diamond db error | Versão diamond ≠ db | Recrie .dmnd na mesma versão |
```

---

## 4. HOOKS DE AUTOMAÇÃO

### 4.1 `.claude/hooks/pre-commit.sh`

```bash
#!/usr/bin/env bash
set -euo pipefail

echo "🔍 Pre-commit checks..."

# 1. Lint Nextflow
if find . -name "*.nf" | head -1 > /dev/null; then
  nf-core lint --release || { echo "❌ nf-core lint failed"; exit 1; }
fi

# 2. Lint Python
if find . -name "*.py" | head -1 > /dev/null; then
  ruff check bin/ tests/ || { echo "❌ ruff failed"; exit 1; }
fi

# 3. Validar YML envs
for env in envs/*.yml; do
  python -c "import yaml; yaml.safe_load(open('$env'))" \
    || { echo "❌ YAML inválido: $env"; exit 1; }
done

# 4. Verificar segredos
if grep -rE "(API_KEY|PASSWORD|TOKEN)=[\"'][^\"']+[\"']" \
     --exclude-dir=.git . > /dev/null 2>&1; then
  echo "❌ Possível secret detectado!"; exit 1
fi

echo "✅ Pre-commit OK"
```

### 4.2 `.claude/hooks/pre-push.sh`

```bash
#!/usr/bin/env bash
set -euo pipefail

echo "🚀 Pre-push: rodando teste stub..."
cd nextflow
nextflow run main.nf -profile test -stub-run --outdir results_stub
rm -rf results_stub
echo "✅ Stub-run passou"
```

### 4.3 `.claude/hooks/post-execute.sh`

```bash
#!/usr/bin/env bash
# Acionado após cada fase concluída
PHASE=$1
echo "📊 Validando fase ${PHASE}..."

# Invoca literature-validator via Claude Code
# (placeholder — Claude Code chama isso automaticamente)
echo "Phase ${PHASE} completed at $(date)" >> LEARNINGS.md
```

---

## 5. COMANDOS SLASH CUSTOMIZADOS

### 5.1 `/run-phase` — `.claude/commands/run-phase.md`

```markdown
---
description: Executa uma fase específica do pipeline com validação
---

Execute a fase $ARGUMENTS do pipeline.

Procedimento:
1. Invoque `bioinformatics-architect` para revisar parâmetros da fase
2. Invoque `nextflow-developer` para garantir que módulos estão prontos
3. Rode: `nextflow run main.nf -profile debian --phase $ARGUMENTS -resume`
4. Após conclusão, invoque `literature-validator`
5. Se aprovado, commit: `git commit -m "feat(phase-$ARGUMENTS): completed"`
6. Push para `dev`
```

### 5.2 `/sync-server`

```markdown
---
description: Sincroniza local → GitHub → Servidor Debian
---

1. `git status` — confirme mudanças
2. Invoque `code-reviewer` se houver código novo
3. `git push origin dev`
4. SSH no servidor: `ssh ${SERVER_USER}@${SERVER_HOST}`
5. `cd ~/trypsin-agemmatalis-structural && git pull`
6. Reporte hash do commit deployado.
```

### 5.3 `/compact-context`

```markdown
---
description: Compacta o contexto preservando estado essencial
---

Invoque `token-optimizer` com instrução de reduzir contexto em 50%
preservando: fase atual, LEARNINGS.md último entry, próximos 3 passos.
```

### 5.4 `/validate`

```markdown
---
description: Roda validação completa: lint + stub-run + literature check
---

1. `.claude/hooks/pre-commit.sh`
2. `.claude/hooks/pre-push.sh`
3. Invoque `literature-validator` com último resultado em /results
4. Reporte status final.
```

---

## 6. PIPELINE NEXTFLOW — ESQUELETO

### 6.1 `nextflow/nextflow.config`

```groovy
// nextflow.config
manifest {
    name            = 'trypsin-agemmatalis-structural'
    author          = 'Eulálio Gutemberg'
    homePage        = 'https://github.com/<user>/trypsin-agemmatalis-structural'
    description     = 'Structural characterization of Anticarsia gemmatalis trypsins'
    mainScript      = 'main.nf'
    nextflowVersion = '>=24.04.0'
    version         = '1.0.0'
}

params {
    // Inputs
    input_fasta       = "${projectDir}/data/raw/trinity_assembly.fasta"
    outdir            = "${projectDir}/results"

    // References
    uniprot_db        = "${projectDir}/data/references/uniprot_sprot.dmnd"
    pfam_db           = "${projectDir}/data/references/Pfam-A.hmm"
    busco_lineage     = "${projectDir}/data/references/insecta_odb10"
    lep_refs          = "${projectDir}/data/references/trypsin_refs_lepidoptera.fasta"

    // Tool thresholds
    diamond_evalue    = 1e-10
    hmmer_evalue      = 1e-10
    min_aa_length     = 220
    blast_coverage    = 0.80
    blast_identity    = 30

    // Phases (1-11; null = run all)
    phase             = null

    // Compute
    max_cpus          = 32
    max_memory        = '128.GB'
    max_time          = '72.h'
    use_gpu           = false
}

profiles {
    test    { includeConfig 'conf/test.config' }
    debian  { includeConfig 'conf/debian.config' }
    gpu     { includeConfig 'conf/gpu.config' }
}

// Conda enabled by default
conda.enabled = true
mamba.enabled = true

// Tracing
timeline.enabled = true
timeline.file    = "${params.outdir}/pipeline_info/timeline.html"
report.enabled   = true
report.file      = "${params.outdir}/pipeline_info/report.html"
trace.enabled    = true
trace.file       = "${params.outdir}/pipeline_info/trace.txt"
dag.enabled      = true
dag.file         = "${params.outdir}/pipeline_info/dag.svg"
```

### 6.2 `nextflow/main.nf` (esqueleto)

```groovy
#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Subworkflows
include { QC_ASSEMBLY      } from './subworkflows/qc_assembly.nf'
include { ORF_PREDICTION   } from './subworkflows/orf_prediction.nf'
include { TRYPSIN_ID       } from './subworkflows/trypsin_id.nf'
include { COMPLETENESS     } from './subworkflows/completeness.nf'
include { PRIMARY_CHAR     } from './subworkflows/primary_char.nf'
include { PHYLOGENY        } from './subworkflows/phylogeny.nf'
include { STRUCTURE_PRED   } from './subworkflows/structure_pred.nf'
include { STRUCT_VALIDATION} from './subworkflows/struct_validation.nf'
include { DOCKING          } from './subworkflows/docking.nf'
include { MD_SIMULATION    } from './subworkflows/md_simulation.nf'
include { REPORT_GEN       } from './subworkflows/report_gen.nf'

workflow {
    ch_input = Channel.fromPath(params.input_fasta, checkIfExists: true)

    // Cada fase recebe a saída da anterior
    QC_ASSEMBLY      ( ch_input )
    ORF_PREDICTION   ( QC_ASSEMBLY.out.clean_fasta )
    TRYPSIN_ID       ( ORF_PREDICTION.out.pep )
    COMPLETENESS     ( TRYPSIN_ID.out.trypsins_confident )
    PRIMARY_CHAR     ( COMPLETENESS.out.complete_trypsins )
    PHYLOGENY        ( COMPLETENESS.out.complete_trypsins, params.lep_refs )
    STRUCTURE_PRED   ( COMPLETENESS.out.complete_trypsins )
    STRUCT_VALIDATION( STRUCTURE_PRED.out.pdbs )
    DOCKING          ( STRUCT_VALIDATION.out.validated_pdbs )
    MD_SIMULATION    ( DOCKING.out.top_complexes )
    REPORT_GEN       (
        QC_ASSEMBLY.out, ORF_PREDICTION.out, TRYPSIN_ID.out,
        COMPLETENESS.out, PRIMARY_CHAR.out, PHYLOGENY.out,
        STRUCTURE_PRED.out, STRUCT_VALIDATION.out,
        DOCKING.out, MD_SIMULATION.out
    )
}

workflow.onComplete {
    log.info "Pipeline finished: ${workflow.success ? 'SUCCESS' : 'FAILED'}"
    log.info "Duration: ${workflow.duration}"
    log.info "Output: ${params.outdir}"
}
```

### 6.3 Exemplo de módulo: `modules/local/alphafold3.nf`

```groovy
process ALPHAFOLD3 {
    tag        "${meta.id}"
    label      'process_gpu'
    publishDir "${params.outdir}/07_structures/alphafold3", mode: 'copy'

    conda 'envs/structure.yml'

    cpus       8
    memory     '64.GB'
    time       '24.h'
    accelerator 1, type: 'nvidia-tesla-a100'

    errorStrategy 'retry'
    maxRetries    2

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}_model_*.pdb"), emit: pdb
    tuple val(meta), path("${meta.id}_confidence.json"), emit: confidence
    tuple val(meta), path("${meta.id}.log"), emit: log

    script:
    """
    # Roda nf-core/proteinfold em modo AF3 single
    run_alphafold3.py \\
        --json_path ${fasta} \\
        --output_dir . \\
        --model_dir \${ALPHAFOLD3_MODELS} \\
        --db_dir \${ALPHAFOLD3_DB} \\
        --jackhmmer_n_cpu ${task.cpus} \\
        --num_recycles 10 \\
        2>&1 | tee ${meta.id}.log

    # Renomear outputs
    for pdb in *_model_*.cif; do
        python -c "from Bio.PDB import MMCIFParser, PDBIO; \\
                   p=MMCIFParser().get_structure('m','\$pdb'); \\
                   io=PDBIO(); io.set_structure(p); \\
                   io.save('\${pdb%.cif}.pdb')"
    done
    """

    stub:
    """
    touch ${meta.id}_model_0.pdb
    echo '{"plddt": 85.0, "ptm": 0.8}' > ${meta.id}_confidence.json
    touch ${meta.id}.log
    """
}
```

---

## 7. AMBIENTES MAMBA

### 7.1 `envs/assembly_qc.yml`

```yaml
name: assembly_qc
channels:
  - conda-forge
  - bioconda
dependencies:
  - python=3.11
  - busco=5.7.1
  - cd-hit=4.8.1
  - seqkit=2.8.2
  - bbmap=39.06
```

### 7.2 `envs/structure.yml` (AlphaFold3 + análise)

```yaml
name: structure
channels:
  - conda-forge
  - bioconda
  - nvidia
dependencies:
  - python=3.11
  - cuda-toolkit=12.4
  - pytorch>=2.3
  - jaxlib
  - biopython=1.83
  - pymol-open-source=3.0
  - foldseek=9
  - tmalign=20220412
  - dssp=4.4
  - pip
  - pip:
    - alphafold3   # via repo oficial
    - colabfold
```

### 7.3 `envs/md.yml`

```yaml
name: md
channels:
  - conda-forge
  - bioconda
dependencies:
  - gromacs=2025.1
  - openmm=8.1
  - mdanalysis=2.7
  - mdtraj=1.10
  - python=3.11
  - numpy
  - matplotlib
  - seaborn
```

---

## 8. FASES DETALHADAS DO PIPELINE

| # | Fase | Subworkflow | Tempo estimado | Ambiente |
|---|---|---|---|---|
| 1 | QC do assembly (BUSCO + CD-HIT-EST) | `qc_assembly.nf` | 1-2h | assembly_qc |
| 2 | Predição de ORFs (TransDecoder + Pfam) | `orf_prediction.nf` | 3-6h | orf_prediction |
| 3 | Identificação de tripsinas (DIAMOND + HMMER) | `trypsin_id.nf` | 30min | annotation |
| 4 | Filtro de completude (tríade + Met + ≥220aa) | `completeness.nf` | 10min | annotation |
| 5 | Caracterização primária (ProtParam + SignalP + InterProScan) | `primary_char.nf` | 1h | annotation |
| 6 | Filogenia (MAFFT + trimAl + IQ-TREE2) | `phylogeny.nf` | 2-4h | phylogeny |
| 7 | **Predição estrutural (AlphaFold3 via nf-core/proteinfold)** | `structure_pred.nf` | 12-48h | structure (GPU) |
| 8 | Validação (MolProbity + ProSA + Foldseek + ConSurf) | `struct_validation.nf` | 2h | structure |
| 9 | Docking (AutoDock Vina + HADDOCK + PLIP) | `docking.nf` | 4-8h | docking |
| 10 | **MD simulations 100ns (GROMACS 2025)** | `md_simulation.nf` | 24-72h | md (GPU) |
| 11 | Geração de relatório/figuras/manuscrito | `report_gen.nf` | 1h | analysis |

---

## 9. ESTRATÉGIA DE EXECUÇÃO REMOTA

### 9.1 Fluxo Local → GitHub → Servidor

```bash
# LOCAL (Claude Code)
git checkout dev
# ... desenvolvimento ...
.claude/hooks/pre-commit.sh
git add . && git commit -m "feat(phase-3): trypsin id with diamond"
.claude/hooks/pre-push.sh
git push origin dev

# SERVIDOR DEBIAN (via SSH ou GitHub Actions self-hosted runner)
ssh user@server
cd ~/trypsin-agemmatalis-structural
git pull origin dev

# Criar envs (uma vez)
for env in envs/*.yml; do mamba env create -f $env; done

# Rodar pipeline
nextflow run main.nf \
    -profile debian,gpu \
    --input_fasta data/raw/trinity_assembly.fasta \
    --outdir results \
    -resume \
    -with-report \
    -with-timeline \
    -with-dag dag.svg
```

### 9.2 Configuração `conf/debian.config`

```groovy
params {
    config_profile_name        = 'Debian server'
    config_profile_description = 'Local Debian 12 with 64 cores, 256GB RAM, 2x A100'
    max_cpus    = 64
    max_memory  = '256.GB'
    max_time    = '120.h'
}

process {
    executor       = 'local'
    cpus           = { check_max( 2     * task.attempt, 'cpus'   ) }
    memory         = { check_max( 8.GB  * task.attempt, 'memory' ) }
    time           = { check_max( 4.h   * task.attempt, 'time'   ) }

    withLabel:process_low    { cpus = 4;  memory = '16.GB'; time = '4.h' }
    withLabel:process_medium { cpus = 16; memory = '64.GB'; time = '12.h' }
    withLabel:process_high   { cpus = 32; memory = '128.GB'; time = '48.h' }
    withLabel:process_gpu    {
        cpus = 16; memory = '128.GB'; time = '72.h'
        accelerator = [request: 1, type: 'nvidia-tesla-a100']
    }
}
```

---

## 10. ARQUIVOS DE ESTADO PERSISTENTE

### 10.1 `CLAUDE.md` (raiz)

```markdown
# Projeto: Tripsinas A. gemmatalis — Caracterização Estrutural

## Estado atual
- Fase atual: <preenchido pelo agente>
- Último commit: <hash>
- Bloqueios: <listar>

## Contexto biológico
- Organismo: Anticarsia gemmatalis (Lepidoptera: Erebidae)
- Tecido: midgut (pH alcalino, ~10-11)
- Tripsinas esperadas: 5-30 isoformas funcionais
- Função: digestão de proteínas vegetais (soja)

## Decisões de design
- AlphaFold3 > AF2 (paper InsectBase 3.0, NAR 2026)
- Foldseek obrigatório (Cell Res 2026)
- GROMACS 100ns mínimo para Q1 journals
- Inibidores de teste: SKTI + BPTI-RCL peptides (TGPCK, AVIMK)

## Para Claude Code
SEMPRE leia este arquivo primeiro. Atualize "Estado atual" ao final de cada sessão.
```

### 10.2 `LEARNINGS.md` (template inicial)

```markdown
# Aprendizados, Erros e Decisões

## Convenção
- ERRO-NNN: erros catalogados
- DEC-NNN: decisões de design
- INS-NNN: insights biológicos

---

## DEC-001 [2026-05-25] AlphaFold3 vs AlphaFold2
- Decisão: usar AF3
- Razão: AF3 prediz interações com ligantes (essencial para docking com SKTI)
- Referência: Abramson et al. 2024 Nature; InsectBase 3.0 (NAR 2026)

## DEC-002 [2026-05-25] Filtro de completude
- Critérios: Met inicial + ≥220 aa + tríade His-Asp-Ser por alinhamento + Pfam PF00089 completo
- Referência: revisão IntechOpen 2022 (Trypsins Lepidoptera)
```

### 10.3 `ISSUES.md`

```markdown
# Issues Conhecidos e Workarounds

## Aberto
(nenhum)

## Resolvido
(será populado durante execução)
```

---

## 11. ECONOMIA DE TOKENS — REGRAS RÍGIDAS

1. **Referencie, não copie:** Use `@nextflow/modules/local/alphafold3.nf:42` em vez de colar o arquivo.
2. **Compact após cada fase:** Acione `token-optimizer` no fim de cada fase concluída.
3. **Outputs longos:** Salve em arquivo + reporte só o caminho + 5 linhas chave.
4. **Logs:** Use `tail -50` e `grep ERROR` ao invés de cat completo.
5. **Histórico de erros:** Está no LEARNINGS.md — não relembre, releia.
6. **Mensagens de commit:** Use `git log --oneline -10` para contexto histórico.
7. **Saídas de Nextflow:** Use `.nextflow.log` cirurgicamente (últimas 100 linhas).

---

## 12. CRITÉRIOS DE QUALIDADE DO ARTIGO FINAL

Para garantir IF > 5 (alvo: *Insect Biochem Mol Biol*, *Int J Biol Macromol*, *Comput Struct Biotechnol J*):

| Critério | Threshold |
|---|---|
| Tripsinas completas identificadas | ≥ 8 isoformas |
| pLDDT médio dos modelos AF3 | > 80 |
| Ramachandran favored | > 95% |
| MD simulations | ≥ 100 ns por complexo, ≥ 3 réplicas |
| Inibidores testados | ≥ 3 (SKTI + 2 BPTI-RCL peptides) |
| Análise filogenética | ≥ 30 sequências de Lepidoptera |
| Foldseek hits | TM-score > 0.7 com referências PDB |
| Reprodutibilidade | Pipeline público no GitHub + Zenodo DOI |

---

## 13. ROADMAP DE EXECUÇÃO (PARA CLAUDE CODE)

```
SESSÃO 1 (setup):
  [ ] /init — criar estrutura de diretórios
  [ ] Invocar `bioinformatics-architect` para validar arquitetura
  [ ] Invocar `mamba-env-manager` para gerar todos os envs
  [ ] git init, criar repo no GitHub, primeiro commit

SESSÃO 2 (módulos base):
  [ ] Invocar `nextflow-developer` para fases 1-4 (módulos + subworkflows)
  [ ] /validate
  [ ] git push dev

SESSÃO 3 (módulos avançados):
  [ ] Invocar `nextflow-developer` para fases 5-7
  [ ] Integrar nf-core/proteinfold como submódulo
  [ ] /validate
  [ ] git push dev

SESSÃO 4 (estrutura e docking):
  [ ] Fases 8-10 (validação estrutural, docking, MD)
  [ ] /validate
  [ ] git push dev

SESSÃO 5 (execução remota):
  [ ] /sync-server
  [ ] No servidor: criar envs Mamba
  [ ] Rodar pipeline com dado mínimo (teste)
  [ ] Diagnosticar e corrigir

SESSÃO 6 (execução completa):
  [ ] Rodar pipeline com Trinity FASTA real
  [ ] Validar cada fase com `literature-validator`
  [ ] Iterar parâmetros se necessário

SESSÃO 7 (manuscrito):
  [ ] Invocar `publication-formatter`
  [ ] Gerar figuras + tabelas
  [ ] Redigir Methods + Results
  [ ] Submeter para revista alvo
```

---

## 14. PRIMEIRO COMANDO PARA O CLAUDE CODE

Cole este prompt inicial após colocar este `.md` no diretório:

```
Leia projeto_caracteriz_trypsin.md INTEIRO. Após ler:

1. Confirme que entendeu o objetivo, arquitetura e fases.
2. Crie a estrutura de diretórios da seção 1.
3. Crie todos os arquivos de agentes em .claude/agents/ (seção 2).
4. Crie skeleton de skills em .claude/skills/ (seção 3).
5. Crie hooks em .claude/hooks/ (seção 4) e dê chmod +x.
6. Crie commands em .claude/commands/ (seção 5).
7. Crie CLAUDE.md, LEARNINGS.md, ISSUES.md (seção 10).
8. Crie nextflow.config + main.nf esqueleto (seção 6).
9. Crie todos os envs/*.yml (seção 7).
10. git init, primeiro commit "chore: project scaffold", criar repo GitHub.

Use o agente `bioinformatics-architect` para revisar a estrutura ANTES de criar.
Use o agente `code-reviewer` ANTES de qualquer commit.

Reporte status no final em formato de checklist.
```

---

## 15. NOTAS FINAIS

- **Backup:** Configure `cron` no servidor Debian para snapshot diário de `/results`.
- **Citação:** Reserve DOI no Zenodo ao publicar v1.0.
- **Reprodutibilidade:** Pin TODAS as versões de software no `manifest`.
- **Ética:** Declare uso de IA (Claude Code + AlphaFold3) no manuscrito.

**FIM DO ARQUIVO DE PROJETO**

> Este documento é a fonte canônica. Qualquer divergência entre código e este arquivo deve ser resolvida atualizando AMBOS conscientemente.
