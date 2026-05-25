---
description: Roda validação completa — lint + stub-run + verificação biológica
---

Execute validação completa do pipeline atual.

## Procedimento

1. **Pre-commit checks:**
   ```bash
   bash .claude/hooks/pre-commit.sh
   ```

2. **Stub-run (sem executar ferramentas):**
   ```bash
   nextflow run nextflow/main.nf \
       -stub-run \
       -profile test \
       --outdir results_stub \
       --input_fasta tests/data/test_assembly.fasta
   rm -rf results_stub work
   ```

3. **Validação biológica (se há resultados):**
   Invoque `literature-validator` com os resultados disponíveis em `results/`

4. **Relatório final:**
   - ✅ Se tudo passou: reporte "VALIDADO — pronto para push"
   - ❌ Se algo falhou: reporte erros específicos com linha/arquivo

## O que é validado
- Sintaxe Nextflow (DSL2, processos, canais)
- Ambientes Mamba (YAML válido, versões fixadas)
- Scripts Python (ruff/pylint)
- Hooks bash (shellcheck se disponível)
- Plausibilidade biológica dos últimos resultados
