---
description: Executa uma fase específica do pipeline com validação completa
---

Execute a fase $ARGUMENTS do pipeline de tripsinas de *A. gemmatalis*.

## Procedimento obrigatório

1. **Consulte LEARNINGS.md** para erros conhecidos na fase $ARGUMENTS
2. **Invoque `bioinformatics-architect`** para revisar parâmetros da fase $ARGUMENTS
3. **Verifique pré-requisitos:** a fase anterior concluiu com sucesso?
4. **Valide ambiente Mamba** correspondente à fase
5. **Execute:**
   ```bash
   nextflow run nextflow/main.nf \
       -profile debian \
       --phase $ARGUMENTS \
       --input_fasta data/raw/trinity_assembly.fasta \
       --outdir results \
       -resume \
       -with-report results/pipeline_info/report_phase$ARGUMENTS.html
   ```
6. **Monitore:** `tail -f .nextflow.log`
7. **Após conclusão:** invoque `literature-validator` com resultados da fase
8. **Se aprovado:** commit:
   ```bash
   git add results/ LEARNINGS.md
   git commit -m "feat(phase-$ARGUMENTS): completed with N outputs"
   git push origin dev
   ```
9. **Execute post-execute hook:**
   ```bash
   .claude/hooks/post-execute.sh $ARGUMENTS
   ```

## Fases disponíveis
| Fase | Nome | Ambiente |
|---|---|---|
| 1 | QC Assembly (BUSCO + CD-HIT-EST) | assembly_qc |
| 2 | ORF Prediction (TransDecoder) | orf_prediction |
| 3 | Trypsin ID (DIAMOND + HMMER) | annotation |
| 4 | Completeness Filter | annotation |
| 5 | Primary Characterization | annotation |
| 6 | Phylogeny (MAFFT + IQ-TREE2) | phylogeny |
| 7 | Structure Prediction (AlphaFold3) | structure [GPU] |
| 8 | Structural Validation | structure |
| 9 | Docking (Vina + HADDOCK) | docking |
| 10 | MD Simulations (GROMACS) | md [GPU] |
| 11 | Report & Manuscript | analysis |
