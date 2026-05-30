---
name: completeness-filter
description: Filtra tripsinas completas por motivo Ser195 (GDSGG/GDXGG) + Met inicial + ≥200 aa. Detecta variantes com Ser195 possivelmente mutado. Trigger: "filtrar tripsinas completas", "tríade catalítica", "sequências completas".
---

# Completeness Filter — Fase 4

## Resultado obtido em A. gemmatalis (2026-05-30)
- Input: 191 confident (DIAMOND∩HMMER)
- **67 completas** (Ser195 intacto, Met, ≥200 aa) → modelagem estrutural
- **6 variantes Ser195** (GDXGG, X≠S) → análise evolutiva
- 55 borderline (Met+comprimento, sem motivo)
- 63 rejeitadas (sem Met ou curtas)
- Comprimento médio (completas): 314 aa (range 209–676 aa)

## Critérios (todos obrigatórios — AND lógico)
1. **Met inicial:** ORF começa com Metionina
2. **Comprimento ≥ 200 aa:** cobre pré-proteína com propeptídeo
3. **Motivo Ser195** (um dos três):
   - Canônico: `GDS[AG]G`
   - Alt conservado: `[GA]DS[GASTVNC]G`
   - Mutante (Ser195→X): `GD[ACTVNILMF]GG` — retido separadamente

## Script de execução
```bash
conda activate orf_prediction
bash scripts/phase4_completeness/run.sh
```

## Outputs
- `results/phase4/complete_trypsins.fasta` — 67 completas (usar para Fase 5+)
- `results/phase4/mutant_ser195_trypsins.fasta` — 6 variantes Ser195
- `results/phase4/borderline_trypsins.fasta` — 55 borderline
- `results/phase4/completeness_report.tsv` — critério por sequência

## Atenção
- Sequências com >500 aa (DN5998: 601 aa, DN11972: 607 aa) são suspeitas de quimera — verificar manualmente antes de modelagem
- Min-len alterado de 220→200 aa para capturar isoformas menores funcionais de Lepidoptera
- A tríade completa His57-Asp102-Ser195 via alinhamento MAFFT ainda não foi verificada (requer `data/references/bovine_trypsin_1TGN.fasta`)

## Expectativa biológica
- *A. gemmatalis*: 5–30 tripsinas distintas após clustering proteico (Fase 5)
- 67 pré-cluster é normal para Trinity (isoformas + variantes alélicas)
- Se < 5: reduzir --min-len para 180 ou incluir suggestive
