---
name: trypsin-identifier
description: Identifica tripsinas em proteoma predito via dupla validação DIAMOND + hmmsearch (PF00089 extraído). Trigger: "identificar tripsinas", "DIAMOND HMMER", após TransDecoder produzir .pep.
---

# Trypsin Identifier (DIAMOND + HMMER) — Fase 3

## Resultado obtido em A. gemmatalis (2026-05-29)
- Input: 17.923 proteínas preditas (TransDecoder)
- DIAMOND vs NR (keywords): **285 candidatos**
- hmmsearch vs PF00089 (cobertura ≥80%): **10.675 candidatos**
- **Confident (∩): 191** — range esperado para Lepidoptera ✅
- Suggestive (∪\∩): 10.578

## CRÍTICO: hmmsearch ≠ hmmscan (colunas diferentes)
Sempre usar **hmmsearch** (query=HMM, target=proteínas). Colunas do domtblout:
```
p[0]  = target name  → ID da proteína ← usar como seq_id
p[3]  = query name   → HMM model (PF00089)
p[5]  = qlen         → comprimento do perfil HMM ← usar para cobertura
p[11] = c-Evalue     → e-value do domínio
p[15] = hmm_from     → início no perfil HMM ← usar para cobertura
p[16] = hmm_to       → fim no perfil HMM   ← usar para cobertura
```
**hmmscan** (query=proteína, target=HMM): p[0]=HMM, p[2]=seq_id, p[17-18]=coords.
Confundir os dois resulta em 0 confident (ERRO-002 no LEARNINGS.md).

## Procedimento correto
```bash
conda activate orf_prediction
bash scripts/phase3_trypsin_id/run.sh
```

O script faz automaticamente:
1. DIAMOND blastp vs NR (`--sensitive -e 1e-10 -k 5 --query-cover 50`)
2. `hmmfetch Pfam-A.hmm PF00089` → extrai modelo isolado
3. `hmmsearch --domtblout` (rápido: minutos vs horas do hmmscan completo)
4. Filtro de cobertura ≥80% com colunas corretas para hmmsearch
5. Interseção Python → confident + suggestive FASTAs

## Bancos disponíveis neste servidor
- NR DIAMOND: `/home/eulalio/databases/nr/nr.dmnd`
- Pfam-A (pressionado): `/home/eulalio/databases/pfam/Pfam-A.hmm`

## Outputs
- `results/phase3/trypsins_confident.fasta` — DIAMOND∩HMMER (191 seqs)
- `results/phase3/trypsins_suggestive.fasta` — só um critério (10.578)
- `results/phase3/identification_report.tsv` — tabela completa
- `results/phase3/hmmer_tryp_spc.domtblout` — domtblout completo do hmmsearch

## Erros comuns
| Erro | Causa | Fix |
|---|---|---|
| 0 confident | Colunas hmmscan usadas para hmmsearch | Corrigir p[0] e p[15-16] |
| hmmsearch lento | Usando Pfam-A.hmm completo | Extrair PF00089 com hmmfetch |
| mamba run falha | exec -- incompatível com bash | Usar conda activate no início do script |
