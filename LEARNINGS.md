# Aprendizados, Erros e Decisões

## Convenção de nomenclatura
- `ERRO-NNN [data]` — bugs catalogados e solucionados
- `DEC-NNN [data]` — decisões de design justificadas
- `INS-NNN [data]` — insights biológicos relevantes

---

## DEC-001 [2026-05-25] AlphaFold3 vs AlphaFold2
- **Decisão:** usar AlphaFold3 (AF3)
- **Razão:** AF3 prediz interações proteína-ligante com precisão atômica — essencial para modelar complexos tripsina-inibidor (SKTI, BPTI). AF2 não lida com ligantes pequenos.
- **Referências:**
  - Abramson et al. 2024 *Nature* 630, 493–500 (AF3 paper)
  - InsectBase 3.0, NAR 2026 (estruturas de insetos anotadas com AF3)
- **Alternativa descartada:** ColabFold (AF2-based) — sem suporte nativo a ligantes

## DEC-002 [2026-05-25] Filtro de completude — critérios
- **Critérios obrigatórios (AND):**
  1. Metionina inicial na posição 1 do ORF predito
  2. Comprimento ≥ 220 aminoácidos (zymogen + propeptídeo)
  3. Tríade catalítica His-Asp-Ser detectada por alinhamento com bovine trypsin (P00760)
  4. Domínio Pfam PF00089 (Tryp_SPc) com cobertura ≥ 80% e E-value ≤ 1e-10
- **Referência:** Lima et al. 2022 *IntechOpen* — Trypsins in Lepidoptera digestive systems

## DEC-003 [2026-05-25] DIAMOND vs BLAST para identificação
- **Decisão:** DIAMOND v2.1+ no lugar de BLAST (50–500x mais rápido, mesma sensibilidade)
- **Parâmetros:** `--evalue 1e-10 --max-target-seqs 5 --sensitive`
- **Referência:** Buchfink et al. 2021 *Nat Methods* 18, 366–368

## DEC-004 [2026-05-25] HMMER dupla validação
- **Estratégia:** interseção DIAMOND + HMMER (PF00089)
  - `confident_trypsins` = hits em AMBOS (alta especificidade)
  - `suggestive_trypsins` = hits em apenas um (para revisão manual)
- **Razão:** falsos positivos de peptidases S1 não tripsinas são reduzidos pela exigência do domínio Tryp_SPc

## DEC-005 [2026-05-25] MD simulations — parâmetros base
- **Force field:** AMBER99SB-ILDN (padrão para proteínas, compatível com GROMACS 2025)
- **Duração:** 100 ns mínimo (critério de Q1 journals)
- **Réplicas:** 3 por complexo (com velocidades iniciais diferentes)
- **Inibidores testados:** SKTI (PDB: 1AVU) + peptídeos BPTI-RCL (TGPCK, AVIMK)
- **Análises:** RMSD, RMSF, Rg, SASA, contatos H-bonds na interface
- **Referência:** Souza et al. 2024 *IJBM* (MD de tripsinas de Spodoptera)

## INS-001 [2026-05-25] Tripsinas de Lepidoptera — contexto
- Lepidoptera larval midgut é altamente alcalino (pH 10-11)
- Tripsinas são a principal enzima digestiva em *A. gemmatalis* (praga da soja)
- Resistência a inibidores (SKTI) pode evoluir via diversificação de isoformas
- *A. gemmatalis* tem ~3-5 genes de tripsina no genoma; Trinity pode detectar múltiplas isoformas e splicing alternativo
- Referências:
  - Brioschi et al. 2007 *Insect Biochem Mol Biol* (tripsinas A. gemmatalis)
  - Brito et al. 2001 *J Insect Physiol* (interação com SKTI)

## ERRO-001 [2026-05-28] mamba run incompatível com bash no servidor Debian
- **Sintoma:** `mamba run -n ENV cmd` falha com `exec: --: invalid option`
- **Causa raiz:** mamba gera script temporário com `exec -- cmd`; bash builtin exec não aceita `--`
- **Workaround tentado:** `conda activate`, PATH injection, subshell+exec — todos falharam para Perl scripts (TransDecoder) por shebang quebrado no env isolado
- **Solução definitiva:** instalar todas as ferramentas diretamente no ambiente base (`mamba install -n base pkg`). mamba_run e check_env viram no-ops em config.sh.
- **Prevenção:** neste servidor, NÃO usar envs conda isolados — instalar tudo no base

## INS-002 [2026-05-25] Assembly Trinity — características
- Input: 42.372 transcritos (assembly transcriptômico do midgut larval)
- Formato TRINITY_DN*_c*_g*_i* — múltiplas isoformas por gene
- CD-HIT-EST 0.95 recomendado antes de TransDecoder (reduz redundância)
- BUSCO insecta_odb10 permite avaliar completude da representação gênica
