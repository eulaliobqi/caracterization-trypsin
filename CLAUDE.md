# Projeto: Tripsinas A. gemmatalis — Caracterização Estrutural

## Estado atual
- Fase atual: **SESSÃO 2 INICIADA** — Bugs corrigidos, skills completas (2026-05-25)
- Último commit: feat(session-2): fix HMMER intersection bug + 5 missing skills + README
- Branch ativa: `dev` — GitHub: https://github.com/eulaliobqi/caracterization-trypsin
- Assembly: `data/raw/trinity_assembly.fasta` (42.372 transcritos TRINITY)
- Bloqueios: nenhum — próximo passo: stub-run de validação no servidor Debian

## Contexto biológico
- Organismo: *Anticarsia gemmatalis* Hübner 1818 (Lepidoptera: Erebidae)
- Tecido alvo: midgut (pH alcalino, ~10-11) — principal sítio de digestão de soja
- Tripsinas esperadas: 5-30 isoformas funcionais (serina-proteases, família S1A)
- Função: digestão de proteínas de reserva de soja (*Glycine max*) — alvo de inibidores SKTI/BPTI
- Tríade catalítica obrigatória: **His57 – Asp102 – Ser195** (numeração quimotripsina)
- Tamanho esperado: 20–35 kDa (220–320 aa como pré-proteína)

## Decisões de design (ver LEARNINGS.md para detalhes)
- **AlphaFold3 > AF2**: AF3 prediz interações com ligantes (essencial para docking SKTI/BPTI)
- **Foldseek obrigatório**: busca estrutural por TM-score contra AFDB + PDB
- **GROMACS 2025.1**: ≥100 ns MD por complexo para Q1 journals; 3 réplicas mínimo
- **Inibidores**: SKTI (soybean Kunitz trypsin inhibitor) + peptídeos BPTI-like (TGPCK, AVIMK)
- **Pipeline Nextflow DSL2**: reprodutibilidade total, execução no servidor Debian

## Para Claude Code
SEMPRE leia este arquivo primeiro.
Atualize "Estado atual" ao final de cada sessão.
Consulte LEARNINGS.md antes de qualquer debug.
Use `bioinformatics-architect` antes de mudar parâmetros.

## Estrutura do pipeline (11 fases)
```
Fase 1:  QC assembly       → BUSCO v5 (insecta_odb10) + CD-HIT-EST (0.95)
Fase 2:  ORFs              → TransDecoder + Pfam + BLAST hint
Fase 3:  ID tripsinas      → DIAMOND (UniProt) ∩ HMMER (PF00089 Tryp_SPc)
Fase 4:  Completude        → His+Asp+Ser (tríade) + Met inicial + ≥220 aa
Fase 5:  Caract. primária  → ProtParam + SignalP6 + InterProScan
Fase 6:  Filogenia         → MAFFT-linsi + trimAl + IQ-TREE2
Fase 7:  Estrutura         → AlphaFold3 via nf-core/proteinfold (GPU)
Fase 8:  Validação struct  → MolProbity + ProSA + Foldseek + ConSurf
Fase 9:  Docking           → AutoDock Vina + HADDOCK + PLIP
Fase 10: MD simulations    → GROMACS 2025.1 100 ns (GPU)
Fase 11: Relatório         → Figuras + Manuscrito
```
