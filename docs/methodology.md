# Metodologia — Caracterização Estrutural de Tripsinas de *Anticarsia gemmatalis*

> Documento atualizado progressivamente a cada fase concluída.
> Última atualização: Fase 1 (2026-05-28)

---

## 2. Material e Métodos

### 2.1 Material biológico e dados transcriptômicos

O transcriptoma de larvas de *Anticarsia gemmatalis* Hübner, 1818 (Lepidoptera: Erebidae) foi obtido a partir de RNA total extraído do intestino médio (midgut) de lagartas em 5º instar alimentadas com folhas de soja (*Glycine max*). O sequenciamento foi realizado em plataforma Illumina e as leituras foram montadas de novo com o programa Trinity (Grabherr et al., 2011), resultando em 42.372 transcritos.

### 2.2 Avaliação de qualidade do transcriptoma (Fase 1)

#### 2.2.1 Completude do assembly — BUSCO

A completude do transcriptoma foi avaliada com o programa BUSCO v6.0.0 (Manni et al., 2021), utilizando o conjunto de genes conservados de insetos (insecta_odb10, n = 1.367 BUSCOs, atualização 2024-01-08). A análise foi executada no modo transcriptoma (`--mode transcriptome`), com busca por exons via MetaEuk (Levy Karin et al., 2020) seguida de validação por HMMER v3.4 (Eddy, 2011).

#### 2.2.2 Remoção de redundância — CD-HIT-EST

Para reduzir isoformas redundantes geradas pela montagem de novo, os transcritos foram agrupados com CD-HIT-EST v4.8.1 (Li & Godzik, 2006; Fu et al., 2012) com limiar de identidade nucleotídica de 95% (`-c 0.95`), tamanho de palavra `-n 10`, utilizando 16 núcleos de processamento (`-T 16`) e 64 GB de memória (`-M 64000`). Os representantes de cada cluster foram utilizados nas etapas subsequentes.

---

<!-- FASE 2 — TransDecoder (a preencher) -->
<!-- FASE 3 — DIAMOND + HMMER (a preencher) -->
<!-- FASE 4 — Completude da tríade catalítica (a preencher) -->
<!-- FASE 5 — Caracterização primária (a preencher) -->
<!-- FASE 6 — Filogenia (a preencher) -->
<!-- FASE 7 — AlphaFold3 (a preencher) -->
<!-- FASE 8 — Validação estrutural (a preencher) -->
<!-- FASE 9 — Docking molecular (a preencher) -->
<!-- FASE 10 — Dinâmica molecular (a preencher) -->

---

## Referências (parcial — será completado)

- Eddy SR (2011) Accelerated profile HMM searches. *PLoS Comput Biol* 7:e1002195.
- Fu L et al. (2012) CD-HIT: accelerated for clustering the next-generation sequencing data. *Bioinformatics* 28:3150–3152.
- Grabherr MG et al. (2011) Full-length transcriptome assembly from RNA-Seq data without a reference genome. *Nat Biotechnol* 29:644–652.
- Levy Karin E et al. (2020) MetaEuk—sensitive, high-throughput gene discovery, and annotation for large-scale eukaryotic metagenomics. *Microbiome* 8:48.
- Li W, Godzik A (2006) Cd-hit: a fast program for clustering and comparing large sets of protein or nucleotide sequences. *Bioinformatics* 22:1658–1659.
- Manni M et al. (2021) BUSCO update: novel and streamlined workflows along with broader and deeper phylogenetic coverage for scoring of eukaryotic, prokaryotic, and viral genomes. *Mol Biol Evol* 38:4647–4654.
