# Materials and Methods
> Draft — updated progressively as each phase is completed.
> Last update: Phase 2 in progress (2026-05-28)

---

## 2. Materials and Methods

### 2.1 Biological material and transcriptomic data

The larval midgut transcriptome of *Anticarsia gemmatalis* Hübner, 1818 (Lepidoptera: Erebidae) was obtained from total RNA extracted from the midgut tissue of fifth-instar larvae reared on soybean (*Glycine max*) leaves. Paired-end sequencing was performed on the Illumina platform, and reads were assembled *de novo* using Trinity v2.x (Grabherr et al., 2011), yielding 42,372 transcripts. The assembly is deposited under [accession number — to be confirmed].

### 2.2 Transcriptome quality assessment (Phase 1)

#### 2.2.1 Completeness evaluation — BUSCO

Transcriptome completeness was assessed using BUSCO v6.0.0 (Manni et al., 2021) against the insect-specific gene set (insecta_odb10; n = 1,367 conserved genes; release 2024-01-08). The analysis was run in transcriptome mode (`--mode transcriptome`), with exon prediction performed by MetaEuk v7 (Levy Karin et al., 2020) and model validation by HMMER v3.4 (Eddy, 2011).

#### 2.2.2 Redundancy reduction — CD-HIT-EST

To collapse redundant isoforms inherent to *de novo* Trinity assemblies, transcripts were clustered with CD-HIT-EST v4.8.1 (Li & Godzik, 2006; Fu et al., 2012) at 95% nucleotide identity (`-c 0.95`, `-n 10`), using 16 CPU threads (`-T 16`) and 64 GB RAM (`-M 64000`). Cluster representatives were retained for all downstream analyses.

### 2.3 ORF prediction (Phase 2)

#### 2.3.1 Long ORF extraction — TransDecoder

Candidate open reading frames were extracted from the non-redundant transcript set using TransDecoder v5.7.1 (Haas et al., 2013). In the first step, `TransDecoder.LongOrfs` was used to identify all ORFs encoding ≥ 100 amino acids. To improve the accuracy of final ORF selection, two independent sources of homology evidence were incorporated as scoring hints.

#### 2.3.2 Homology-based hints — DIAMOND blastp

Translated ORF candidates were searched against the NCBI non-redundant protein database (NR; downloaded [date]) using DIAMOND v2.1.9 (Buchfink et al., 2021) in `blastp` mode (`--evalue 1e-5 --max-target-seqs 1`). The resulting tabular output (outfmt 6) was provided to TransDecoder as BLAST hints (`--retain_blastp_hits`), favouring the retention of ORFs with significant homology to known proteins.

#### 2.3.3 Domain-based hints — HMMER/Pfam

ORF candidates were additionally screened against the Pfam-A database (release [version]; Mistry et al., 2021) using `hmmscan` (HMMER v3.4; `--cpu 16 --domtblout`). Domain hits were provided as Pfam-based scoring hints (`--retain_pfam_hits`), prioritising ORFs containing known protein domains.

#### 2.3.4 ORF selection — TransDecoder.Predict

Final ORF predictions were generated with `TransDecoder.Predict` using the `--single_best_only` flag, which retains a single optimal ORF per transcript to minimise redundancy in the predicted proteome. The output proteome (`.pep` file) was used for all subsequent analyses.

---

### 2.4 Trypsin identification — dual validation strategy (Phase 3)

#### 2.4.1 Sequence homology — DIAMOND blastp

The predicted proteome (17,923 sequences) was searched against the NCBI non-redundant protein database (NR; DIAMOND v2.1.9; Buchfink et al., 2021) using `blastp` in sensitive mode (`--sensitive --evalue 1e-10 --query-cover 50 --max-target-seqs 5 --outfmt 6`). Hits whose description contained the keywords *trypsin*, *serine protease*, *chymotrypsin*, *trypsinogen*, or *serine endopeptidase* (case-insensitive) were retained as DIAMOND candidates.

#### 2.4.2 Domain detection — HMMER hmmsearch vs PF00089

The Tryp_SPc domain model (PF00089) was extracted from Pfam-A (release [version]; Mistry et al., 2021) using `hmmfetch`, and the predicted proteome was searched with `hmmsearch` (HMMER v3.4; `--cpu 16 -E 1e-10`). Hits with domain coverage ≥ 80% of the PF00089 profile length were retained as HMMER candidates. This approach — searching with a single extracted HMM rather than scanning against all Pfam — substantially reduces computation time while maintaining sensitivity for the target domain.

#### 2.4.3 Dual-validation intersection

Trypsin candidates were classified into two tiers based on the overlap between DIAMOND and HMMER evidence: (i) *confident* — sequences positive in both DIAMOND and HMMER searches (high specificity, low false-positive rate); and (ii) *suggestive* — sequences positive in only one search (retained for manual review). Only confident candidates were advanced to subsequent analyses.

### 2.5 Completeness filter (Phase 4)

Full-length trypsin candidates were identified among the 191 confident sequences using three mandatory criteria applied in conjunction: (i) N-terminal methionine at position 1; (ii) minimum length of 200 amino acids; and (iii) presence of a Ser195 active-site motif. Three motif classes were considered: canonical GDSGG (regex `GDS[AG]G`); conservative alternative `[GA]DS[GASTVNC]G`; and putative Ser195 variants `GD[ACTVNILMF]GG | [GA]D[ACTVNILMF][GA]G`, in which a non-serine residue occupies the Ser195 position but the flanking GDXGG context is preserved. Sequences meeting all three criteria with an intact serine motif were designated *complete* and advanced to primary characterisation; those with a putative Ser195 substitution were retained separately for evolutionary analysis.

### 2.6 Primary characterisation (Phase 5)

#### 2.6.1 Protein-level redundancy reduction — CD-HIT

Prior to characterisation, the 67 complete trypsin sequences were clustered with CD-HIT v4.8.1 (`-c 0.90 -n 5`) at 90% amino acid identity to remove near-identical isoforms generated by Trinity, retaining one representative per cluster. Sequences with atypical lengths (> 500 aa) were inspected manually for possible chimeric assemblies.

#### 2.6.2 Physicochemical properties — ProtParam

Physicochemical properties of non-redundant trypsin sequences were computed with Biopython's `ProteinAnalysis` module (Cock et al., 2009), which implements the ProtParam algorithm (Gasteiger et al., 2005). Parameters calculated included: molecular weight (kDa), theoretical isoelectric point (pI), instability index, aliphatic index, grand average of hydropathicity (GRAVY), and molar extinction coefficient.

#### 2.6.3 Signal peptide prediction — SignalP-6.0

The presence of N-terminal signal peptides — mandatory for secreted digestive enzymes — was predicted with SignalP-6.0 (Teufel et al., 2022), a transformer-based model that supports eukaryotic organisms. Sequences with predicted Sec/SPI signal peptides were considered consistent with secretory pathway localisation appropriate for midgut digestive enzymes.

#### 2.6.4 Functional annotation — eggNOG-mapper

Functional annotation of the non-redundant trypsin set was performed using eggNOG-mapper v2.1.13 (Cantalapiedra et al., 2021) against the eggNOG 5.0 database (Huerta-Cepas et al., 2019). Protein sequences were searched against the eggNOG diamond database (`eggnog_proteins.dmnd`) with DIAMOND v2.1 in sensitive mode (`--sensitive --iterate`; E-value ≤ 0.001), and the top three hits were used for orthology assignment. COG functional categories, GO terms, KEGG pathway identifiers, and ortholog descriptions were extracted from the resulting annotation table.

### 2.7 Phylogenetic analysis (Phase 6)

#### 2.7.1 Reference sequences

Reference trypsin and serine protease sequences were retrieved from NCBI Protein and UniProtKB databases via programmatic queries. Published *A. gemmatalis* sequences (including AWL83213.1 and others) were supplemented by a dynamic NCBI search (organism: *Anticarsia gemmatalis*; title: trypsin; up to 30 results). Additional Lepidoptera references were retrieved for *S. frugiperda*, *H. armigera*, *T. ni*, *M. sexta*, *B. mori*, and *L. glycinivorella* using curated accession numbers. *Bos taurus* trypsin (UniProt P00760) served as outgroup for tree rooting.

#### 2.7.2 Multiple sequence alignment and trimming

Sequences were combined (48 *A. gemmatalis* + 57 references = 105 total; 1 identical pair collapsed) and aligned using MAFFT v7 with the L-INS-i algorithm (`--localpair --maxiterate 1000`), which applies local pairwise alignment and iterative refinement for maximum accuracy on divergent sequences (Katoh & Standley, 2013). Poorly aligned and gap-rich columns were removed with trimAl v1.4.1 (`-gappyout`), retaining 255 of 1,923 alignment positions (13.3%).

#### 2.7.3 Phylogenetic inference

Maximum-likelihood trees were inferred with IQ-TREE v3.1.2 (Wong et al., 2025) using automatic model selection via ModelFinder (Kalyaanamoorthy et al., 2017). The best-fit substitution model (selected by BIC) was **Q.PFAM+G4** (Γ shape α = 1.198), calibrated for Pfam protein family data. Branch support was assessed with 1,000 ultrafast bootstrap replicates (`-B 1000`) and 1,000 SH-aLRT replicates (`-alrt 1000`). Branches with UFBoot ≥ 95% and SH-aLRT ≥ 80% were considered strongly supported.

<!-- PHASE 7 — Structure: AlphaFold3 via nf-core/proteinfold (to fill) -->
<!-- PHASE 8 — Structural validation: MolProbity + ProSA + Foldseek + ConSurf (to fill) -->
<!-- PHASE 9 — Docking: AutoDock Vina + HADDOCK + PLIP (to fill) -->
<!-- PHASE 10 — MD simulations: GROMACS 2025.1 100 ns (to fill) -->

---

## References (partial)

- Buchfink B et al. (2021) Sensitive protein alignments at tree-of-life scale using DIAMOND. *Nat Methods* 18:366–368.
- Eddy SR (2011) Accelerated profile HMM searches. *PLoS Comput Biol* 7:e1002195.
- Fu L et al. (2012) CD-HIT: accelerated for clustering the next-generation sequencing data. *Bioinformatics* 28:3150–3152.
- Grabherr MG et al. (2011) Full-length transcriptome assembly from RNA-Seq data without a reference genome. *Nat Biotechnol* 29:644–652.
- Haas BJ et al. (2013) *De novo* transcript sequence reconstruction from RNA-seq using the Trinity platform for reference generation and analysis. *Nat Protoc* 8:1494–1512.
- Levy Karin E et al. (2020) MetaEuk — sensitive, high-throughput gene discovery, and annotation for large-scale eukaryotic metagenomics. *Microbiome* 8:48.
- Li W, Godzik A (2006) Cd-hit: a fast program for clustering and comparing large sets of protein or nucleotide sequences. *Bioinformatics* 22:1658–1659.
- Manni M et al. (2021) BUSCO update: novel and streamlined workflows along with broader and deeper phylogenetic coverage for scoring of eukaryotic, prokaryotic, and viral genomes. *Mol Biol Evol* 38:4647–4654.
- Mistry J et al. (2021) Pfam: The protein families database in 2021. *Nucleic Acids Res* 49:D412–D419.
