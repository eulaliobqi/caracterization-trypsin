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

<!-- PHASE 3 — Trypsin identification: DIAMOND + HMMER (to fill) -->
<!-- PHASE 4 — Completeness filter: catalytic triad + Met + length (to fill) -->
<!-- PHASE 5 — Primary characterisation: ProtParam + SignalP6 + InterProScan (to fill) -->
<!-- PHASE 6 — Phylogeny: MAFFT + trimAl + IQ-TREE2 (to fill) -->
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
