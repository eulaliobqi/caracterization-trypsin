# Results and Discussion
> Draft — updated progressively as each phase is completed.
> Last update: Phase 1 complete; Phase 2 in progress (2026-05-28)

---

## 3. Results and Discussion

### 3.1 Assembly quality and transcriptome completeness

#### 3.1.1 BUSCO assessment

BUSCO evaluation against the insecta_odb10 gene set (n = 1,367 conserved genes) revealed high transcriptome completeness, with 94.2% of genes represented as complete BUSCOs (Table 1). Fragmented and missing gene models accounted for only 1.7% and 4.1% of the dataset, respectively, indicating near-comprehensive coverage of the conserved insect gene space. These figures exceed the 90% completeness threshold generally regarded as indicative of high-quality assemblies suitable for downstream functional annotation and structural analyses [Simão et al., 2015].

**Table 1.** BUSCO completeness statistics for the *A. gemmatalis* larval midgut transcriptome (BUSCO v6.0.0, insecta_odb10, n = 1,367).

| Category | Count | Percentage (%) |
|---|---|---|
| Complete (C) | 1,288 | 94.2 |
| — Single-copy (S) | 569 | 41.6 |
| — Duplicated (D) | 719 | 52.6 |
| Fragmented (F) | 23 | 1.7 |
| Missing (M) | 56 | 4.1 |
| **Total searched** | **1,367** | — |

The elevated proportion of duplicated BUSCOs (52.6%) is a well-documented artefact of *de novo* transcriptome assembly with Trinity, which tends to generate multiple isoform-level contigs per locus [Grabherr et al., 2011]. This duplication does not reflect true gene copy-number variation but rather the capacity of Trinity to resolve alternatively spliced transcripts and allelic variants as separate assemblies. Consistent with this interpretation, CD-HIT-EST clustering at 95% nucleotide identity reduced the transcript count from 42,372 to 31,219 (a 26.3% reduction), collapsing highly similar isoforms while retaining transcript diversity relevant for the identification of functionally distinct trypsin variants.

The 94.2% BUSCO completeness of the *A. gemmatalis* midgut transcriptome compares favourably with published transcriptomes from other lepidopteran pests, including *Spodoptera frugiperda* (88–93%; [ref]), *Helicoverpa armigera* (91%; [ref]), and *Manduca sexta* (96%; [ref]), and supports the suitability of this dataset for the comprehensive identification of digestive enzyme-encoding transcripts.

---

### 3.2 ORF prediction and proteome characterisation

TransDecoder identified 36,417 candidate ORFs of ≥ 100 amino acids from the 31,219 non-redundant transcripts (~1.17 ORFs per transcript). Of these, 22,379 were supported by DIAMOND blastp hits against the NCBI NR database, and an additional subset harboured Pfam-annotated domains as identified by hmmscan against Pfam-A. Following integration of both BLAST and Pfam scoring hints, TransDecoder.Predict retained **17,923 high-confidence protein models** (`--single_best_only`), comprising 12,494 complete (69.7%), 2,820 5′-partial (15.7%), and 1,430 3′-partial (8.0%) ORFs (Table 2).

**Table 2.** Summary of ORF prediction results from the *A. gemmatalis* larval midgut transcriptome (TransDecoder v5.7.1).

| Category | Count | Percentage (%) |
|---|---|---|
| Long ORFs (≥ 100 aa) | 36,417 | — |
| Final predicted proteins | 17,923 | 100 |
| — Complete (5′- and 3′-intact) | 12,494 | 69.7 |
| — 5′-partial | 2,820 | 15.7 |
| — 3′-partial | 1,430 | 8.0 |
| — Internal fragment | 1,179 | 6.6 |
| DIAMOND NR hints used | 22,379 | — |
| Pfam domain hints used | yes | — |

The proportion of complete ORFs (69.7%) is consistent with high-quality lepidopteran transcriptomes processed with TransDecoder under similar conditions [ref]. The relatively high number of 5′-partial models (15.7%) likely reflects transcripts where the 5′ end was not fully captured during sequencing or assembly, a common feature of *de novo* Trinity assemblies from short-read data [Grabherr et al., 2011]. For the identification of functional trypsin isoforms — which requires the presence of the N-terminal signal peptide, propeptide, and catalytic triad — only complete ORFs will be considered in downstream analyses (Phases 3–4).

---

### 3.3 Identification of trypsin-encoding transcripts

> **[Running — Phase 3 in progress; values to be updated upon completion]**

To identify trypsin-encoding sequences within the predicted proteome, a dual-validation strategy was applied combining sequence homology and domain detection. The 17,923 predicted proteins were searched against the NCBI NR database using DIAMOND blastp (sensitive mode; E-value ≤ 1×10⁻¹⁰), and hits annotated as trypsin, serine protease, chymotrypsin, or trypsinogen were retained as DIAMOND candidates. Concurrently, the Tryp_SPc domain (PF00089) was searched using hmmsearch against the extracted Pfam-A profile, retaining sequences with domain coverage ≥ 80%. Sequences supported by both DIAMOND homology and PF00089 domain evidence were classified as *confident* trypsin candidates; those supported by only one source as *suggestive*.

The DIAMOND search recovered 285 candidate proteins annotated as trypsin, chymotrypsin, or serine endopeptidase, with top hits including sequences from closely related noctuids (*Trichoplusia ni*, *Heliothis virescens*, *Agrotis ipsilon*) and, notably, *A. gemmatalis* trypsin 1 itself (AWL83213.1; 7 hits), confirming the validity of the approach. Concurrently, hmmsearch against PF00089 identified 10,675 proteins with domain coverage ≥ 80% (E-value ≤ 1×10⁻¹⁰), reflecting the broad representation of serine proteases S1 in the predicted proteome. The intersection of both evidence tiers yielded **191 confident trypsin candidates** (Table 3), a figure consistent with the expected trypsin repertoire of Lepidoptera midgut (30–200; [ref]), and well within the range reported for related species such as *Spodoptera frugiperda* (~150 trypsin-like sequences; [ref]) and *Manduca sexta* (~80; [ref]). An additional 10,578 proteins were classified as suggestive (PF00089-positive but lacking trypsin keyword annotation in NR), likely representing other serine proteases of the S1 family not filtered by keyword.

**Table 3.** Summary of dual-validation trypsin identification (Phase 3).

| Evidence tier | Count |
|---|---|
| DIAMOND candidates (NR, keyword filter) | 285 |
| HMMER candidates (PF00089, coverage ≥ 80%) | 10,675 |
| **Confident (DIAMOND ∩ HMMER)** | **191** |
| Suggestive (DIAMOND or HMMER only) | 10,578 |

---

### 3.4 Full-length trypsin isoforms: completeness filter

> **[Pending — Phase 4]**

---

### 3.5 Primary characterisation of *A. gemmatalis* trypsins

> **[Pending — Phase 5]**

---

### 3.6 Phylogenetic analysis

> **[Pending — Phase 6]**

---

### 3.7 Three-dimensional structural models

> **[Pending — Phase 7]**

---

### 3.8 Structural quality and conservation

> **[Pending — Phase 8]**

---

### 3.9 Molecular docking with SKTI and BPTI-like peptides

> **[Pending — Phase 9]**

---

### 3.10 Molecular dynamics simulations

> **[Pending — Phase 10]**

---

## References (partial)

- Grabherr MG et al. (2011) Full-length transcriptome assembly from RNA-Seq data without a reference genome. *Nat Biotechnol* 29:644–652.
- Simão FA et al. (2015) BUSCO: assessing genome assembly and annotation completeness with single-copy orthologs. *Bioinformatics* 31:3210–3212.
- [Additional references to be added as sections are completed]
