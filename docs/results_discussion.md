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

> **[In progress — Phase 2 running]**
> *To be updated upon completion of TransDecoder.Predict + validation.*

TransDecoder identified [X] candidate ORFs of ≥ 100 amino acids from the 31,219 non-redundant transcripts, of which [Y] were supported by DIAMOND homology hits against the NCBI NR database and [Z] harboured Pfam-annotated domains. Following integration of these homology and domain-based scoring hints, TransDecoder.Predict retained [N] high-confidence protein models (`--single_best_only`), comprising [n_complete] complete (5′- and 3′-intact), [n_5] 5′-partial, and [n_3] 3′-partial ORFs.

---

### 3.3 Identification of trypsin-encoding transcripts

> **[Pending — Phase 3]**

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
