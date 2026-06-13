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

The 191 confident trypsin candidates were subjected to a three-criterion completeness filter: (i) N-terminal methionine (complete ORF start), (ii) minimum length of 200 amino acids, and (iii) presence of the Ser195 active-site motif GDSGG or its conservative variants (Table 4). Of the 191 sequences, 132 (69.1%) contained an initiator methionine and 184 (96.3%) met the length criterion, indicating that the majority of DIAMOND∩HMMER candidates represent full-length or near-full-length ORFs. The canonical GDSGG motif or a conservative variant was detected in 100 sequences (52.4%), with an additional seven sequences (3.7%) displaying a GDXGG pattern (X ≠ Ser) suggestive of a Ser195 substitution.

Following application of all three criteria, **67 complete trypsin sequences** were retained, with a mean length of 314 amino acids (range: 209–676 aa), consistent with the expected size range for lepidopteran trypsin zymogens including signal peptide and propeptide (250–350 aa; [ref]). An additional **6 sequences** were identified as putative Ser195 variants (GDXGG motif preserved), representing a potentially inactive or neofunctionalised trypsin-like enzyme class of evolutionary interest. Fifty-five sequences were classified as borderline (Met + length criteria met but lacking a detectable Ser195 context), and 63 were rejected due to absence of an initiator methionine or insufficient length.

The 67-sequence set is expected to include isoforms arising from alternative splicing and allelic variants, inherent to *de novo* Trinity assemblies. Protein-level sequence clustering will be applied in Phase 5 to derive a non-redundant set of functionally distinct trypsin isoforms for structural modelling.

**Table 4.** Completeness filter results for *A. gemmatalis* trypsin candidates (Phase 4).

| Category | Count | % of 191 |
|---|---|---|
| Sequences analysed | 191 | 100 |
| Met initial | 132 | 69.1 |
| ≥ 200 aa | 184 | 96.3 |
| GDSGG motif (canonical) | 100 | 52.4 |
| GDXGG motif (Ser195 variant) | 7 | 3.7 |
| **Complete (strict)** | **67** | **35.1** |
| Putative Ser195 variants | 6 | 3.1 |
| Borderline | 55 | 28.8 |
| Rejected | 63 | 33.0 |
| Mean length (complete) | 314 aa | — |
| Length range (complete) | 209–676 aa | — |

---

### 3.5 Primary characterisation of *A. gemmatalis* trypsins

#### 3.5.1 Non-redundant trypsin set

Protein-level clustering of the 67 complete trypsin sequences with CD-HIT (≥90% amino acid identity) yielded **52 non-redundant sequences**, collapsing 15 near-identical isoforms (22.4%). Automated inspection for atypically long sequences identified **four putative chimeras** (> 500 aa: DN5998, 601 aa; DN11972, 607 aa; DN8708, 676 aa; DN4838, 551 aa), which were retained in a separate set for manual BLAST verification but excluded from primary characterisation. The final working set comprised **48 trypsin sequences**.

#### 3.5.2 Physicochemical properties

Physicochemical properties of the 48 non-redundant, non-chimeric trypsin sequences were computed with Biopython's ProteinAnalysis module (Table 5). Molecular weights ranged from 22.82 to 46.19 kDa (mean 31.9 kDa), with 93.8% (45/48) of sequences falling within the expected range for Lepidoptera trypsin zymogens (20–45 kDa). The mean molecular weight of 31.9 kDa is consistent with the pre-processed form inclusive of N-terminal signal peptide and propeptide, as previously reported for *S. frugiperda* (30–33 kDa; [ref]) and *H. armigera* midgut trypsins [ref]. Three sequences exceeded 45 kDa (46.2, 45.1, and 45.6 kDa) and likely represent larger isoforms or fusion proteins warranting further investigation.

Theoretical isoelectric points ranged broadly from 4.10 to 9.58 (mean 6.8), reflecting the isoform diversity of the digestive trypsin complement. The bimodal pI distribution — with acidic isoforms (pI 4.1–5.5, n = ~15) and basic isoforms (pI 7.0–9.6, n = ~20) — is consistent with the alkaline midgut environment (pH 10–11) of *A. gemmatalis* larvae, in which alkaline-adapted trypsin isoforms predominate physiologically [Brito et al., 2001]. Protein sequences were predominantly hydrophilic (mean GRAVY = −0.019; range: −0.668 to +0.351), consistent with secreted enzymes functioning in aqueous digestive milieu. The aliphatic index ranged from 63.2 to 100.5 (mean ~84), indicating moderate thermostability for all isoforms. Based on the instability index criterion (< 40 = stable *in vitro*), 31 of 48 sequences (64.6%) were classified as stable; the 17 unstable sequences (35.4%; instability index range: 40.5–86.1) may represent isoforms with different turnover or regulatory properties.

**Table 5.** Summary of physicochemical properties of *A. gemmatalis* midgut trypsin isoforms (Phase 5; n = 48).

| Parameter | Min | Max | Mean |
|---|---|---|---|
| Molecular weight (kDa) | 22.82 | 46.19 | 31.9 |
| Isoelectric point (pI) | 4.10 | 9.58 | 6.8 |
| Sequence length (aa) | 217 | 420 | 295 |
| Instability index | 17.7 | 86.1 | — |
| GRAVY | −0.668 | +0.351 | −0.019 |
| Aliphatic index | 63.2 | 100.5 | ~84 |
| Within expected MW (20–45 kDa) | — | — | 45/48 (93.8%) |
| Stable (instability index < 40) | — | — | 31/48 (64.6%) |

#### 3.5.3 Signal peptide prediction

Signal peptide prediction with SignalP-6.0 was not available at the time of analysis (academic licence pending). This analysis will be completed prior to submission. Secretory signal peptides are expected in the majority of trypsin isoforms, consistent with their function as secreted digestive enzymes in the midgut lumen [ref].

#### 3.5.4 Functional annotation (EggNOG-mapper)

Functional annotation using eggNOG-mapper v2.1.13 against the eggNOG 5.0 database yielded annotations for all 48 sequences (100% coverage). All sequences were assigned to **COG category O** (Post-translational modification, protein turnover, chaperones — encompassing serine proteases in the eggNOG classification), COG group COG5640 (Serine protease, subtilisin/trypsin superfamily), and were classified within the Lepidoptera-level orthologous group 44642@7088 (Lepidoptera). The predominant functional annotation was "Trypsin-like serine protease", with preferred Pfam annotation **Trypsin** (PF00089).

Enzyme Commission assignments mapped to **EC 3.4.21.4** (trypsin) for the majority of sequences, with a subset assigned to EC 3.4.21.1 (chymotrypsin-like). KEGG Orthology assignments included **ko:K01310** (trypsin-1) and **ko:K01312** (trypsin-2), and KEGG Pathway annotation placed sequences in the "Protein digestion and absorption" pathway (ko04974) — directly consistent with the biological role of these enzymes in midgut-mediated soybean protein hydrolysis. Key Gene Ontology terms assigned included GO:0004252 (serine-type endopeptidase activity), GO:0006508 (proteolysis), and GO:0017171 (serine hydrolase activity). Seed orthologs were predominantly from *Bombyx mori* (Lepidoptera: Bombycidae) and *Leguminivora glycinivorella* (Lepidoptera: Tortricidae), indicating strong conservation of trypsin structure within Lepidoptera.

---

### 3.6 Phylogenetic analysis

The 48 non-redundant *A. gemmatalis* trypsin sequences were subjected to phylogenetic analysis alongside 57 reference sequences comprising published trypsin and serine protease sequences from *A. gemmatalis* (including AWL83213.1 and additional sequences retrieved by NCBI programmatic search), other Noctuidae/Erebidae (*Heliothis virescens*, *Agrotis ipsilon*, *Mamestra configurata*), and representatives of broader Lepidoptera clades (*Spodoptera frugiperda*, *Helicoverpa armigera*, *Trichoplusia ni*, *Manduca sexta*, *Bombyx mori*, *Leguminivora glycinivorella*), plus bovine trypsin (P00760) as outgroup. After MAFFT-linsi alignment (1,923 positions) and trimAl (-gappyout; 255 positions retained, 13.3%), the dataset comprised 104 unique sequences.

Maximum-likelihood phylogenetic inference with IQ-TREE v3.1.2 (Wong et al., 2025) selected **Q.PFAM+G4** as the best-fit model by BIC (BIC score = 74,347.748; Γ shape α = 1.198), substantially outperforming LG+G4 (ΔBIC = 88.4) and all other tested models. Q.PFAM is calibrated on Pfam protein family alignments, making it particularly appropriate for trypsin/serine protease phylogenomics. Tree search converged after 561 iterations (final log-likelihood = −36,563.837; s.e. = 765.3; ultrafast bootstrap correlation coefficient = 0.996), with 1,000 ultrafast bootstrap replicates (UFBoot) and 1,000 SH-aLRT replicates used to assess branch support. All 255 alignment sites were parsimony-informative (254/255; 99.6%), indicating high phylogenetic signal despite the conservative trimming. Total tree length was 87.42, with internal branches accounting for 23.31 (26.7%), consistent with deep divergence among the included taxa.

[**Figure X.** Maximum-likelihood phylogenetic tree — to be generated in Phase 11 with annotated clades and bootstrap values]

Preliminary inspection of the topology reveals that *A. gemmatalis* trypsin sequences group into multiple distinct clades, consistent with the diversification of the digestive trypsin complement reported in other Lepidoptera [ref]. Seven sequences with > 50% alignment gaps were retained for the initial analysis but will be excluded in the final publication-quality tree.

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
