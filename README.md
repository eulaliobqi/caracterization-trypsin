# Trypsin-Agemmatalis-Structural

**Structural characterization of serine proteases (trypsins) from *Anticarsia gemmatalis***

[![Nextflow](https://img.shields.io/badge/Nextflow-%E2%89%A524.04-brightgreen)](https://www.nextflow.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

This Nextflow DSL2 pipeline performs end-to-end structural characterization of trypsin isoforms
from a Trinity transcriptome assembly of *A. gemmatalis* (velvetbean caterpillar) larval midgut.

**Organism:** *Anticarsia gemmatalis* Hübner 1818 (Lepidoptera: Erebidae)  
**Target:** Serine proteases — trypsins (family S1A, Pfam PF00089)  
**Goal:** Identify complete trypsin isoforms → predict 3D structures (AlphaFold3) → validate → dock with soy inhibitors (SKTI, BPTI-like peptides) → 100 ns MD simulations

## Pipeline summary

```
Trinity assembly (42,372 transcripts)
    │
    ▼
Phase 1 ─ QC Assembly (BUSCO + CD-HIT-EST)
    │
    ▼
Phase 2 ─ ORF Prediction (TransDecoder + Pfam hints)
    │
    ▼
Phase 3 ─ Trypsin Identification (DIAMOND ∩ HMMER PF00089)
    │
    ▼
Phase 4 ─ Completeness Filter (catalytic triad + Met-start + ≥220 aa)
    │
    ▼
Phase 5 ─ Primary Characterization (ProtParam + SignalP6 + InterProScan)
    │
    ▼
Phase 6 ─ Phylogeny (MAFFT-linsi + trimAl + IQ-TREE2)
    │
    ▼
Phase 7 ─ Structure Prediction (AlphaFold3 via nf-core/proteinfold) [GPU]
    │
    ▼
Phase 8 ─ Structural Validation (MolProbity + ProSA + Foldseek + ConSurf)
    │
    ▼
Phase 9 ─ Docking (AutoDock Vina + HADDOCK + PLIP)
    │
    ▼
Phase 10 ─ MD Simulations 100 ns (GROMACS 2025.1) [GPU]
    │
    ▼
Phase 11 ─ Report & Manuscript generation
```

## Quick start

```bash
# Clone repository
git clone https://github.com/eulaliobqi/trypsin-agemmatalis-structural
cd trypsin-agemmatalis-structural

# Place Trinity assembly
cp your_assembly.fasta data/raw/trinity_assembly.fasta

# Create Mamba environments (on Debian server)
for env in envs/*.yml; do mamba env create -f "$env"; done

# Run full pipeline
nextflow run nextflow/main.nf \
    -profile debian \
    --input_fasta data/raw/trinity_assembly.fasta \
    --outdir results \
    -resume

# Run with GPU (phases 7 and 10)
nextflow run nextflow/main.nf \
    -profile debian,gpu \
    --input_fasta data/raw/trinity_assembly.fasta \
    --outdir results \
    -resume
```

## Requirements

- Nextflow ≥ 24.04
- Mamba/Conda
- CUDA 12.4+ (for AlphaFold3 and GROMACS GPU)
- Debian 12+ (kernel 6.x) on execution server

## Reference databases (download separately)

| Database | Path | Used in |
|---|---|---|
| UniProt/Swiss-Prot | `data/references/uniprot_sprot.fasta` | Phase 3 |
| Pfam-A.hmm (≥36.0) | `data/references/Pfam-A.hmm` | Phases 2, 3 |
| BUSCO insecta_odb10 | `data/references/insecta_odb10/` | Phase 1 |
| Lepidoptera trypsin refs | `data/references/trypsin_refs_lepidoptera.fasta` | Phase 6 |

## Citation

If you use this pipeline, please cite:
> Gutemberg E. *et al.* (2026). Structural characterization of digestive trypsins from *Anticarsia gemmatalis*.
> GitHub: https://github.com/eulaliobqi/trypsin-agemmatalis-structural
> DOI: (Zenodo — to be assigned upon publication)

## License

MIT License — see [LICENSE](LICENSE)
