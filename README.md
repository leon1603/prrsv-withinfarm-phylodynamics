# Reconstructing Within-Farm Transmission of a PRRSV-1 Outbreak Using Structured and Unstructured Phylodynamic Models

This repository contains the computational pipeline, raw data, and analytical scripts used to quantify the transmission dynamics of Porcine Reproductive and Respiratory Syndrome Virus 1 (PRRSV-1) within a single swine herd. This research leverages Bayesian phylodynamics to distinguish between spatial levels of transmission (within-pen vs. between-pen) and evaluate the impact of genomic resolution (ORF5 vs. WGS) on epidemiological inference.

## Project Overview
This project utilizes genomic data from a longitudinal outbreak study (Clilverd et al., 2023) on a 300-sow Spanish farm to:
1.  **Reconstruct the epidemic trajectory** over a 12-month period using Birth-Death Skyline (BDSKY) models.
2.  **Quantify fine-scale spatial dynamics** using structured Birth-Death with Migration Models (BDMM).
3.  **Evaluate genomic resolution** by comparing Whole Genome Sequences (WGS) against the traditional ORF5 marker.

## Repository Structure

```text
├── data/
│   ├── 01_raw/                 # Raw FASTA sequences (WGS & ORF5)
│   ├── 02_processed/           # Processed files (Aligned, Formatted, Split)
│   ├── 03_Metadata/            # Metadata csv files
│   └── 04_InfectionClassification/ # Filtered datasets (independent transmission events)
├── scripts/
│   ├── 01_Analysis/            # Validation, correlation, and log summarization
│   ├── 02_Process/             # Data cleaning, MSA, and header reformatting
│   ├── 03_Metadata/            # Script for merging structural, infection classification and sequence similarity data
│   └── 04_Visualization/       # Plotting scripts (BDSKY (Error bar and Skyline plots), BDMM (Error bar plot), Structurally Annotated Phylograms)
├── results/
│   ├── logs/                   # Summarized BEAST2 .log files (TSV format)
│   ├── tables/                 # Sequence validation and classification tables
│   └── trees/                  # ML trees from IQ-TREE (unstructured and structured)
├── figures/                    # Generated PNG/PDF figures
├── models/                     # BEAST2 XML configurations (and raw log files if beast analysis is performed)
├── Main.qmd                    # Primary Quarto document controlling the workflow 
└── Thesis                      # Digital copy of thesis and supplementary files
```

## Workflow & Execution

The analysis is designed to be executed via `Main.qmd` in RStudio or via the Quarto CLI. The pipeline follows these sequential steps:

### 1. Data Validation & Preprocessing
* **Sequence Validation:** Automated screening for duplicates, zero-length entries, and excessive IUPAC ambiguity.
* **MSA:** Multiple Sequence Alignment using the `DECIPHER` package.
* **Date Correction:** Metadata-driven correction of sampling dates (specifically Batch 1, Week 9) to ensure accurate molecular clock calibration.

### 2. Infection Classification
To ensure the model is informed by independent transmission events, sequences are classified into:
* **Initial Infections/Reinfections:** Retained for phylodynamic analysis.
* **Persistent Infections:** Within-host evolution variants identified via hierarchical clustering (5-clade threshold), which are excluded to prevent upward bias in growth rate estimates.

### 3. Unstructured Modeling (BDSKY)
Estimates the Effective Reproductive Number ($R_e$) across the full study period.
* Models configurations: 3-, 4-, and 5-epoch strategies.
* Purpose: Recovery of the outbreak, dormancy, and resurgence phases.

### 4. Structured Modeling (BDMM)
Quantifies transmission between specific farm compartments (Farrowing vs. Nursery).
* **Spatial Levels:** Within-pen, within-room, and between-room $R_e$.
* **Migration:** Modeled as discrete movements reflecting the production cycle.
* **Ghost Deme:** Utilizes an unsampled reservoir deme to represent the sow-herd source.

## Requirements & Software

To run this analysis, ensure the following software is installed:
* **R (≥ 4.0)** with packages: `DECIPHER`, `ape`, `tidyverse`, `coda`, `ggtree`, `patchwork`.
* **Quarto:** For executing the `.qmd` pipeline.
* **IQ-TREE 3:** For Maximum Likelihood phylogenetic reconstruction.
* **BEAST v2.7.7:** Including packages `bdmm-prime`, `BDSKY`, `SampledAncestors`, and `bModelTest`.

## How to Cite
If you use this code or the findings from this study, please refer to the final report located in the `/thesis` directory:

> Balthaus, L. (2026). *Reconstructing Within-Farm Transmission of a PRRSV-1 Outbreak Using Structured and Unstructured Phylodynamic Models*. Master's Thesis, VU Amsterdam.

---
**Contact:** Leon Balthaus — [GitHub](https://github.com/leon1603)

**Citation**
Clilverd, H., Martín-Valls, G., Li, Y., Martín, M., Cortey, M., & Mateu, E. (2023). Infection dynamics, transmission, and evolution after an outbreak of porcine reproductive and respiratory syndrome virus. Frontiers in microbiology, 14, 1109881. https://doi.org/10.3389/fmicb.2023.1109881