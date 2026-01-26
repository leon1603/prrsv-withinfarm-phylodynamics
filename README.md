# Reconstructing Within-Farm Transmission of a PRRSV-1 Outbreak Using Structured and Unstructured Phylodynamic Models

This repository contains the computational pipeline, raw data, and analytical scripts used to quantify the transmission dynamics of Porcine Reproductive and Respiratory Syndrome Virus 1 (PRRSV-1) within a single swine herd. This research leverages Bayesian phylodynamics to distinguish between spatial levels of transmission (within-pen vs. between-pen) and evaluate the impact of genomic resolution (ORF5 vs. WGS) on epidemiological inference.

## Project Overview

This project utilizes genomic data from a longitudinal outbreak study (Clilverd et al., 2023) on a 300-sow Spanish farm to: 1. **Reconstruct the epidemic trajectory** over a 12-month period using Birth-Death Skyline (BDSKY) models. 2. **Quantify fine-scale spatial dynamics** using structured Birth-Death with Migration Models (BDMM). 3. **Evaluate genomic resolution** by comparing Whole Genome Sequences (WGS) against the traditional ORF5 marker.

## Repository Structure

``` text
├── data/
│   ├── 01_raw/                 # Raw FASTA sequences (WGS & ORF5)
│   ├── 02_processed/           # Processed files (Aligned, Formatted, Split)
│   ├── 03_Metadata/            # Metadata csv files
│   └── 04_InfectionClassification/ # Filtered datasets (independent transmission events)
├── scripts/
│   ├── 01_Analysis/            # Validation, correlation, and log summarization
│   ├── 02_Process/             # Data cleaning, MSA, and header reformatting
│   ├── 03_Metadata/            # Script for merging structural, infection classification and sequence similarity data
│   └── 04_Visualization/       # Plotting scripts (BDSKY (Error bar and Skyline), BDMM (Error bar), Structurally Annotated Phylograms)
├── results/
│   ├── logs/                   # Summarized BEAST2 .log files (TSV format)
│   ├── tables/                 # Sequence validation and classification tables
│   └── trees/                  # ML trees from IQ-TREE (unstructured and structured)
├── figures/                    # Generated PNG/PDF figures
├── models/                     # BEAST2 XML configurations
├── Main.qmd                    # Primary Quarto document controlling the workflow 
└── Thesis                      # Digital copy of thesis and supplementary files
```

## Workflow & Execution

The analysis is designed to be executed via `Main.qmd` in RStudio. The pipeline contains these components:

### 1. Data Validation & Preprocessing

-   **Sequence Validation:** Automated screening for duplicates, zero-length entries, and excessive IUPAC ambiguity.
-   **MSA:** Multiple Sequence Alignment.
-   **Date Correction:** Metadata-driven correction of sampling dates (specifically Batch 1, Week 9).

### 2. Maximum Likelihood Phylogenetic Reconstruction (IQ-TREE)

IQ-TREE is utilized at two distinct stages of the pipeline:

-   **First Pass (Unstructured Data):** Generated from alignments to provide the foundation for:

    -   **Clade-based Infection Classification:** Identifying persistent vs. independent infections.

    -   **Lateral Introduction Analysis:** Screening for external virus entry.

    -   **Clock Signal Validation:** Performing Root-to-Tip regression to assess temporal signal.

-   **Second Pass (Structured Data):** Generated using sequence headers annotated with structural metadata. These trees facilitate the qualitative assessment of spatial clustering via structurally annotated Phylograms.

### 3. Infection Classification

To ensure the model is informed by independent transmission events, sequences are classified into: Initial Infections/Reinfections: Retained for phylodynamic analysis and Persistent Infections: Within-host evolution variants identified via hierarchical clustering (5-clade threshold), which are excluded to prevent upward bias in growth rate estimates.

### 4. Bayesian Phylodynamic Modeling (BEAST2)

The core analysis utilizes both Birth-Death Skyline (BDSKY) (Unstructured) and Birth Death with Migration Model (BDMM) (Structured) models.

> [!I MPORTANT] **Data Availability Note:** Due to their significant size, the raw `.log` and `.trees` output files from BEAST2 are not included in this repository.
>
> -   **To Reproduce Results:** Execute the `.xml` files located in the `models/` directory using BEAST2.
>
> -   **To Summarize Data:** Feed the resulting raw log files into the `summarize_beast_logs` function within `scripts/01_Analysis/04_summariseLog.R` (as automated in `Main.qmd`).
>
> -   **Pre-calculated Results:** The summarized TSV versions of these logs are available in `results/logs/`for immediate use in the visualization scripts.

## Requirements & Software

To run this analysis, ensure the following core software is installed:

-   **R (≥ 4.0)**

-   **IQ-TREE 3 (Wong et al., 2025):** For Maximum Likelihood phylogenetic reconstruction.

-   **BEAST v2.7.7 (Bouckaert et al., 2019):** Ensure the following BEAST-specific packages are installed via the BEAST2 Package Manager:`bdmm-prime`, `BDSKY`, `SampledAncestors`, and `bModelTest`.

#### Managing R Dependencies with `renv`

This project uses the `renv` package to ensure a reproducible environment.  You do not need to install R packages manually.To set up the environment:

1.  Open the project in RStudio (or set your working directory to the project root).

2.  Install `renv` if you haven't already: `install.packages("renv")`.

3.  Run the following command to restore the local library:

``` r
renv::restore()
```

This will automatically install the correct versions of all required R packages (including `DECIPHER`, `ape`, `tidyverse`, `coda`,`ggtree`, `patchwork`, etc.) as specified in the `renv.lock` file.

## AI Acknowledgement

Google Gemini was used as a computational resource during this project. Specifically, Gemini assisted in generating the basic layout for the computational code and served as a technical resource for debugging throughout the data analysis process.

## **Citations**

Bouckaert, R., Vaughan, T. G., Barido-Sottani, J., Duchêne, S., Fourment, M., Gavryushkina, A., Heled, J., Jones, G., Kühnert, D., De Maio, N., Matschiner, M., Mendes, F. K., Müller, N. F., Ogilvie, H. A., du Plessis, L., Popinga, A., Rambaut, A., Rasmussen, D., Siveroni, I., Suchard, M. A., … Drummond, A. J. (2019). BEAST 2.5: An advanced software platform for Bayesian evolutionary analysis.  *PLoS computational biology*, *15*(4), e1006650. <https://doi.org/10.1371/journal.pcbi.1006650>

Clilverd, H., Martín-Valls, G., Li, Y., Martín, M., Cortey, M., & Mateu, E. (2023). Infection dynamics, transmission, and evolution after an outbreak of porcine reproductive and respiratory syndrome virus. Frontiers in microbiology, 14, 1109881. <https://doi.org/10.3389/fmicb.2023.1109881>

Wong, T. K., Ly-Trong, N., Ren, H., Baños, H., Roger, A. J., Susko, E., ... & Minh, B. Q. (2025). IQ-TREE 3: Phylogenomic Inference Software using Complex Evolutionary Models. <https://doi.org/10.32942/X2P62N>
