# prrsv-withinfarm-phylodynamics

## Project Overview
This repository contains the data, code, and phylodynamic model configurations for the analysis of PRRSV-1 within-farm transmission dynamics. 

## Methodology
The analysis is divided into two modeling approaches as described in the thesis:
1. **Long-term Dynamics (BDSKYSA):** Unstructured analysis of the full 12-month period to reconstruct outbreak phases.
2. **Fine-Scale Dynamics (SMTBD):** Structured Multi-type Birth-Death analysis of the initial outbreak (Batch 1) to resolve spatial transmission between farrowing and nursery units.

## Data Source
Data derived from the longitudinal study by Clilverd et al. (2023).

## Directory Structure
* `models/`: Contains BEAST 2 XML configurations.
* `scripts/`: R scripts for pre-processing and post-processing.
* `data/`: Processed sequence data and metadata.

