# manuscript-code-2026
Code and analysis scripts for the study titled "Cross-Omics Integration of Metabolomic and Transcriptomic Networks reveals Key Molecular Pathways in Non-Obstructive Azoospermia"
# Cross-Omics Integration of Metabolomic and Transcriptomic Networks in NOA

Code and analysis scripts for the study: 
**"Cross-Omics Integration of Metabolomic and Transcriptomic Networks reveals Key Molecular Pathways in Non-Obstructive Azoospermia"** (2026)

## Pipeline Overview

The analysis is organized into sequential steps:

1.  **`01_preprocessing_alignment_counting.sh`**: Shell script for raw FASTQ processing, alignment (STAR/HiSAT2), and feature counting.
2.  **`02_rnaseq_de_pipeline.r`**: Differential Expression Analysis (DEA) using Limma-Voom.
3.  **`03_gsea_fgsea.R`**: Pathway enrichment analysis using the `fgsea` package.
4.  **`04_gsea_msea_gene_overlap.py`**: Python script to identify shared gene sets/pathways between Transcriptomics (GSEA) and Metabolomics (MSEA).
5.  **`05_transcriptomics_metabolomics_integration.R`**: Integration logic and visualization for multi-omics data.
6.  **`06_qpcr_relative_expression_ddct.R`**: Validation analysis of key target genes (ACAN, CYP2A13, UGT2B28, UGT2B11) using the 2^-ΔΔCt method.

## Data Structure
To run these scripts, maintain the following directory structure:
- `data/`: Raw/Processed input files (counts, metabolite lists, Ct values).
- `results/`: Output tables and publication-quality figures.
- `scripts/`: Source code (as listed above).

## Requirements
- **R**: `tidyverse`, `limma`, `fgsea`, `ggpubr`, `edgeR`.
- **Python**: `pandas`, `matplotlib_venn`, `seaborn`.

## Usage
Clone the repository and run scripts in order. For example, to run the qPCR validation:
```bash
Rscript 06_qpcr_relative_expression_ddct.R
