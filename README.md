manuscript-code-2026

Code and analysis scripts for the study:

Cross-Omics Integration of Metabolomic and Transcriptomic Networks Reveals Key Molecular Pathways in Non-Obstructive Azoospermia (2026)

Overview

This repository contains the analysis pipeline used for transcriptomic, metabolomic, and qPCR validation in non-obstructive azoospermia (NOA). The workflow includes RNA-seq preprocessing, differential expression analysis, pathway enrichment, cross-omics overlap analysis, multi-omics integration, and qPCR validation.

Repository Structure

manuscript-code-2026/
├── data/       # input data files
├── results/    # output tables and figures
├── scripts/    # analysis scripts
└── README.md


Pipeline

01_preprocessing_alignment_counting.sh
Raw FASTQ quality control, alignment, and gene-level counting.

02_rnaseq_de_pipeline.R
Differential expression analysis using edgeR and limma-voom.

03_gsea_fgsea.R
Gene set enrichment analysis using fgsea.

04_venn_gene_overlap_analysis.py
Identification of shared genes/pathways across transcriptomic and metabolomic enrichment results, with Venn diagram visualization.

05_transcriptomics-metqabolomics-integraqtion.R
Integration of transcriptomic and metabolomic findings to prioritize key genes and pathways.

06_qpcr_ddct_analysis.R
qPCR validation of selected target genes using the 2^-ΔΔCt method.

Input Data

Place required input files in the data/ directory. These may include:
RNA-seq count matrices
sample metadata
GSEA/MSEA pathway results
metabolomics-derived candidate lists
qPCR Ct tables

Output
Results will be written to the results/ directory, including:
differential expression tables
enrichment results
overlap summary tables
Venn diagrams
integrated candidate gene tables
qPCR summary statistics and plots

Requirements

R packages:
tidyverse
edgeR
limma
fgsea
ggplot2
ggpubr

Python packages:
pandas
matplotlib
matplotlib-venn
seaborn


Scripts are organized to match the manuscript workflow.

Citation

If you use this code, please cite the associated manuscript.
