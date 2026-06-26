#!/usr/bin/env python3
"""
04_gsea_msea_gene_overlap.py

Gene overlap analysis between transcriptomics-derived gene sets
(from GO, KEGG, and Reactome enrichment outputs) and metabolomics-
derived gene sets (from MSEA results).

This script:
1. Loads enrichment result CSV files
2. Extracts all non-empty string values from each file as gene-like entries
3. Builds unique gene sets for GO, KEGG, Reactome, and MSEA
4. Calculates shared genes across transcriptomic databases
5. Calculates shared genes between transcriptomics and metabolomics
6. Saves summary CSV outputs
7. Generates Venn diagrams

Author: Your Name
Project: Multi-omics integration
"""

from pathlib import Path
import logging
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3


# ---------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------
DATA_DIR = Path("data/integration")
OUTPUT_DIR = Path("results/integration")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

INPUT_FILES = {
    "go_up": DATA_DIR / "GO_UpRegulated.csv",
    "go_down": DATA_DIR / "GO_DownRegulated.csv",
    "kegg_up": DATA_DIR / "KEGG_UpRegulated.csv",
    "kegg_down": DATA_DIR / "KEGG_DownRegulated.csv",
    "reactome_up": DATA_DIR / "REACTOME_UpRegulated.csv",
    "reactome_down": DATA_DIR / "REACTOME_DownRegulated.csv",
    "msea": DATA_DIR / "MSEA_Test.csv",
}


# ---------------------------------------------------------------------
# Logging setup
# ---------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s"
)


# ---------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------
def check_input_files(file_dict):
    """Check whether all required input files exist."""
    missing_files = [str(path) for path in file_dict.values() if not path.exists()]
    if missing_files:
        raise FileNotFoundError(
            "The following required input files are missing:\n" + "\n".join(missing_files)
        )


def extract_all_genes_keep_numbers(df):
    """
    Extract all non-empty string values from all columns of a DataFrame.

    Note:
        This keeps the original notebook logic, where all non-null cell values
        are treated as gene-like entries. If your CSV files contain non-gene
        columns (e.g., p-values, descriptions, pathway names), consider replacing
        this with a column-specific extraction function.
    """
    gene_set = set()
    for col in df.columns:
        values = df[col].dropna().astype(str).str.strip()
        values = values[values != ""]
        gene_set.update(values.tolist())
    return gene_set


def save_venn2(set1, set2, label1, label2, colors, title, output_file):
    """Create and save a 2-set Venn diagram."""
    plt.figure(figsize=(9, 9))
    venn = venn2(
        [set1, set2],
        set_labels=(label1, label2),
        set_colors=colors,
        alpha=0.6,
    )

    if venn.set_labels:
        for text in venn.set_labels:
            if text:
                text.set_fontsize(14)

    if venn.subset_labels:
        for text in venn.subset_labels:
            if text:
                text.set_fontsize(14)

    plt.title(title, fontsize=16, fontweight="bold")
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close()


def save_venn3(set1, set2, set3, labels, colors, title, output_file):
    """Create and save a 3-set Venn diagram."""
    plt.figure(figsize=(9, 9))
    venn = venn3(
        [set1, set2, set3],
        set_labels=labels,
        set_colors=colors,
        alpha=0.6,
    )

    if venn.set_labels:
        for text in venn.set_labels:
            if text:
                text.set_fontsize(14)

    if venn.subset_labels:
        for text in venn.subset_labels:
            if text:
                text.set_fontsize(12)

    plt.title(title, fontsize=16, fontweight="bold")
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close()


# ---------------------------------------------------------------------
# Main workflow
# ---------------------------------------------------------------------
def main():
    logging.info("Starting gene overlap analysis.")

    check_input_files(INPUT_FILES)
    logging.info("All input files found.")

    # Load data
    go_up = pd.read_csv(INPUT_FILES["go_up"])
    go_down = pd.read_csv(INPUT_FILES["go_down"])
    kegg_up = pd.read_csv(INPUT_FILES["kegg_up"])
    kegg_down = pd.read_csv(INPUT_FILES["kegg_down"])
    reactome_up = pd.read_csv(INPUT_FILES["reactome_up"])
    reactome_down = pd.read_csv(INPUT_FILES["reactome_down"])
    msea_df = pd.read_csv(INPUT_FILES["msea"])

    logging.info("Input CSV files loaded successfully.")

    # Extract gene sets
    go_genes = extract_all_genes_keep_numbers(go_up) | extract_all_genes_keep_numbers(go_down)
    kegg_genes = extract_all_genes_keep_numbers(kegg_up) | extract_all_genes_keep_numbers(kegg_down)
    reactome_genes = (
        extract_all_genes_keep_numbers(reactome_up)
        | extract_all_genes_keep_numbers(reactome_down)
    )
    metabo_genes = extract_all_genes_keep_numbers(msea_df)

    all_transcriptomic_genes = go_genes | kegg_genes | reactome_genes
    shared_genes_across_databases = go_genes & kegg_genes & reactome_genes
    shared_transcriptomics_metabolomics = all_transcriptomic_genes & metabo_genes

    logging.info("Gene sets extracted.")
    logging.info("GO genes: %d", len(go_genes))
    logging.info("KEGG genes: %d", len(kegg_genes))
    logging.info("Reactome genes: %d", len(reactome_genes))
    logging.info("Metabolomics genes: %d", len(metabo_genes))
    logging.info("Shared genes across GO/KEGG/Reactome: %d", len(shared_genes_across_databases))
    logging.info(
        "Shared genes between transcriptomics and metabolomics: %d",
        len(shared_transcriptomics_metabolomics)
    )

    # Save combined gene table
    go_list = sorted(go_genes)
    kegg_list = sorted(kegg_genes)
    reactome_list = sorted(reactome_genes)
    metabo_list = sorted(metabo_genes)

    max_len = max(len(go_list), len(kegg_list), len(reactome_list), len(metabo_list))

    go_list += [None] * (max_len - len(go_list))
    kegg_list += [None] * (max_len - len(kegg_list))
    reactome_list += [None] * (max_len - len(reactome_list))
    metabo_list += [None] * (max_len - len(metabo_list))

    combined_df = pd.DataFrame({
        "GO_Genes": go_list,
        "KEGG_Genes": kegg_list,
        "Reactome_Genes": reactome_list,
        "Metabolomics_Genes": metabo_list,
    })
    combined_df.to_csv(OUTPUT_DIR / "All_Gene_Lists.csv", index=False)

    # Save shared gene outputs
    pd.DataFrame({
        "Shared_Genes_GO_KEGG_Reactome": sorted(shared_genes_across_databases)
    }).to_csv(OUTPUT_DIR / "Shared_Genes_GO_KEGG_Reactome.csv", index=False)

    pd.DataFrame({
        "Shared_Genes_Transcriptomics_Metabolomics": sorted(shared_transcriptomics_metabolomics)
    }).to_csv(
        OUTPUT_DIR / "Shared_Genes_Transcriptomics_Metabolomics.csv",
        index=False
    )

    logging.info("CSV outputs saved.")

    # Save plots
    save_venn2(
        go_genes,
        metabo_genes,
        "GO (Transcriptomics)",
        "MSEA (Metabolomics)",
        ("palegreen", "lightcoral"),
        "Venn Diagram: GO vs MSEA Genes",
        OUTPUT_DIR / "venn_GO_vs_MSEA.png",
    )

    save_venn2(
        kegg_genes,
        metabo_genes,
        "KEGG (Transcriptomics)",
        "MSEA (Metabolomics)",
        ("skyblue", "lightcoral"),
        "Venn Diagram: KEGG vs MSEA Genes",
        OUTPUT_DIR / "venn_KEGG_vs_MSEA.png",
    )

    save_venn2(
        reactome_genes,
        metabo_genes,
        "Reactome (Transcriptomics)",
        "MSEA (Metabolomics)",
        ("mediumorchid", "lightcoral"),
        "Venn Diagram: Reactome vs MSEA Genes",
        OUTPUT_DIR / "venn_Reactome_vs_MSEA.png",
    )

    save_venn3(
        go_genes,
        kegg_genes,
        reactome_genes,
        ("GO", "KEGG", "Reactome"),
        ("palegreen", "lightskyblue", "mediumorchid"),
        "Venn Diagram: GO vs KEGG vs Reactome Genes",
        OUTPUT_DIR / "venn_GO_KEGG_Reactome.png",
    )

    save_venn2(
        all_transcriptomic_genes,
        metabo_genes,
        "Transcriptomics",
        "Metabolomics",
        ("slateblue", "lightcoral"),
        "Venn Diagram: Transcriptomics vs Metabolomics Genes",
        OUTPUT_DIR / "venn_Transcriptomics_vs_Metabolomics.png",
    )

    logging.info("Venn diagrams saved.")
    logging.info("Gene overlap analysis completed successfully.")


if __name__ == "__main__":
    main()
