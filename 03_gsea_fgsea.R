# 03_gsea_fgsea.R
# Gene Set Enrichment Analysis (GSEA) using fgsea
# Databases: GO, REACTOME, KEGG

rm(list = ls(all.names = TRUE))
gc()

options(
  max.print = .Machine$integer.max,
  scipen = 999,
  stringsAsFactors = FALSE,
  dplyr.summarise.inform = FALSE
)

suppressPackageStartupMessages({
  library(tidyverse)
  library(fgsea)
  library(data.table)
  library(stringr)
})

# -------------------------------------------------------------------
# Paths
# -------------------------------------------------------------------
infile <- "data/pval_fc_table.csv"
gmt_dir <- "data/Background_genes"
out_dir <- "results"

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# -------------------------------------------------------------------
# Helper functions
# -------------------------------------------------------------------
matrix_to_list <- function(pws) {
  pws_list <- list()
  for (pw in colnames(pws)) {
    pws_list[[pw]] <- rownames(pws)[as.logical(pws[, pw])]
  }
  return(pws_list)
}

prepare_gmt <- function(gmt_file, genes_in_data) {
  gmt <- gmtPathways(gmt_file)
  all_genes <- unique(unlist(gmt))

  mat <- matrix(
    NA,
    dimnames = list(all_genes, names(gmt)),
    nrow = length(all_genes),
    ncol = length(gmt)
  )

  for (i in seq_len(ncol(mat))) {
    mat[, i] <- as.numeric(all_genes %in% gmt[[i]])
  }

  genes_overlap <- intersect(genes_in_data, all_genes)
  mat <- mat[
    genes_overlap,
    colnames(mat)[colSums(mat[genes_overlap, , drop = FALSE]) > 5],
    drop = FALSE
  ]

  matrix_to_list(mat)
}

make_rankings <- function(df) {
  rankings <- sign(df$logFC) * (-log10(df$P.Value))
  names(rankings) <- df$gene_symbol

  rankings <- rankings[!is.na(names(rankings))]
  rankings <- rankings[!is.na(rankings)]
  rankings <- rankings[names(rankings) != "None"]
  rankings <- rankings[names(rankings) != ""]

  rankings <- sort(rankings, decreasing = TRUE)
  rankings <- rankings[!duplicated(names(rankings))]

  set.seed(123)
  rankings <- rankings + rnorm(length(rankings), mean = 0, sd = 1e-5)
  rankings <- sort(rankings, decreasing = TRUE)

  return(rankings)
}

clean_pathway_names <- function(pathway_vec, db_name) {
  if (db_name == "GO") {
    pathway_vec <- gsub("^GOBP_", "", pathway_vec)
  } else if (db_name == "KEGG") {
    pathway_vec <- gsub("^KEGG_MEDICUS_REFERENCE_", "", pathway_vec)
    pathway_vec <- gsub("^KEGG_", "", pathway_vec)
  } else if (db_name == "REACTOME") {
    pathway_vec <- gsub("^REACTOME_", "", pathway_vec)
  }

  pathway_vec <- gsub("_", " ", pathway_vec)
  return(pathway_vec)
}

run_gsea_for_db <- function(db_name, gmt_file, df, out_dir) {
  message("Running GSEA for: ", db_name)

  genes_in_data <- df$gene_symbol
  pathways <- prepare_gmt(gmt_file, genes_in_data)
  rankings <- make_rankings(df)

  pathway_gene_counts <- sapply(pathways, function(genes) {
    sum(genes %in% names(rankings))
  })

  pathways <- pathways[pathway_gene_counts >= 10]

  gsea_res <- fgsea(
    pathways = pathways,
    stats = rankings,
    scoreType = "std",
    minSize = 10,
    maxSize = 500,
    nproc = 1,
    nPermSimple = 100000
  )

  gsea_res <- as.data.frame(gsea_res) %>%
    arrange(padj, pval)

  # save table
  fwrite(
    gsea_res,
    file = file.path(out_dir, paste0("NOAvsControl_", db_name, "_gsea_results.tsv")),
    sep = "\t",
    sep2 = c("", " ", "")
  )

  write.csv(
    gsea_res,
    file = file.path(out_dir, paste0("NOAvsControl_", db_name, "_gsea_results.csv")),
    row.names = FALSE
  )

  # top up/down barplot
  top_up <- gsea_res %>%
    filter(ES > 0) %>%
    arrange(padj, pval) %>%
    head(10)

  top_down <- gsea_res %>%
    filter(ES < 0) %>%
    arrange(padj, pval) %>%
    head(10)

  top_pathways <- bind_rows(top_up, top_down) %>%
    mutate(pathway_clean = clean_pathway_names(pathway, db_name))

  p_bar <- ggplot(top_pathways,
                  aes(x = reorder(pathway_clean, NES), y = NES, fill = NES)) +
    geom_col() +
    coord_flip() +
    theme_bw() +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
    scale_x_discrete(labels = function(x) str_trunc(x, width = 55)) +
    labs(
      title = paste0(db_name, ": Top 10 Upregulated & Downregulated Pathways"),
      x = "Pathway",
      y = "Normalized Enrichment Score (NES)"
    ) +
    theme(
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 10)
    )

  ggsave(
    filename = file.path(out_dir, paste0("barplot_up_down_", db_name, ".png")),
    plot = p_bar,
    width = 15,
    height = 8,
    dpi = 600
  )

  # top pathway enrichment plot
  top_pathway <- gsea_res$pathway[which.min(gsea_res$padj)]

  p_enrich <- plotEnrichment(
    pathways[[top_pathway]],
    rankings
  ) + labs(title = top_pathway)

  ggsave(
    filename = file.path(out_dir, paste0("enrichment_top_", db_name, ".png")),
    plot = p_enrich,
    width = 10,
    height = 6,
    dpi = 600
  )

  return(gsea_res)
}

# -------------------------------------------------------------------
# Load differential expression table
# -------------------------------------------------------------------
dat <- read.csv(infile, row.names = 1)

df <- dat %>%
  dplyr::select(gene_symbol, P.Value, adj.P.Val, logFC)

# -------------------------------------------------------------------
# Database mapping
# -------------------------------------------------------------------
gmt_files <- c(
  GO = file.path(gmt_dir, "GO.gmt"),
  REACTOME = file.path(gmt_dir, "REACTOME.gmt"),
  KEGG = file.path(gmt_dir, "KEGG.gmt")
)

# -------------------------------------------------------------------
# Run analysis
# -------------------------------------------------------------------
all_results <- lapply(names(gmt_files), function(db) {
  run_gsea_for_db(
    db_name = db,
    gmt_file = gmt_files[[db]],
    df = df,
    out_dir = out_dir
  )
})

names(all_results) <- names(gmt_files)
