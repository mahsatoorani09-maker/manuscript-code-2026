counts_file <- "data/count_matrix.csv"
meta_file   <- "data/sample_metadata.csv"

# Load required libraries
library(edgeR)
library(org.Hs.eg.db)
library(dplyr)

# Seq data loading (already done in your original code)
seqdata <- read.csv("merge.csv")
colnames(seqdata) <- substr(colnames(seqdata), start = 1, stop = 11)
countdata <- seqdata[,-(1)]
rownames(countdata) <- seqdata[,1]
head(countdata)
countdata <- countdata[, -(12)]
countdata
# Sample information (already done in your original code)
sampleinfo <- read.csv("Metadata.csv")

sampleinfo <- sampleinfo %>%
  mutate(Group = recode(Group,
                        "SP_NOA" = "Azoospermia",
                        "rST_NOA" = "Azoospermia",
                        "OA" = "Control"))
sampleinfo
sampleinfo <- sampleinfo [-(12), ]
sampleinfo

write_csv(sampleinfo, "updated_Merge.csv")

# Convert counts to DGEList object
y <- DGEList(countdata)

group <- paste(sampleinfo$Group)
group
group<- factor(group)

y$samples$group <- group
head(y)
# Add annotation usgroup# Add annotation using org.Hs.eg.db package
columns(org.Hs.eg.db)  # Check available columns in org.Hs.eg.db

# Perform annotation


library(AnnotationDbi)
library(org.Hs.eg.db)

# Perform annotation
ann <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = rownames(y$counts),     # These are Ensembl IDs
  keytype = "ENSEMBL",           # Keytype is "ENSEMBL"
  columns = c("ENTREZID", "SYMBOL", "GENENAME")  # Add desired annotation columns
)





# Match annotation to counts data
ann_matched <- ann[match(rownames(y$counts), ann$ENSEMBL), ]
y$genes <- ann_matched
head(y$genes)
# Inspect the annotated data
head(y$genes)

# Filter out low-abundance genes (optional, but recommended)
myCPM <- cpm(countdata)
thresh <- myCPM > 1

keep <- rowSums(thresh) >= 2
countdata_filtered <- countdata[keep, ]



# Convert filtered counts to DGEList object
y <- DGEList(countdata_filtered)

y$samples$group <- factor(sampleinfo$Group)

# Normalize using TMM (Trimmed Mean of M-values)
y <- calcNormFactors(y, method = "TMM")
y
# Visualize library sizes before and after normalization

colors <- ifelse(y$samples$group == "Control", "purple", "yellow")
# Control = purple, Azoo = yellow
barplot(y$samples$lib.size / 1e06, names = colnames(y), las = 2, col = colors,
        main = "Library Sizes", ylab = "Library Size (Millions)")

colors <- ifelse(y$samples$group == "Control", "purple", "yellow")


barplot(y$samples$lib.size * y$samples$norm.factors / 1e06, 
        names = colnames(y), 
        las = 2, 
        col = colors,  # رنگ‌ها براساس گروه
        main = "Library Sizes (normalized)", 
        ylab = "Normalized Library Size (Millions)")




# Get log2 counts per million
logcounts <- cpm(y,log=TRUE)
# Check distributions of samples using boxplots
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
title("logCPMs (unnormalised)")
par(mfrow = c(1, 1)) 
par(mar = c(5, 5, 4, 2))  


sampleinfo$Group <- as.factor(sampleinfo$Group)
col.cell <- c("purple", "orange")[sampleinfo$Group]


head(y$counts)


logcounts <- cpm(y, log=TRUE, prior.count=1) 


summary(as.vector(logcounts))

keep <- filterByExpr(y)
y_filtered <- y[keep, , keep.lib.sizes=FALSE]
logcounts_filtered <- cpm(y_filtered, log=TRUE)
summary(as.vector(logcounts_filtered))
#get log 2 counts per million
logcounts <- cpm (y, log=TRUE)

head(logcounts)
write.csv(logcounts, "logCMP.csv")
# Create the design matrix to model differences between experimental groups
design <- model.matrix(~ y$samples$group)

# Estimate dispersion for each gene to account for biological variability
y <- estimateDisp(y, design)

# Fit a Quasi-Likelihood Generalized Linear Model (GLM) to the data
fit <- glmQLFit(y, design)

# Conduct the Quasi-Likelihood F-test to identify differentially expressed genes
qlf <- glmQLFTest(fit)

# Extract all genes with their respective statistics (logFC, p-value, FDR)
res <- topTags(qlf, n = Inf)

# Print results summary
print(res)
summary(res)

# Filter results based on the significance threshold (FDR < 0.01)
res_filtered <- res$table[res$table$FDR < 0.01, ]
print(res_filtered)

# ==============================================================================
# Data Saving and Exploratory Data Analysis (PCA)
# ==============================================================================

# Save normalized and annotated objects for future downstream analysis
save(y, logcounts, sampleinfo, ann_matched, file = "normalized_annotated_data.RData")

# Perform Principal Component Analysis (PCA) on log-transformed data
pca_result <- prcomp(t(logcounts), scale. = TRUE)

# Plot the PCA results (PC1 vs PC2)
plot(pca_result$x[, 1:2], 
     col = as.numeric(y$samples$group), 
     pch = 19, 
     xlab = "PC1", 
     ylab = "PC2", 
     main = "PCA Plot")

# Add sample labels to the points on the plot
text(pca_result$x[, 1], pca_result$x[, 2], 
     labels = rownames(y$samples), 
     pos = 3, cex = 0.8, col = "black")


# Multidimensional scaling (MDS) plot
sampleinfo$Group <- as.factor(sampleinfo$Group)

col.cell <- c("purple", "orange")[sampleinfo$Group]

plotMDS(y,
        col = col.cell,
        cex = 1.2,             
        cex.axis = 1.2,        
        cex.lab = 1.4,         
        xlab = "Leading logFC dim 1 (49%)", 
)



title("condition", cex.main = 1.5) 
legend("topleft", fill = c("purple", "orange"),
       legend = levels(sampleinfo$Group),
       cex = 1.2, bty = "n")

library(RColorBrewer)
library(gplots)
var_genes <- apply(logcounts, 1, var)
head(var_genes)
# Get the gene names for the top 500 most variable genes
select_var <- names(sort(var_genes, decreasing=TRUE))[1:500]
head(select_var)
# Subset logcounts matrix
highly_variable_lcpm <- logcounts[select_var,]
dim(highly_variable_lcpm)
head(highly_variable_lcpm)




## Get some nicer colours
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
# Set up colour vector for celltype variable
col.cell <- c("purple","orange")[sampleinfo$Group]


  # Create sparse gene labels (show only some genes)
  gene_labels <- rownames(highly_variable_lcpm)
  
  sparse_labels <- rep("", length(gene_labels))
  
  # show one gene name every 15 genes
  idx <- seq(1, length(gene_labels), by = 1)
  
  sparse_labels[idx] <- gene_labels[idx]
  
  # Save heatmap
  png("C:/Users/ASUS/Desktop/heat_plot_clean24.png",
      width = 7000,
      height = 5000,
      res = 400)
  
  heatmap.2(
    highly_variable_lcpm,
    
    col = rev(morecols(50)),
    trace = "none",
    
    main = "Top 500 most variable genes across samples",
    
    ColSideColors = col.cell,
    
    scale = "row",
    
    margins = c(18, 12),
    
    cexCol = 1.3,
    cexRow = 1.3,
    
    srtCol = 90,
    
    labCol = colnames(highly_variable_lcpm),
    
    labRow = sparse_labels
  )
  
  

dev.off()
#ُSVG
# Save heatmap as SVG
svg("C:/Users/ASUS/Desktop/heat_plot_clean24.svg",
    width = 14,
    height = 10)

heatmap.2(
  highly_variable_lcpm,
  
  col = rev(morecols(50)),
  trace = "none",
  
  main = "Top 500 most variable genes across samples",
  
  ColSideColors = col.cell,
  
  scale = "row",
  
  margins = c(18, 12),
  
  cexCol = 1.3,
  cexRow = 1.3,
  
  srtCol = 90,
  
  labCol = colnames(highly_variable_lcpm),
  
  labRow = sparse_labels
)

dev.off()
#checking
dim(highly_variable_lcpm)
colnames(highly_variable_lcpm)



# Save the heatmap
png(file="High_var_genes.heatmap.png")
heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none", main="Top 500 most variable genes across samples",ColSideColors=col.cell,scale="row")
pdf("High_var_genes.heatmap15.pdf", width = 10, height = 15)
dev.off()
png("C:/Users/ASUS/Desktop/High_var_genes_heatmap.png")
dev.off()

# Apply normalisation to DGEList object
y <- calcNormFactors(y)
y$samples

par(mfrow=c(1,2))
plotMD(logcounts,column = 13)
abline(h=0,col="grey")
plotMD(logcounts,column = 13)
abline(h=0,col="grey")

par(mfrow=c(1,2))
plotMD(y,column = 3)
abline(h=0,col="grey")
plotMD(y,column = 13)
abline(h=0,col="grey")
group <- factor(sampleinfo$Group)
save(group, y, logcounts, sampleinfo, file = "day1objects.Rdata")  


#Differential expression with limma-voom

# load("day1objects.Rdata")
# objects()
# Look at group variable again
group
# Specify a design matrix without an intercept term
design <- model.matrix(~ 0 + group)

design
## Make the column names of the design matrix a bit nicer
colnames(design) <- levels(group)
design

par(mfrow=c(1,1))
v <- voom(y,design,plot = TRUE)
v
# What is contained in this object?
names(v)

#comparison
par(mfrow=c(1,2))
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2,main="Unnormalised logCPM")
## Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
boxplot(v$E, xlab="", ylab="Log2 counts per million",las=2,main="Voom transformed logCPM")
## Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(v$E),col="blue")

# Fit the linear model
fit <- lmFit(v)
names(fit)
head(fit)

colnames(design)

cont.matrix <- makeContrasts(
  Disease_vs_Control = Azoospermia - Control,
  levels = design
)

cont.matrix

fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
head(fit.cont)
dim(fit.cont)
summa.fit <- decideTests(fit.cont)
summary(summa.fit)


up_down_genes <- rownames(summa.fit)[summa.fit[, 1] != 0]  


up_down_status <- summa.fit[summa.fit[, 1] != 0, 1]

#
results_df <- data.frame(Gene = up_down_genes, Status = up_down_status)

# 
write.csv(results_df, "C:/Users/ASUS/Desktop/Up_Down_Genes.csv", row.names = FALSE)
read.csv("Up_Down_Genes.csv")

head(results_df)

# pvalue & gene to dataframe
head(fit.cont$p.value)
head(fit.cont$F)
head(fit.cont$F.p.value)

#saave
# Convert p-values to dataframe
pval_df <- data.frame(
  Gene = rownames(fit.cont$p.value),
  P_Value = fit.cont$p.value[,1]
)

# Save to Desktop
write.csv(
  pval_df,
  "C:/Users/ASUS/Desktop/pvalues.csv",
  row.names = FALSE
)

# We want to highlight the significant genes. We can get this from decideTests.
par(mfrow=c(1,1))
plotMD(fit.cont,coef=1,status=summa.fit[,"Disease_vs_Control"], values = c(-1, 1), hl.col=c("blue","red"))
svg("C:/Users/ASUS/Desktop/MD_plot.svg",
    width = 10,
    height = 8)

plotMD(
  fit.cont,
  coef = 1,
  status = summa.fit[,"Disease_vs_Control"],
  values = c(-1, 1),
  hl.col = c("blue", "red")
)

dev.off()

# For the volcano plot we have to specify how many of the top genes to highlight.
# We can also specify that we want to plot the gene symbol for the highlighted genes.
# let's highlight the top 100 most DE genes


p_ma <- ggplot(dis.vs.c2, aes(x = AveExpr, y = logFC, color = threshold)) +

  geom_point(alpha = 0.9, size = 1.8) + 
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.8) +
  
  #  (Darker Red and Blue)
  scale_color_manual(values = c("NS" = "grey80", "Up" = "firebrick3", "Down" = "dodgerblue4")) +
  
  labs(
    x = "Average log-expression (A)",
    y = expression(log[2]~"fold change (M)")
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",

    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.text = element_text(color = "black", face = "bold")
  )


print(p_ma)


ggsave("Figure8a_MA_v3_bold.pdf", p_ma, width = 6, height = 5)



####

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("EnhancedVolcano")



library(EnhancedVolcano)


#### another way to get data include pvalue & logfc & convert to data_frame
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggrepel)

# Extract full differential expression results
dis.vs.c2 <- topTable(
  fit.cont,
  coef = "Disease_vs_Control",
  number = Inf,
  sort.by = "P"
)

# Add gene identifiers and annotation
dis.vs.c2$ENSEMBL <- rownames(dis.vs.c2)
dis.vs.c2$ENTREZID <- y$genes$ENTREZID[match(dis.vs.c2$ENSEMBL, y$genes$ENSEMBL)]
dis.vs.c2$SYMBOL <- y$genes$SYMBOL[match(dis.vs.c2$ENSEMBL, y$genes$ENSEMBL)]
dis.vs.c2$GENENAME <- y$genes$GENENAME[match(dis.vs.c2$ENSEMBL, y$genes$ENSEMBL)]

# Replace missing gene symbols with ENSEMBL IDs
dis.vs.c2$SYMBOL[is.na(dis.vs.c2$SYMBOL) | dis.vs.c2$SYMBOL == ""] <- 
  dis.vs.c2$ENSEMBL[is.na(dis.vs.c2$SYMBOL) | dis.vs.c2$SYMBOL == ""]

# Save full result table
write.csv(dis.vs.c2, "limma_voom_DE_results_full.csv", row.names = FALSE)

# Filter significant DEGs using the final study threshold
sig_genes <- dis.vs.c2 %>%
  filter(adj.P.Val < 0.01 & abs(logFC) > 1)

# Save significant DEGs
write.csv(sig_genes, "limma_voom_DEGs_FDR0.01_logFC1.csv", row.names = FALSE)

# Check DEG count
nrow(sig_genes)

# ==============================================================================
# Expression summary for genes of interest
# ==============================================================================

genes_of_interest <- c(
  "ENSG00000243955",
  "ENSG00000244067",
  "ENSG00000174156",
  "ENSG00000242366",
  "ENSG00000241119",
  "ENSG00000242515",
  "ENSG00000137364",
  "ENSG00000197838",
  "ENSG00000198077",
  "ENSG00000205277",
  "ENSG00000157766",
  "ENSG00000197558",
  "ENSG00000181143"
)

# Use voom-transformed expression values
expr_matrix <- v$E

# Match sample information to expression matrix columns
rownames(sampleinfo) <- colnames(expr_matrix)
sampleinfo$Group <- factor(sampleinfo$Group)

# Subset expression matrix for genes of interest
expr_selected <- expr_matrix[rownames(expr_matrix) %in% genes_of_interest, , drop = FALSE]

# Define sample groups
samples_ctrl <- rownames(sampleinfo)[sampleinfo$Group == "Control"]
samples_case <- rownames(sampleinfo)[sampleinfo$Group == "Azoospermia"]

# Compute mean expression per group
expr_control <- expr_selected[, samples_ctrl, drop = FALSE]
expr_case <- expr_selected[, samples_case, drop = FALSE]

mean_ctrl <- rowMeans(expr_control)
mean_case <- rowMeans(expr_case)

# Build expression summary table
expression_summary <- data.frame(
  ENSEMBL = rownames(expr_selected),
  Mean_Control = round(mean_ctrl, 2),
  Mean_Azoospermia = round(mean_case, 2)
)

# Add annotation and DE statistics
expression_summary$SYMBOL <- y$genes$SYMBOL[
  match(expression_summary$ENSEMBL, y$genes$ENSEMBL)
]

expression_summary$logFC <- dis.vs.c2$logFC[
  match(expression_summary$ENSEMBL, dis.vs.c2$ENSEMBL)
]

expression_summary$adj.P.Val <- dis.vs.c2$adj.P.Val[
  match(expression_summary$ENSEMBL, dis.vs.c2$ENSEMBL)
]

# Reorder columns
expression_summary <- expression_summary %>%
  select(ENSEMBL, SYMBOL, Mean_Control, Mean_Azoospermia, logFC, adj.P.Val)

# Save group-level summary
write.csv(expression_summary, "Gene_Expression_by_Group.csv", row.names = FALSE)

# ==============================================================================
# Long-format expression table for per-sample visualization
# ==============================================================================

long_expr <- expr_selected %>%
  as.data.frame() %>%
  rownames_to_column("ENSEMBL") %>%
  pivot_longer(
    cols = -ENSEMBL,
    names_to = "Sample_ID",
    values_to = "Expression"
  ) %>%
  mutate(
    Expression = round(Expression, 2),
    Group = sampleinfo[Sample_ID, "Group"],
    SYMBOL = y$genes$SYMBOL[match(ENSEMBL, y$genes$ENSEMBL)]
  ) %>%
  select(ENSEMBL, SYMBOL, Sample_ID, Group, Expression)

# Save per-sample expression table
write.csv(long_expr, "Gene_Expression_per_Sample.csv", row.names = FALSE)

# Preview
head(long_expr)

# ==============================================================================
# Volcano plot
# ==============================================================================

# Prepare volcano plot data
dis.vs.c2 <- dis.vs.c2 %>%
  mutate(
    negLog10FDR = -log10(adj.P.Val),
    threshold = case_when(
      adj.P.Val < 0.01 & logFC > 1  ~ "Up",
      adj.P.Val < 0.01 & logFC < -1 ~ "Down",
      TRUE ~ "NS"
    )
  )

# Select top genes for labeling
label_genes <- dis.vs.c2 %>%
  filter(threshold != "NS") %>%
  group_by(threshold) %>%
  arrange(adj.P.Val, desc(abs(logFC)), .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  ungroup()

# Create volcano plot
p_volcano <- ggplot(dis.vs.c2, aes(x = logFC, y = negLog10FDR, color = threshold)) +
  geom_point(alpha = 0.75, size = 1.8) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "grey40") +
  geom_text_repel(
    data = label_genes,
    aes(label = SYMBOL),
    size = 3.5,
    box.padding = 0.35,
    point.padding = 0.25,
    segment.color = "grey50",
    max.overlaps = Inf
  ) +
  scale_color_manual(values = c(
    "NS" = "grey75",
    "Up" = "firebrick2",
    "Down" = "dodgerblue3"
  )) +
  labs(
    title = "Volcano plot: Azoospermia vs Control",
    x = expression(log[2]~fold~change),
    y = expression(-log[10]~FDR),
    caption = "Cut-offs: FDR < 0.01 and |log2FC| > 1"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(face = "bold")
  )

# Show plot
p_volcano

# Save volcano plot
ggsave("Figure8_volcano.pdf", p_volcano, width = 7, height = 5)
ggsave("Figure8_volcano.png", p_volcano, width = 7, height = 5, dpi = 600)

# Count genes by category
table(dis.vs.c2$threshold)
