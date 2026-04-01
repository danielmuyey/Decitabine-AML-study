# AML-transcriptomics-analysis
# Create a minimal, reusable analysis workflow template (R)
# This writes a script you can upload to Zenodo as the reproducible pipeline entry point.

# Cell 1: load libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(pheatmap)
})

cat('Libraries loaded.\
')
# Cell 1 (retry): load libraries that are guaranteed, install missing as needed
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
})

if (!requireNamespace('pheatmap', quietly = TRUE)) {
  cat('Installing pheatmap...\
')
  install.packages('pheatmap', repos = 'https://cloud.r-project.org')
}

suppressPackageStartupMessages({
  library(pheatmap)
})

cat('Libraries loaded.\
')
# Cell 3: Build analysis-ready count matrices + minimal sample metadata

make_counts_from_GSE317967 <- function(df) {
  gene_col <- 'Gene_ID'
  count_cols <- grep('read count', colnames(df), value = TRUE)
  counts <- df[, c(gene_col, count_cols)]
  counts[[gene_col]] <- as.character(counts[[gene_col]])
  # coerce to numeric
  counts[count_cols] <- lapply(counts[count_cols], function(x) as.numeric(as.character(x)))
  mat <- as.matrix(counts[, count_cols, drop = FALSE])
  rownames(mat) <- counts[[gene_col]]
  
  sample_names <- sub(':read count', '', count_cols, fixed = TRUE)
  colnames(mat) <- sample_names
  
  # crude grouping: AML vs ctrl based on substring
  group <- ifelse(grepl('ctrl', sample_names, ignore.case = TRUE), 'CTRL', 'AML')
  meta <- data.frame(sample = sample_names, group = factor(group, levels = c('CTRL','AML')))
  list(counts = mat, meta = meta)
}

make_counts_from_first_col <- function(df, gene_colname = NULL) {
  if (is.null(gene_colname)) gene_colname <- colnames(df)[1]
  gene <- as.character(df[[gene_colname]])
  out <- df
  out[[gene_colname]] <- NULL
  # keep only numeric-ish columns
  out <- out[, sapply(out, function(x) {
    x2 <- suppressWarnings(as.numeric(as.character(x)))
    sum(!is.na(x2)) > 0.9 * length(x2)
  }), drop = FALSE]
  out[] <- lapply(out, function(x) as.numeric(as.character(x)))
  mat <- as.matrix(out)
  rownames(mat) <- gene
  list(counts = mat)
}

# GSE317967 (AML vs control)
obj_317967 <- make_counts_from_GSE317967(GSE317967_raw)
counts_317967 <- obj_317967$counts
meta_317967 <- obj_317967$meta

cat('GSE317967 count matrix dims (genes x samples):\
')
print(dim(counts_317967))
cat('\
GSE317967 metadata head:\
')
print(head(meta_317967))

# GSE289928 (Decitabine vs DMSO in KG1A - inferred from column names)
# This table has annotation columns; counts start at columns with patterns hsa_
count_cols_289928 <- grep('^hsa_', colnames(GSE289928_raw), value = TRUE)
counts_289928_df <- GSE289928_raw[, c('symbol', count_cols_289928), drop = FALSE]
counts_289928_df[count_cols_289928] <- lapply(counts_289928_df[count_cols_289928], function(x) as.numeric(as.character(x)))
counts_289928 <- as.matrix(counts_289928_df[, count_cols_289928, drop = FALSE])
rownames(counts_289928) <- as.character(counts_289928_df$symbol)
colnames(counts_289928) <- count_cols_289928

meta_289928 <- data.frame(
  sample = count_cols_289928,
  treatment = ifelse(grepl('Deci', count_cols_289928, ignore.case = TRUE), 'Decitabine', 'DMSO')
) %>%
  mutate(treatment = factor(treatment, levels = c('DMSO','Decitabine')))

cat('\
GSE289928 count matrix dims:\
')
print(dim(counts_289928))
cat('\
GSE289928 metadata:\
')
print(meta_289928)

# GSE313457 (mouse macrophage polarization style dataset; first column is gene symbol)
obj_313457 <- make_counts_from_first_col(GSE313457_raw, gene_colname = colnames(GSE313457_raw)[1])
counts_313457 <- obj_313457$counts
cat('\
GSE313457 count matrix dims:\
')
print(dim(counts_313457))

# GSE272757 validation cohort (AML only; first column is ID)
obj_272757 <- make_counts_from_first_col(GSE272757_raw, gene_colname = 'ID')
counts_272757 <- obj_272757$counts
cat('\
GSE272757 count matrix dims:\
')
print(dim(counts_272757))

# quick sanity checks
cat('\
Library size summary (GSE317967):\
')
print(summary(colSums(counts_317967)))
# Cell 4: Differential expression with edgeR (fast, no need for DESeq2) + plots

suppressPackageStartupMessages({
  library(edgeR)
})

run_edger_de <- function(count_mat, group_factor, group_name = 'group', ref_level = NULL) {
  stopifnot(ncol(count_mat) == length(group_factor))
  if (!is.null(ref_level)) group_factor <- relevel(factor(group_factor), ref = ref_level)
  
  y <- DGEList(counts = count_mat, group = group_factor)
  keep <- filterByExpr(y, group = group_factor)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  design <- model.matrix(~ group_factor)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, coef = 2)
  tt <- topTags(qlf, n = Inf)$table
  tt$gene <- rownames(tt)
  tt
}

# GSE317967: AML vs CTRL
res_317967 <- run_edger_de(counts_317967, meta_317967$group, ref_level = 'CTRL')

cat('GSE317967 DE result head:\
')
print(head(res_317967))

# Volcano plot
vol_317967 <- res_317967 %>%
  mutate(neglogFDR = -log10(FDR), sig = FDR < 0.05 & abs(logFC) >= 1)

p_vol_317967 <- ggplot(vol_317967, aes(x = logFC, y = neglogFDR)) +
  geom_point(aes(color = sig), alpha = 0.6, size = 1.2) +
  scale_color_manual(values = c('FALSE' = 'grey70', 'TRUE' = '#d62728')) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = 'grey50') +
  geom_vline(xintercept = c(-1, 1), linetype = 2, color = 'grey50') +
  theme_minimal(base_size = 12) +
  labs(title = 'GSE317967: AML vs CTRL (edgeR QL)', x = 'log2 fold-change (AML/CTRL)', y = '-log10(FDR)', color = 'Significant')

print(p_vol_317967)

# PCA on logCPM
y_317967 <- DGEList(counts = counts_317967)
keep_317967 <- filterByExpr(y_317967, group = meta_317967$group)
y_317967 <- y_317967[keep_317967, , keep.lib.sizes = FALSE]
y_317967 <- calcNormFactors(y_317967)
logcpm_317967 <- cpm(y_317967, log = TRUE, prior.count = 1)

pca_317967 <- prcomp(t(logcpm_317967), scale. = TRUE)
pca_df_317967 <- data.frame(
  sample = rownames(pca_317967$x),
  PC1 = pca_317967$x[,1],
  PC2 = pca_317967$x[,2]
) %>% left_join(meta_317967, by = 'sample')

p_pca_317967 <- ggplot(pca_df_317967, aes(PC1, PC2, color = group, label = sample)) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(size = 3, max.overlaps = 20) +
  theme_minimal(base_size = 12) +
  labs(title = 'GSE317967: PCA on logCPM', color = 'Group')

print(p_pca_317967)

# GSE289928: Decitabine vs DMSO
res_289928 <- run_edger_de(counts_289928, meta_289928$treatment, ref_level = 'DMSO')
cat('\
GSE289928 DE result head:\
')
print(head(res_289928))

vol_289928 <- res_289928 %>% mutate(neglogFDR = -log10(FDR), sig = FDR < 0.05 & abs(logFC) >= 1)

p_vol_289928 <- ggplot(vol_289928, aes(x = logFC, y = neglogFDR)) +
  geom_point(aes(color = sig), alpha = 0.7, size = 1.4) +
  scale_color_manual(values = c('FALSE' = 'grey70', 'TRUE' = '#1f77b4')) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = 'grey50') +
  geom_vline(xintercept = c(-1, 1), linetype = 2, color = 'grey50') +
  theme_minimal(base_size = 12) +
  labs(title = 'GSE289928: Decitabine vs DMSO (edgeR QL)', x = 'log2 fold-change (Decitabine/DMSO)', y = '-log10(FDR)', color = 'Significant')

print(p_vol_289928)

# Heatmap of top DE genes for each (top 30 by FDR)
plot_top_heatmap <- function(logcpm_mat, meta_df, group_col, title) {
  top_genes <- meta_df$gene[1:min(30, nrow(meta_df))]
  hm <- logcpm_mat[top_genes, , drop = FALSE]
  # z-score by gene
  hm_z <- t(scale(t(hm)))
  ann <- data.frame(group = meta_df2[[group_col]])
}

# Simple heatmap for GSE317967
sig_317967 <- res_317967 %>% arrange(FDR)
sel_317967 <- sig_317967$gene[1:min(30, nrow(sig_317967))]
hm_317967 <- logcpm_317967[sel_317967, , drop = FALSE]
hm_317967_z <- t(scale(t(hm_317967)))
ann_col_317967 <- data.frame(group = meta_317967$group)
rownames(ann_col_317967) <- meta_317967$sample

pheatmap(hm_317967_z,
         annotation_col = ann_col_317967,
         show_colnames = TRUE,
         show_rownames = FALSE,
         main = 'GSE317967: Top 30 DE genes (z-scored logCPM)')

# Simple heatmap for GSE289928
y_289928 <- DGEList(counts = counts_289928)
keep_289928 <- filterByExpr(y_289928, group = meta_289928$treatment)
y_289928 <- y_289928[keep_289928, , keep.lib.sizes = FALSE]
y_289928 <- calcNormFactors(y_289928)
logcpm_289928 <- cpm(y_289928, log = TRUE, prior.count = 1)

sig_289928 <- res_289928 %>% arrange(FDR)
sel_289928 <- sig_289928$gene[1:min(30, nrow(sig_289928))]
hm_289928 <- logcpm_289928[sel_289928, , drop = FALSE]
hm_289928_z <- t(scale(t(hm_289928)))
ann_col_289928 <- data.frame(treatment = meta_289928$treatment)
rownames(ann_col_289928) <- meta_289928$sample

pheatmap(hm_289928_z,
         annotation_col = ann_col_289928,
         show_colnames = TRUE,
         show_rownames = FALSE,
         main = 'GSE289928: Top 30 DE genes (z-scored logCPM)')
         

# Cell 5: Install missing Bioconductor packages (edgeR) and plotting helpers

if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')

pkgs_needed <- c('edgeR','pheatmap','ggrepel')
for (p in pkgs_needed) {
  if (!requireNamespace(p, quietly = TRUE)) {
    BiocManager::install(p, ask = FALSE, update = FALSE)
  }
}

suppressPackageStartupMessages({
  library(edgeR)
  library(pheatmap)
  library(ggrepel)
})

cat('Packages loaded: edgeR, pheatmap, ggrepel\
')

# Cell 6: edgeR differential expression + PCA + heatmaps (using already-loaded matrices)

suppressPackageStartupMessages({
  library(edgeR)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ggrepel)
  library(pheatmap)
})

run_edger_ql <- function(count_mat, group_factor, ref_level) {
  group_factor <- relevel(factor(group_factor), ref = ref_level)
  y <- DGEList(counts = count_mat, group = group_factor)
  keep <- filterByExpr(y, group = group_factor)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  design <- model.matrix(~ group_factor)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, coef = 2)
  tab <- topTags(qlf, n = Inf)$table
  tab$gene <- rownames(tab)
  list(tab = tab, y = y)
}

make_volcano <- function(res_tab, title_text, fc_cut = 1, fdr_cut = 0.05) {
  df <- res_tab %>%
    mutate(neglogFDR = -log10(FDR), sig = (FDR < fdr_cut) & (abs(logFC) >= fc_cut))
  ggplot(df, aes(logFC, neglogFDR)) +
    geom_point(aes(color = sig), alpha = 0.65, size = 1.2) +
    scale_color_manual(values = c('FALSE' = 'grey75', 'TRUE' = '#d62728')) +
    geom_hline(yintercept = -log10(fdr_cut), linetype = 2, color = 'grey55') +
    geom_vline(xintercept = c(-fc_cut, fc_cut), linetype = 2, color = 'grey55') +
    theme_minimal(base_size = 12) +
    labs(title = title_text, x = 'log2 fold-change', y = '-log10(FDR)', color = paste0('FDR<',fdr_cut,' & |logFC|≥',fc_cut))
}

make_pca_plot <- function(y_obj, meta_df, group_col, title_text) {
  logcpm <- cpm(y_obj, log = TRUE, prior.count = 1)
  pca <- prcomp(t(logcpm), scale. = TRUE)
  pc_df <- data.frame(sample = rownames(pca$x), PC1 = pca$x[,1], PC2 = pca$x[,2])
  pc_df <- pc_df %>% left_join(meta_df, by = 'sample')
  ggplot(pc_df, aes(PC1, PC2, color = .data[[group_col]], label = sample)) +
    geom_point(size = 3) +
    geom_text_repel(size = 3, max.overlaps = 30) +
    theme_minimal(base_size = 12) +
    labs(title = title_text, color = group_col)
}

make_top_heatmap <- function(y_obj, res_tab, meta_df, group_col, title_text, top_n = 30) {
  logcpm <- cpm(y_obj, log = TRUE, prior.count = 1)
  top_genes <- res_tab %>% arrange(FDR) %>% slice_head(n = min(top_n, n())) %>% pull(gene)
  hm <- logcpm[top_genes, , drop = FALSE]
  hm_z <- t(scale(t(hm)))
  ann_col <- data.frame(group = meta_df[[group_col]])
  rownames(ann_col) <- meta_df$sample
  pheatmap(hm_z,
           annotation_col = ann_col,
           show_colnames = TRUE,
           show_rownames = FALSE,
           fontsize_col = 10,
           main = title_text)
}

# --- GSE317967 (AML vs CTRL) ---
edger_317967 <- run_edger_ql(counts_317967, meta_317967$group, ref_level = 'CTRL')
res_317967 <- edger_317967$tab

cat('GSE317967 edgeR results (head):\
')
print(head(res_317967, 10))

# Summary counts
sig_counts_317967 <- res_317967 %>%
  summarize(
    tested_genes = n(),
    fdr_lt_0.05 = sum(FDR < 0.05),
    fdr_lt_0.05_abslogfc_ge_1 = sum(FDR < 0.05 & abs(logFC) >= 1),
    up = sum(FDR < 0.05 & logFC >= 1),
    down = sum(FDR < 0.05 & logFC <= -1)
  )
cat('\
GSE317967 significant-gene summary:\
')
print(sig_counts_317967)

p_vol_317967 <- make_volcano(res_317967, 'GSE317967: AML vs CTRL (edgeR QL)')
print(p_vol_317967)

p_pca_317967 <- make_pca_plot(edger_317967$y, meta_317967, 'group', 'GSE317967: PCA (logCPM)')
print(p_pca_317967)

make_top_heatmap(edger_317967$y, res_317967, meta_317967, 'group', 'GSE317967: Top 30 DE genes (z-scored logCPM)')

# --- GSE289928 (Decitabine vs DMSO) ---
edger_289928 <- run_edger_ql(counts_289928, meta_289928$treatment, ref_level = 'DMSO')
res_289928 <- edger_289928$tab

cat('\
GSE289928 edgeR results (head):\
')
print(head(res_289928, 10))

sig_counts_289928 <- res_289928 %>%
  summarize(
    tested_genes = n(),
    fdr_lt_0.05 = sum(FDR < 0.05),
    fdr_lt_0.05_abslogfc_ge_1 = sum(FDR < 0.05 & abs(logFC) >= 1),
    up = sum(FDR < 0.05 & logFC >= 1),
    down = sum(FDR < 0.05 & logFC <= -1)
  )
cat('\
GSE289928 significant-gene summary:\
')
print(sig_counts_289928)

p_vol_289928 <- make_volcano(res_289928, 'GSE289928: Decitabine vs DMSO (edgeR QL)')
print(p_vol_289928)

p_pca_289928 <- make_pca_plot(edger_289928$y, meta_289928, 'treatment', 'GSE289928: PCA (logCPM)')
print(p_pca_289928)

make_top_heatmap(edger_289928$y, res_289928, meta_289928, 'treatment', 'GSE289928: Top 30 DE genes (z-scored logCPM)')
# Cell 7: Fix heatmap helper and rerun remaining visuals + complete second dataset outputs

suppressPackageStartupMessages({
  library(edgeR)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ggrepel)
  library(pheatmap)
})

make_top_heatmap <- function(y_obj, res_tab, meta_df, group_col, title_text, top_n = 30) {
  logcpm <- cpm(y_obj, log = TRUE, prior.count = 1)
  top_n2 <- min(top_n, nrow(res_tab))
  top_genes <- res_tab %>% arrange(FDR) %>% slice_head(n = top_n2) %>% pull(gene)
  hm <- logcpm[top_genes, , drop = FALSE]
  hm_z <- t(scale(t(hm)))
  ann_col <- data.frame(group = meta_df[[group_col]])
  rownames(ann_col) <- meta_df$sample
  pheatmap(hm_z,
           annotation_col = ann_col,
           show_colnames = TRUE,
           show_rownames = FALSE,
           fontsize_col = 10,
           main = title_text)
}

# Recreate the two edgeR result objects quickly (counts/meta already in memory)
run_edger_ql <- function(count_mat, group_factor, ref_level) {
  group_factor <- relevel(factor(group_factor), ref = ref_level)
  y <- DGEList(counts = count_mat, group = group_factor)
  keep <- filterByExpr(y, group = group_factor)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  design <- model.matrix(~ group_factor)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, coef = 2)
  tab <- topTags(qlf, n = Inf)$table
  tab$gene <- rownames(tab)
  list(tab = tab, y = y)
}

make_volcano <- function(res_tab, title_text, fc_cut = 1, fdr_cut = 0.05) {
  df <- res_tab %>% mutate(neglogFDR = -log10(FDR), sig = (FDR < fdr_cut) & (abs(logFC) >= fc_cut))
  ggplot(df, aes(logFC, neglogFDR)) +
    geom_point(aes(color = sig), alpha = 0.65, size = 1.2) +
    scale_color_manual(values = c('FALSE' = 'grey75', 'TRUE' = '#d62728')) +
    geom_hline(yintercept = -log10(fdr_cut), linetype = 2, color = 'grey55') +
    geom_vline(xintercept = c(-fc_cut, fc_cut), linetype = 2, color = 'grey55') +
    theme_minimal(base_size = 12) +
    labs(title = title_text, x = 'log2 fold-change', y = '-log10(FDR)', color = paste0('FDR<',fdr_cut,' & |logFC|≥',fc_cut))
}

make_pca_plot <- function(y_obj, meta_df, group_col, title_text) {
  logcpm <- cpm(y_obj, log = TRUE, prior.count = 1)
  pca <- prcomp(t(logcpm), scale. = TRUE)
  pc_df <- data.frame(sample = rownames(pca$x), PC1 = pca$x[,1], PC2 = pca$x[,2])
  pc_df <- pc_df %>% left_join(meta_df, by = 'sample')
  ggplot(pc_df, aes(PC1, PC2, color = .data[[group_col]], label = sample)) +
    geom_point(size = 3) +
    geom_text_repel(size = 3, max.overlaps = 30) +
    theme_minimal(base_size = 12) +
    labs(title = title_text, color = group_col)
}

# GSE317967
edger_317967 <- run_edger_ql(counts_317967, meta_317967$group, ref_level = 'CTRL')
res_317967 <- edger_317967$tab

# Heatmap (previously failed)
make_top_heatmap(edger_317967$y, res_317967, meta_317967, 'group', 'GSE317967: Top 30 DE genes (z-scored logCPM)')

# GSE289928
edger_289928 <- run_edger_ql(counts_289928, meta_289928$treatment, ref_level = 'DMSO')
res_289928 <- edger_289928$tab

cat('GSE289928 edgeR results (head):\
')
print(head(res_289928, 10))

sig_counts_289928 <- res_289928 %>%
  summarize(
    tested_genes = n(),
    fdr_lt_0.05 = sum(FDR < 0.05),
    fdr_lt_0.05_abslogfc_ge_1 = sum(FDR < 0.05 & abs(logFC) >= 1),
    up = sum(FDR < 0.05 & logFC >= 1),
    down = sum(FDR < 0.05 & logFC <= -1)
  )
cat('\
GSE289928 significant-gene summary:\
')
print(sig_counts_289928)

print(make_volcano(res_289928, 'GSE289928: Decitabine vs DMSO (edgeR QL)'))
print(make_pca_plot(edger_289928$y, meta_289928, 'treatment', 'GSE289928: PCA (logCPM)'))

make_top_heatmap(edger_289928$y, res_289928, meta_289928, 'treatment', 'GSE289928: Top 30 DE genes (z-scored logCPM)')


workflow_lines <- c(
  "#!/usr/bin/env Rscript",
  "",
  "suppressPackageStartupMessages({",
  "  library(tidyverse)",
  "  # library(DESeq2)   # uncomment if you run DESeq2",
  "  # library(immunedeconv)  # optional wrapper for CIBERSORT/other methods",
  "})",
  "",
  "message('Reproducible AML immune-context workflow starting...')",
  "",
  "# ---- helper readers ----",
  "read_counts_table <- function(path) {",
  "  ext <- tolower(tools::file_ext(path))",
  "  if (ext == 'csv') {",
  "    df <- read.csv(path, check.names = FALSE)",
  "  } else {",
  "    df <- read.table(path, header = TRUE, sep = '\\	', check.names = FALSE, quote = '', comment.char = '')",
  "  }",
  "  df",
  "}",
  "",
  "# ---- inputs (edit paths as needed) ----",
  "in_files <- list(",
  "  GSE317967 = 'GSE317967_AML-ctrl.csv',",
  "  GSE313457 = 'GSE313457_deSeq2_counts.txt',",
  "  GSE289928 = 'GSE289928_read_counts.txt',",
  "  GSE272757 = 'GSE272757_RJAML_147_Rawcount.txt'",
  ")",
  "",
  "# ---- load dataset (example: GSE317967) ----",
  "stopifnot(file.exists(in_files$GSE317967))",
  "gse317967_raw <- read.csv(in_files$GSE317967, check.names = FALSE)",
  "message('Loaded: ', in_files$GSE317967)",
  "print(head(gse317967_raw[, 1:min(8, ncol(gse317967_raw))]))",
  "",
  "# ---- build a counts matrix (edit this block to match your column naming) ----",
  "count_cols <- grep('read count', colnames(gse317967_raw), value = TRUE)",
  "stopifnot(length(count_cols) > 0)",
  "",
  "counts_tbl <- gse317967_raw %>%",
  "  select(Gene_ID, all_of(count_cols)) %>%",
  "  rename(gene_id = Gene_ID)",
  "",
  "counts_mat <- counts_tbl %>%",
  "  column_to_rownames('gene_id') %>%",
  "  as.matrix()",
  "storage.mode(counts_mat) <- 'numeric'",
  "message('Counts matrix dims: ', paste(dim(counts_mat), collapse = ' x '))",
  "",
  "# ---- QC: library sizes ----",
  "lib_sizes <- colSums(counts_mat)",
  "qc_df <- tibble(sample = names(lib_sizes), library_size = as.numeric(lib_sizes))",
  "print(qc_df %>% arrange(desc(library_size)) %>% head(10))",
  "",
  "# ---- (placeholder) differential expression ----",
  "# Define sample metadata with group labels (AML vs control) and covariates.",
  "# Then run DESeq2 or edgeR and export a DE table.",
  "",
  "# ---- immune deconvolution (CIBERSORT/LM22) ----",
  "# Recommended: run CIBERSORT externally (or via an R wrapper) to obtain LM22 fractions and a P-value per sample.",
  "# Save a table with samples as rows and 22 cell fractions as columns plus a column named 'P.value'.",
  "# Then filter P.value < 0.05 and correlate candidate genes with immune fractions.",
  "",
  "# Example expected file:",
  "# cibersort_res <- read.csv('cibersort_LM22_results.csv', check.names = FALSE)",
  "# cibersort_filt <- cibersort_res %>% filter(P.value < 0.05)",
  "",
  "# ---- save minimal outputs ----",
  "dir.create('results', showWarnings = FALSE)",
  "write.csv(qc_df, file = file.path('results', 'qc_library_sizes.csv'), row.names = FALSE)",
  "message('Done. Outputs written to results/.')",
  ""
)

writeLines(workflow_lines, 'zenodo_workflow_template.R')
cat('Wrote file: zenodo_workflow_template.R\
')


