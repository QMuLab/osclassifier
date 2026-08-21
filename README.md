# osclassifier <img src="https://img.shields.io/badge/R-package-blue.svg" align="right" height="30" />

A lightweight R package to classify osteosarcoma samples into molecular subtypes using gene module scores. Designed for bulk RNA-seq data with CPM, TPM, RPKM or FPKM values.

It also serves as a general sample state classifier, allowing users to define their own gene sets to classify custom samples, compute module scores, identify mixed transcriptional states, and visualize the results as a heatmap.

```r
# Installation

# If devtools not installed:
install.packages("devtools")

# Then install from GitHub
devtools::install_github("QMuLab/osclassifier")


# Example

library(osclassifier)

data(gene_list)

# ---- Optional: Use your own gene list ----
# Must be a named list: names = module names,
# elements = character vectors of gene symbols
# Example:
# my_gene_list <- list(
#  "SubtypeA" = c("Gene1", "Gene2", "Gene3"),
#  "SubtypeB" = c("Gene4", "Gene5", "Gene6"),
#  "SubtypeC" = c("Gene7", "Gene8", "Gene9"),
#  "SubtypeD" = c("Gene10", "Gene11", "Gene12"))



# Simulated expression matrix:
# rows = gene symbols, columns = samples (for demonstration only)
# Input must be CPM, TPM, RPKM or FPKM values
expr <- matrix(
  runif(length(unique(unlist(gene_list))) * 100),
  nrow = length(unique(unlist(gene_list))),
  ncol = 100,
  dimnames = list(
    unique(unlist(gene_list)),
    paste0("Sample", 1:100)
  )
)

# ---- Step 0: Normalize ----
# If your data has already been log-transformed and normalized,
# skip this step.
expr_norm <- log1p(expr)
expr_norm <- scale(t(scale(t(expr_norm))))
expr_norm[is.na(expr_norm)] <- 0


# ---- Step 1: Compute module scores and classify samples ----
res_scores <- compute_module_scores(
  expr = expr_norm,
  gene_list = gene_list,
  module_order = c(
    "Proliferating-like",
    "Osteoblast-like",
    "Chondroblast-like",
    "Fibroblast-like"
  ),
  mixed_threshold = 0.1
)

# Relative dominance is calculated as:
# (Top1 - Top2) / (Top1 - MinScore)
#
# Samples with RelativeDominance <= mixed_threshold
# are classified as "Mixed".
#
# mixed_threshold is adjustable.

module_scores <- res_scores$scores

# SimplicityScore is calculated internally using the gap method.
# It is used only for ordering samples within each subtype
# and is not displayed in the heatmap.

table(module_scores$TopCluster)


# ---- Step 2: Order samples and prepare heatmap inputs ----
hp <- order_and_prepare_heatmap(
  scores = module_scores,
  modules_final = res_scores$modules_final
)

# hp$mat: modules × samples matrix for plotting
# hp$anno_col / hp$ann_colors: subtype annotation and colors
# hp$order: sample order

module_scores <- hp$scores_ordered


# ---- Step 3: Draw heatmap ----
plot_module_heatmap(
  hp,
  show_colnames = TRUE, fontsize_col = 6
)
