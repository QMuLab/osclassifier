
# osclassifier <img src="https://img.shields.io/badge/R-package-blue.svg" align="right" height="30" />

A lightweight R package to classify osteosarcoma samples into molecular subtypes using gene module scores. Designed for bulk RNA-seq data with CPM, TPM, RPKM or FPKM values.
It also serves as a general sample state classifier, allowing users to define their own gene sets to classify custom samples, compute similarity scores, and visualize the results as a heatmap.
```r
# Installation

# If devtools not installed:
install.packages("devtools")

# Then install from local directory
devtools::install_github("QMuLab/osclassifier")


# Example

library(osclassifier)

data(gene_list)

# ---- Optional: Use your own gene list ----
# Must be a named list: names = module names, elements = character vectors of gene symbols
# Example:
# my_gene_list <- list(
#   "SubtypeA" = c("Gene1", "Gene2", "Gene3"),
#   "SubtypeB" = c("Gene4", "Gene5", "Gene6")
# )
# Ensure that gene symbols match row names of your expression matrix
# gene_list <- my_gene_list  # Uncomment to use your own list

# Simulated expression matrix: rows = gene symbols, columns = samples (for demonstration only)
# Input must be raw CPM, TPM, RPKM or FPKM values (NOT log-transformed or normalized)
expr <- matrix(
  runif(length(unique(unlist(gene_list))) * 100),
  nrow = length(unique(unlist(gene_list))),
  ncol = 100,
  dimnames = list(unique(unlist(gene_list)), paste0("Sample", 1:100))
)

# ---- Step 0: Normalize ----
# If your data has already been log-transformed and normalization, skip this step.
expr_norm <- log1p(expr)                  # log1p transform
expr_norm <- scale(t(scale(t(expr_norm)))       #Normalization 
expr_norm[is.na(expr_norm)] <- 0           # replace NA from zero-variance genes

# ---- Step 1: Compute module scores & TopCluster ----
res_scores <- compute_module_scores(
  expr = expr_norm,
  gene_list = gene_list,
  module_order = c("Proliferating-like","Osteoblast-like","Chondroblast-like","Fibroblast-like")
)
# res_scores$scores: data.frame with per-module scores + TopCluster
# res_scores$modules_final: final module order

# ---- Step 2: Add SimplicityScore ----
module_scores <- add_simplicity_scores(
  scores        = res_scores$scores,
  modules_final = res_scores$modules_final,
  method        = "gap"   # or "entropy"
)
#"gap" (gap-based): Measures how dominant the top module score is compared to others. Large when one module is much higher than the rest.
#"entropy" (purity-based): Uses Shannon entropy of normalized module scores (softmax). Large when one module dominates, small when scores are evenly distributed.

# ---- Step 3: Order samples & prepare heatmap inputs ----
hp <- order_and_prepare_heatmap(
  scores = module_scores,
  modules_final = res_scores$modules_final
)
# hp$mat:   modules × samples matrix for plotting
# hp$anno_col / hp$ann_colors: column annotations and colors
# hp$order: sample order

# Replace module_scores with the ordered version
module_scores <- hp$scores_ordered

# ---- Step 4: Draw heatmap (pure plotting API) ----
plot_module_heatmap(hp, show_colnames = TRUE)


