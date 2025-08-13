#' Gene Modules for Osteosarcoma Subtypes
#'
#' A named list of curated gene sets representing four osteosarcoma molecular subtypes:
#' - Proliferating-like
#' - Osteoblast-like
#' - Chondroblast-like
#' - Fibroblast-like
#'
#' These gene sets are used by `classify_plot_modules()` to compute module scores and assign samples to the most representative subtype.
#' Genes were derived from consensus non-negative matrix factorization (cNMF) analysis, capturing subtype-specific expression programs.
#'
#' @format A named list with 4 elements. Each element is a character vector of gene symbols.
#' @usage data(gene_list)
#' @source Internal curation.
"gene_list"
