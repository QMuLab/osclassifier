#' Gene Modules for Osteosarcoma Transcriptional Programs
#'
#' A named list of curated gene sets representing four osteosarcoma
#' malignant-cell transcriptional programs:
#' - Proliferating-like
#' - Osteoblast-like
#' - Chondroblast-like
#' - Fibroblast-like
#'
#' These gene sets are used by `compute_module_scores()` to calculate
#' module scores and classify samples according to their dominant
#' transcriptional program.
#'
#' The gene sets were derived from consensus non-negative matrix
#' factorization (cNMF) analysis and represent recurrent transcriptional
#' programs identified in osteosarcoma malignant cells.
#'
#' @format A named list with 4 elements. Each element is a character vector
#'   of gene symbols.
#' @usage data(gene_list)
#' @source Derived from cNMF analysis of osteosarcoma single-cell RNA-seq data.
"gene_list"
