#' Compute module scores and finalize module order
#' @param expr matrix (genes x samples), assumed already log1p + row-wise z-scored
#' @param gene_list named list of gene vectors
#' @param module_order optional desired order of modules
#' @return list(scores=data.frame, modules_final=character)
#' @export
compute_module_scores <- function(expr, gene_list,
                                  module_order = c('Proliferating-like','Osteoblast-like','Chondroblast-like','Fibroblast-like')) {
  if (!is.matrix(expr) || !is.numeric(expr)) stop("'expr' must be a numeric matrix (genes x samples).")
  if (is.null(rownames(expr)) || is.null(colnames(expr))) stop("'expr' must have row and column names.")
  if (!is.list(gene_list) || is.null(names(gene_list))) stop("'gene_list' must be a NAMED list of gene vectors.")
  modules <- names(gene_list)
  scores <- data.frame(row.names = colnames(expr))
  for (m in modules) {
    g <- intersect(gene_list[[m]], rownames(expr))
    scores[[m]] <- if (length(g) > 0) colMeans(expr[g, , drop = FALSE], na.rm = TRUE) else NA_real_
  }
  ordered_mods  <- intersect(module_order, modules)
  modules_final <- c(ordered_mods, setdiff(modules, ordered_mods))
  scores$TopCluster <- apply(scores[, modules_final, drop = FALSE], 1, function(v) names(v)[which.max(v)])
  list(scores = scores, modules_final = modules_final)
}

#' Add SimplicityScore
#'
#' @param scores data.frame from compute_module_scores()
#' @param module_order character vector of module names
#' @param method character, one of:
#'   \itemize{
#'     \item \code{"gap"}    Gap-based method (legacy). Measures the difference between
#'           the dominant module score and others (ADDS), subtracts internal differences
#'           among non-dominant modules (ADNS), and applies a range-based correction.
#'     \item \code{"entropy"} Entropy-based method. Normalizes scores to probabilities,
#'           computes Shannon entropy, and converts it to a normalized "purity" score:
#'           \eqn{1 - H / \log(N)}, where \eqn{N} = number of modules.
#'   }
#' @return data.frame with an added \code{SimplicityScore} column
#' @export
add_simplicity_scores <- function(scores, modules_final, method = c("gap", "entropy")) {
  method <- match.arg(method)

  if (method == "gap") {
    calc_simplicity <- function(x) {
      r <- sort(as.numeric(x), decreasing = TRUE)
      N <- length(r)
      if (N < 3) return(max(r) - min(r))
      adds <- sum(r[1] - r[2:N])
      mlen <- N - 2
      if (mlen <= 0) {
        adns <- 0
      } else {
        M1 <- matrix(rep(r[2:(N-1)], each = mlen), nrow = mlen)
        M2 <- matrix(rep(r[3:N], times = mlen), nrow = mlen)
        ut <- upper.tri(matrix(0, mlen, mlen), diag = TRUE)
        adns <- sum((M1 - M2)[ut])
      }
      correction <- if (N > 1) (r[1] - r[N]) / (N - 1) else 1
      (adds - adns) * correction
    }
    scores$SimplicityScore <- apply(scores[, modules_final, drop = FALSE], 1, calc_simplicity)

  } else if (method == "entropy") {
    calc_entropy_simplicity <- function(x) {
      x_shift <- x - min(x, na.rm = TRUE)
      if (sum(x_shift) == 0) return(0)
      p <- x_shift / sum(x_shift)
      p <- p[p > 0]  # avoid log(0)
      H <- -sum(p * log(p))
      1 - H / log(length(modules_final))  # normalized simplicity
    }
    scores$SimplicityScore <- apply(scores[, modules_final, drop = FALSE], 1, calc_entropy_simplicity)
  }

  scores
}

#' Order samples and assemble heatmap inputs
#' @param scores data.frame after add_simplicity_scores()
#' @param modules_final character vector of module names
#' @return list(mat, anno_col, ann_colors, order, scores_ordered)
#' @importFrom grDevices colorRampPalette
#' @export
order_and_prepare_heatmap <- function(scores, modules_final) {
  ord   <- order(factor(scores$TopCluster, levels = modules_final), -scores$SimplicityScore)
  s_ord <- scores[ord, , drop = FALSE]
  mat   <- t(s_ord[, modules_final, drop = FALSE])
  simp  <- s_ord$SimplicityScore
  rng   <- range(simp, na.rm = TRUE)
  simp01 <- if (diff(rng) > 0) (simp - rng[1]) / diff(rng) else rep(0, length(simp))
  anno_col   <- data.frame(`Simplicity Score` = simp01, check.names = FALSE, row.names = colnames(mat))
  ann_colors <- list(`Simplicity Score` = grDevices::colorRampPalette(c("#E0F3FF", "deepskyblue4"))(100))
  list(mat = mat, anno_col = anno_col, ann_colors = ann_colors, order = ord, scores_ordered = s_ord)
}

#' Plot heatmap for module classification
#'
#' Pure plotting function. It expects the pre-computed heatmap inputs
#' from `order_and_prepare_heatmap()`, and draws the figure only.
#'
#' @param hp list returned by order_and_prepare_heatmap()
#' @param show_colnames logical, whether to show column names
#' @param scale one of "none", "row", "column" passed to pheatmap (default: "column")
#' @param border_color cell border color (default: "white")
#' @param cellwidth cell width in pixels (default: 6)
#' @param cluster_cols logical, whether to cluster columns (default: FALSE)
#' @param cluster_rows logical, whether to cluster rows (default: FALSE)
#' @param legend logical, whether to show legend (default: TRUE)
#' @param ... additional args passed to pheatmap::pheatmap()
#' @return the pheatmap object (invisibly)
#' @import pheatmap
#' @export
plot_module_heatmap <- function(
    hp,
    show_colnames = TRUE,
    scale = "column",
    border_color = "white",
    cellwidth = 6,
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    legend = TRUE,
    ...
) {
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("Package 'pheatmap' is required. Please install it.")
  }
  ph <- pheatmap::pheatmap(
    hp$mat,
    cluster_cols = cluster_cols,
    cluster_rows = cluster_rows,
    scale = scale,
    show_colnames = show_colnames,
    legend = legend,
    annotation_col = hp$anno_col,
    annotation_colors = hp$ann_colors,
    border_color = border_color,
    cellwidth = cellwidth,
    ...
  )
  invisible(ph)
}


