#' Compute module scores and classify samples
#'
#' Calculates module scores for each sample, assigns the dominant
#' transcriptional subtype, identifies Mixed samples based on relative
#' dominance, and calculates a gap-based SimplicityScore for sample ordering.
#'
#' @param expr Numeric matrix (genes x samples), assumed already normalized
#'   and row-wise z-scored.
#' @param gene_list Named list of gene vectors defining transcriptional modules.
#' @param module_order Character vector specifying the preferred module order.
#' @param mixed_threshold Numeric threshold for defining Mixed samples based on
#'   RelativeDominance. Default is 0.1.
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{scores}: data.frame containing module scores, Top1, Top2,
#'     MinScore, RelativeDominance, Mixed, TopCluster, and SimplicityScore.
#'     \item \code{modules_final}: final module order.
#'   }
#'
#' @export
compute_module_scores <- function(
    expr,
    gene_list,
    module_order = c(
      "Proliferating-like",
      "Osteoblast-like",
      "Chondroblast-like",
      "Fibroblast-like"
    ),
    mixed_threshold = 0.1
) {

  if (!is.matrix(expr) && !is.data.frame(expr)) {
    stop("'expr' must be a numeric matrix or data.frame.")
  }

  expr <- as.matrix(expr)

  if (!is.numeric(expr)) {
    stop("'expr' must contain numeric values.")
  }

  if (is.null(rownames(expr)) || is.null(colnames(expr))) {
    stop("'expr' must have gene names as row names and sample names as column names.")
  }

  if (!is.list(gene_list) || is.null(names(gene_list))) {
    stop("'gene_list' must be a named list of gene vectors.")
  }

  modules <- names(gene_list)

  if (length(modules) < 2) {
    stop("At least two gene modules are required.")
  }

  if (!is.numeric(mixed_threshold) ||
      length(mixed_threshold) != 1 ||
      is.na(mixed_threshold) ||
      mixed_threshold < 0 ||
      mixed_threshold > 1) {
    stop("'mixed_threshold' must be a numeric value between 0 and 1.")
  }

  scores <- data.frame(row.names = colnames(expr))

  # Module scores
  for (m in modules) {
    genes <- intersect(gene_list[[m]], rownames(expr))

    if (length(genes) > 0) {
      scores[[m]] <- colMeans(
        expr[genes, , drop = FALSE],
        na.rm = TRUE
      )
    } else {
      warning(
        sprintf(
          "No genes from module '%s' were found in the expression matrix.",
          m
        )
      )
      scores[[m]] <- NA_real_
    }
  }

  modules_final <- c(
    intersect(module_order, modules),
    setdiff(modules, module_order)
  )

  score_mat <- as.matrix(
    scores[, modules_final, drop = FALSE]
  )

  n_valid <- rowSums(is.finite(score_mat))

  # Top1, Top2, and minimum score
  top2 <- t(
    apply(
      score_mat,
      1,
      function(x) {
        x <- x[is.finite(x)]

        if (length(x) == 0) {
          return(c(NA_real_, NA_real_))
        }

        x <- sort(x, decreasing = TRUE)

        if (length(x) == 1) {
          return(c(x[1], NA_real_))
        }

        x[1:2]
      }
    )
  )

  scores$Top1 <- top2[, 1]
  scores$Top2 <- top2[, 2]

  scores$MinScore <- apply(
    score_mat,
    1,
    function(x) {
      x <- x[is.finite(x)]

      if (length(x) == 0) {
        NA_real_
      } else {
        min(x)
      }
    }
  )

  # Relative dominance
  scores$RelativeDominance <- with(
    scores,
    (Top1 - Top2) / (Top1 - MinScore)
  )

  same_score <- (
    n_valid >= 2 &
      is.finite(scores$Top1) &
      is.finite(scores$Top2) &
      is.finite(scores$MinScore) &
      scores$Top1 == scores$MinScore
  )

  scores$RelativeDominance[same_score] <- 0
  scores$RelativeDominance[n_valid < 2] <- NA_real_

  # Mixed classification
  scores$Mixed <- ifelse(
    is.na(scores$RelativeDominance),
    NA,
    scores$RelativeDominance <= mixed_threshold
  )

  # Dominant subtype
  scores$TopCluster <- apply(
    score_mat,
    1,
    function(x) {
      if (all(!is.finite(x))) {
        return(NA_character_)
      }

      x[!is.finite(x)] <- -Inf

      modules_final[which.max(x)]
    }
  )

  scores$TopCluster[
    !is.na(scores$Mixed) & scores$Mixed
  ] <- "Mixed"

  scores$TopCluster <- factor(
    scores$TopCluster,
    levels = c(modules_final, "Mixed")
  )

  # Gap-based SimplicityScore, used only for ordering
  calc_simplicity <- function(x) {
    x <- x[is.finite(x)]

    if (length(x) == 0) {
      return(NA_real_)
    }

    r <- sort(
      as.numeric(x),
      decreasing = TRUE
    )

    N <- length(r)

    if (N == 1) {
      return(0)
    }

    if (N == 2) {
      return(max(r) - min(r))
    }

    adds <- sum(r[1] - r[2:N])
    m <- N - 2

    if (m > 0) {
      M1 <- matrix(
        rep(r[2:(N - 1)], each = m),
        nrow = m
      )

      M2 <- matrix(
        rep(r[3:N], times = m),
        nrow = m
      )

      adns <- sum(
        (M1 - M2)[
          upper.tri(
            matrix(0, m, m),
            diag = TRUE
          )
        ]
      )
    } else {
      adns <- 0
    }

    (adds - adns) *
      (r[1] - r[N]) /
      (N - 1)
  }

  scores$SimplicityScore <- apply(
    score_mat,
    1,
    calc_simplicity
  )

  list(
    scores = scores,
    modules_final = modules_final
  )
}


#' Order samples and prepare heatmap inputs
#'
#' Samples are grouped by transcriptional subtype and ordered by decreasing
#' SimplicityScore within each subtype. SimplicityScore is used only for
#' ordering and is not displayed in the heatmap.
#'
#' @param scores Data.frame returned in the \code{scores} element of
#'   \code{compute_module_scores()}.
#' @param modules_final Character vector specifying the module order.
#'
#' @return A list containing the heatmap matrix, subtype annotation,
#'   annotation colors, sample order, and ordered score table.
#'
#' @export
order_and_prepare_heatmap <- function(scores, modules_final) {

  subtype_levels <- c(
    modules_final,
    "Mixed"
  )

  # Order by subtype and then SimplicityScore
  ord <- order(
    factor(
      as.character(scores$TopCluster),
      levels = subtype_levels
    ),
    -scores$SimplicityScore,
    na.last = TRUE
  )

  s_ord <- scores[
    ord,
    ,
    drop = FALSE
  ]

  mat <- t(
    s_ord[
      ,
      modules_final,
      drop = FALSE
    ]
  )

  # Display subtype only
  anno_col <- data.frame(
    Subtype = factor(
      as.character(s_ord$TopCluster),
      levels = subtype_levels
    ),
    row.names = colnames(mat),
    check.names = FALSE
  )

  cluster_colors <- c(
    "Proliferating-like" = "#FFFFB3",
    "Osteoblast-like"    = "#8DD3C7",
    "Chondroblast-like"  = "#FB8072",
    "Fibroblast-like"    = "#BEBADA",
    "Mixed"              = "gray60"
  )

  # Generate colors for custom module names
  missing_modules <- setdiff(
    subtype_levels,
    names(cluster_colors)
  )

  if (length(missing_modules) > 0) {
    extra_colors <- grDevices::rainbow(
      length(missing_modules)
    )

    names(extra_colors) <- missing_modules

    cluster_colors <- c(
      cluster_colors,
      extra_colors
    )
  }

  ann_colors <- list(
    Subtype = cluster_colors[subtype_levels]
  )

  list(
    mat = mat,
    anno_col = anno_col,
    ann_colors = ann_colors,
    order = ord,
    scores_ordered = s_ord
  )
}


#' Plot module classification heatmap
#'
#' Draws a heatmap using the prepared inputs from
#' \code{order_and_prepare_heatmap()}.
#'
#' @param hp List returned by \code{order_and_prepare_heatmap()}.
#' @param show_colnames Logical; whether to show sample names.
#' @param scale Character; scaling method passed to \code{pheatmap}.
#'   One of \code{"none"}, \code{"row"}, or \code{"column"}.
#' @param border_color Cell border color.
#' @param cellwidth Cell width in pixels.
#' @param cluster_cols Logical; whether to cluster columns.
#' @param cluster_rows Logical; whether to cluster rows.
#' @param legend Logical; whether to show the heatmap legend.
#' @param ... Additional arguments passed to \code{pheatmap::pheatmap()}.
#'
#' @return The pheatmap object, returned invisibly.
#'
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
