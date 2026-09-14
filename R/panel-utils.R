#' Build contiguous panel group boundaries from an id/time index
#'
#' Sorts `data` by (id, time) and returns the sorted data plus 0-based
#' group_start / group_size integer vectors describing where each panel
#' unit's rows live in the sorted data -- the layout the C++ backends
#' expect.
#'
#' @param data a data.frame
#' @param index length-2 character vector: c(id_col, time_col)
#' @keywords internal
#' @noRd
build_panel_index <- function(data, index) {
  df <- as.data.frame(data)
  id_col <- index[1]
  time_col <- if (length(index) > 1) index[2] else NULL

  if (!id_col %in% names(df)) {
    stop("index column '", id_col, "' not found in data", call. = FALSE)
  }

  ord <- if (!is.null(time_col) && time_col %in% names(df)) {
    order(df[[id_col]], df[[time_col]])
  } else {
    order(df[[id_col]])
  }
  df <- df[ord, , drop = FALSE]

  if (!is.null(time_col) && anyDuplicated(df[, c(id_col, time_col), drop = FALSE])) {
    stop("each individual-time index pair must identify a unique observation",
         call. = FALSE)
  }

  # base::rle() on the already-sorted id column groups contiguous runs
  # directly (no as.character() coercion needed for an atomic id column, so
  # none of the R-level overhead the earlier data.table::rleidv() swap
  # avoided applies here).
  group_size <- as.integer(rle(df[[id_col]])$lengths)
  group_start <- as.integer(c(0L, cumsum(group_size)[-length(group_size)]))
  group_values <- df[[id_col]][group_start + 1L]

  list(data = df, order = ord, group_start = group_start,
       group_size = group_size, group_id = group_values,
       n_groups = length(group_size))
}

#' 0-based dense rank of a vector (ties get the same rank, ranks are
#' consecutive integers with no gaps).
#'
#' @param x a vector
#' @keywords internal
#' @noRd
dense_rank0 <- function(x) {
  as.integer(match(x, sort(unique(x))) - 1L)
}

#' Gauss-Hermite quadrature nodes and weights (physicists' convention)
#'
#' Computed via the Golub-Welsch algorithm (eigendecomposition of the
#' Jacobi tridiagonal matrix for Hermite polynomials), so no dependency
#' on the `statmod` package is needed. Weights are normalized by 1/sqrt(pi)
#' so that `sum(weights * f(nodes))` approximates `E[f(u)]` for
#' `u ~ N(0,1)` after the usual `sqrt(2)*sigma*node` change of variable.
#'
#' @param n number of quadrature points
#' @keywords internal
#' @noRd
gauss_hermite_quadrature <- function(n) {
  if (n < 1) stop("n must be >= 1", call. = FALSE)
  i <- seq_len(n - 1)
  b <- sqrt(i / 2)
  J <- matrix(0, n, n)
  if (n > 1) {
    J[cbind(i, i + 1)] <- b
    J[cbind(i + 1, i)] <- b
  }
  e <- eigen(J, symmetric = TRUE)
  nodes <- e$values
  # raw Gauss-Hermite weight is sqrt(pi) * v[1,]^2; dividing by sqrt(pi)
  # here gives weights already normalized for a standard-normal expectation
  weights <- e$vectors[1, ]^2
  ord <- order(nodes)
  list(nodes = nodes[ord], weights = weights[ord])
}
