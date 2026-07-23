#' Two-way (individual + time) demeaning via alternating projections
#'
#' Iteratively demeans by `id` then by `time` until both group means are
#' simultaneously ~0 (the standard Gauss-Seidel / method-of-alternating-
#' projections algorithm for exact two-way fixed-effects demeaning; see
#' Guimaraes & Portugal 2010). Exact for both balanced and unbalanced
#' panels. The group-sum reduction and the demeaning step are parallelized
#' across observations in C++ (RcppParallel), the same hot-path treatment
#' as the rest of the package's estimators.
#'
#' @keywords internal
#' @noRd
demean_twoway <- function(X, y, group_start, group_size, time, tol = 1e-10, maxit = 10000) {
  # id groups are already contiguous/sorted (from build_panel_index), so no
  # factor()/match() is needed to derive their codes; only the time
  # dimension needs a grouping pass, done via data.table's fast dense rank.
  n_id <- length(group_size)
  id_code <- rep.int(seq_len(n_id) - 1L, group_size)
  time_code <- as.integer(data.table::frankv(time, ties.method = "dense") - 1L)
  n_time <- max(time_code) + 1L

  M <- cbind(y, X)
  dm <- twoway_demean_cpp(M, id_code, time_code, n_id, n_time, maxit, tol)

  list(y = dm$M[, 1], X = dm$M[, -1, drop = FALSE], iterations = dm$iterations, n_time = n_time)
}

fit_within_twoways_gaussian <- function(X, y, group_start, group_size, time, maxit, tol) {
  if (is.null(time)) stop("effect = 'twoways' requires a two-element index = c(id, time)", call. = FALSE)

  dm <- demean_twoway(X, y, group_start, group_size, time, tol = tol)
  res <- irls_fit_cpp(dm$X, dm$y, 0L, 0L, maxit, tol)

  n <- nrow(X); k <- ncol(X)
  n_id <- length(group_size); n_time <- dm$n_time
  # two-way within df: n - k - (n_id - 1) - (n_time - 1) - 1, matching plm's
  # model = "within", effect = "twoways" convention
  df_resid <- n - k - (n_id - 1) - (n_time - 1) - 1

  resid <- dm$y - dm$X %*% res$coefficients
  sigma2 <- sum(resid^2) / max(1, df_resid)
  vcov <- res$vcov_unscaled * sigma2
  bread <- res$vcov_unscaled
  coefs <- as.numeric(res$coefficients)
  names(coefs) <- colnames(X)
  dimnames(vcov) <- dimnames(bread) <- list(colnames(X), colnames(X))

  list(coefficients = coefs, vcov = vcov, bread = bread, fitted.values = NULL,
       loglik = NA_real_, dispersion = sigma2, df.residual = df_resid,
       iterations = res$iterations, demean_iterations = dm$iterations)
}
