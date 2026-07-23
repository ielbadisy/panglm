#' Two-way (individual + time) demeaning via alternating projections
#'
#' Iteratively demeans by `id` then by `time` until both group means are
#' simultaneously ~0 (the standard Gauss-Seidel / method-of-alternating-
#' projections algorithm for exact two-way fixed-effects demeaning; see
#' Guimaraes & Portugal 2010). Exact for both balanced and unbalanced
#' panels.
#'
#' @keywords internal
#' @noRd
demean_twoway <- function(X, y, id, time, tol = 1e-10, maxit = 10000) {
  xnames <- colnames(X)
  cols <- c(".y", xnames)
  dt <- data.table::data.table(.id = id, .time = time, .y = y)
  dt[, (xnames) := as.data.frame(X)]

  for (iter in seq_len(maxit)) {
    before <- as.matrix(dt[, ..cols])
    dt[, (cols) := lapply(.SD, function(v) v - mean(v)), by = .id, .SDcols = cols]
    dt[, (cols) := lapply(.SD, function(v) v - mean(v)), by = .time, .SDcols = cols]
    after <- as.matrix(dt[, ..cols])
    if (max(abs(after - before)) < tol) break
  }

  list(y = dt$.y, X = as.matrix(dt[, ..xnames]), iterations = iter)
}

fit_within_twoways_gaussian <- function(X, y, id, time, maxit, tol) {
  if (is.null(time)) stop("effect = 'twoways' requires a two-element index = c(id, time)", call. = FALSE)

  dm <- demean_twoway(X, y, id, time, tol = tol)
  res <- irls_fit_cpp(dm$X, dm$y, 0L, 0L, maxit, tol)

  n <- nrow(X); k <- ncol(X)
  n_id <- length(unique(id)); n_time <- length(unique(time))
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
