#' Recompute the coefficient covariance matrix for a panglm fit
#'
#' `type = "classical"` returns the model-based covariance the object was
#' fit with. `type = "HC1"` (heteroskedasticity-robust / White) and
#' `type = "cluster"` (cluster-robust sandwich, clustered by the panel
#' individual unless `cluster` is given) are only available for
#' `model = "pooling"` and `model = "within"`; random-effects models
#' already integrate out the individual-level correlation via the
#' likelihood, so they report model-based (information-matrix) standard
#' errors, matching the convention used by `lme4`/`glmmTMB`.
#'
#' @param object a `"panglm"` fit
#' @param type one of `"classical"`, `"HC1"`, `"cluster"`
#' @param cluster optional vector giving the clustering variable (one value
#'   per original observation, i.e. `nrow(data)` long); defaults to the
#'   panel individual identifier used to fit the model
#' @param ... unused
#' @export
vcov.panglm <- function(object, type = c("classical", "HC1", "cluster"), cluster = NULL, ...) {
  type <- match.arg(type)
  if (type == "classical") return(object$vcov)

  if (!object$model %in% c("pooling", "within")) {
    stop("type = '", type, "' is only available for model = 'pooling' or 'within'; ",
         "random-effects models already integrate out the individual-level ",
         "correlation via the likelihood, so a sandwich correction on top of ",
         "that isn't a well-motivated correction here -- they report ",
         "model-based (information-matrix) standard errors, matching the ",
         "convention used by lme4/glmmTMB", call. = FALSE)
  }

  if (object$model == "within" && object$family$family == "negbin") {
    return(robust_vcov_within_negbin(object, type, cluster))
  }

  score <- panglm_score(object)
  bread <- object$bread

  if (type == "HC1") {
    n <- nrow(score); k <- ncol(score)
    meat <- crossprod(score)
    v <- bread %*% meat %*% bread
    v <- v * n / (n - k)
  } else {
    cl <- if (is.null(cluster)) object$cluster_id else cluster
    if (length(cl) != nrow(score)) stop("'cluster' must have one value per observation in the fitted data", call. = FALSE)
    score_sum <- rowsum(score, cl)
    G <- nrow(score_sum); n <- nrow(score); k <- ncol(score)
    meat <- crossprod(score_sum)
    v <- bread %*% meat %*% bread
    v <- v * (G / (G - 1)) * ((n - 1) / (n - k))
  }
  dimnames(v) <- dimnames(object$vcov)
  v
}

#' Per-observation score (estimating-function) contributions
#'
#' The building block for sandwich/cluster-robust vcov: for pooled GLMs this
#' is the usual `(y - mu) * dmu/deta / V(mu) * X` score; for gaussian
#' "within" it is the OLS score on demeaned data; for poisson "within" it is
#' the score of the conditional (concentrated) likelihood, evaluated at the
#' original (non-demeaned) `X`.
#'
#' @keywords internal
#' @noRd
panglm_score <- function(object) {
  X <- object$X
  y <- object$y
  beta <- object$coefficients

  if (object$model == "pooling") {
    eta <- as.numeric(X %*% beta)
    mu <- linkinv_r(eta, object$family$link_id)
    return((y - mu) * mu_eta_r(eta, object$family$link_id) / variance_r(mu, object$family$family_id) * X)
  }

  if (object$model == "within" && object$family$family == "gaussian" && object$effect == "twoways") {
    dm <- demean_twoway(X, y, object$group_start, object$group_size, object$time_id, tol = 1e-10)
    resid <- as.numeric(dm$y - dm$X %*% beta)
    return(resid * dm$X)
  }

  if (object$model == "within" && object$family$family == "gaussian") {
    dm <- within_demean_cpp(X, y, object$group_start, object$group_size)
    resid <- as.numeric(dm$y - dm$X %*% beta)
    return(resid * dm$X)
  }

  if (object$model == "within" && object$family$family == "poisson") {
    eta <- as.numeric(X %*% beta)
    lit <- exp(eta)
    group <- rep(seq_along(object$group_size), object$group_size)
    Li <- as.numeric(rowsum(lit, group))
    Yi <- as.numeric(rowsum(y, group))
    w_obs <- (Yi / Li)[group]
    gradi <- y - w_obs * lit
    return(gradi * X)
  }

  stop("robust/cluster vcov is not implemented for this model/family combination", call. = FALSE)
}

#' Sandwich vcov for the Allison-Waterman dummy-variable FE-NB2 estimator
#'
#' The shared covariate coefficients and the per-group dummy intercepts are
#' estimated jointly, so a proper sandwich correction needs the score and
#' bread of the *full* (covariates + dummies) parameter vector, not just the
#' covariate block -- otherwise the correction ignores the sampling
#' variability the dummy intercepts contribute. This computes the full
#' (k+G)-dimensional sandwich and returns only the top-left k x k block
#' (the covariates), which is the only part `panglm` reports coefficients
#' for.
#'
#' @keywords internal
#' @noRd
robust_vcov_within_negbin <- function(object, type, cluster) {
  if (is.null(object$X_aug)) {
    stop("robust/cluster vcov is unavailable: this fit predates X_aug being stored ",
         "(re-fit the model)", call. = FALSE)
  }
  Xa <- object$X_aug
  keep_rows <- object$keep_rows
  y <- object$y[keep_rows]
  k <- ncol(object$X)

  beta_full <- c(object$coefficients, object$individual_effects)
  eta <- as.numeric(Xa %*% beta_full)
  mu <- exp(eta)
  theta <- object$theta
  w <- theta / (theta + mu)
  score <- ((y - mu) * w) * Xa

  bread_full <- object$vcov_full_aug

  if (type == "HC1") {
    n <- nrow(score); kk <- ncol(score)
    meat <- crossprod(score)
    v_full <- bread_full %*% meat %*% bread_full
    v_full <- v_full * n / (n - kk)
  } else {
    cl <- if (is.null(cluster)) object$cluster_id else cluster
    if (length(cl) != length(keep_rows)) {
      stop("'cluster' must have one value per observation in the fitted data", call. = FALSE)
    }
    cl <- cl[keep_rows]
    score_sum <- rowsum(score, cl)
    G <- nrow(score_sum); n <- nrow(score); kk <- ncol(score)
    meat <- crossprod(score_sum)
    v_full <- bread_full %*% meat %*% bread_full
    v_full <- v_full * (G / (G - 1)) * ((n - 1) / (n - kk))
  }

  v <- v_full[seq_len(k), seq_len(k), drop = FALSE]
  dimnames(v) <- dimnames(object$vcov)
  v
}

linkinv_r <- function(eta, link_id) {
  switch(link_id + 1L,
    eta,
    exp(eta),
    1 / (1 + exp(-eta)),
    stats::pnorm(eta)
  )
}

mu_eta_r <- function(eta, link_id) {
  switch(link_id + 1L,
    rep(1, length(eta)),
    exp(eta),
    { p <- 1 / (1 + exp(-eta)); p * (1 - p) },
    stats::dnorm(eta)
  )
}

variance_r <- function(mu, family_id) {
  switch(family_id + 1L,
    rep(1, length(mu)),
    mu,
    mu * (1 - mu)
  )
}
