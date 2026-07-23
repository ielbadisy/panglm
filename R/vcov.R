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
         "random-effects models report model-based standard errors", call. = FALSE)
  }
  if (object$effect == "twoways") {
    stop("robust/cluster vcov for effect = 'twoways' is not yet implemented", call. = FALSE)
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

  if (object$model == "within" && object$family$family == "gaussian") {
    dm <- within_demean_cpp(X, y, object$group_start, object$group_size)
    resid <- as.numeric(dm$y - dm$X %*% beta)
    return(resid * dm$X)
  }

  if (object$model == "within" && object$family$family %in% c("poisson", "negbin")) {
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
