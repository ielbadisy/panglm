#' Generic maximum likelihood helpers
#'
#' A likelihood is described by a contribution function
#' `contrib(theta, score)` returning a list with `ll` (one log-likelihood
#' contribution per independent unit: observation for pooled models,
#' individual for random-effects models) and, when `score = TRUE`, `score`
#' (one row of derivatives per unit). `w` holds one weight per unit.
#'
#' @keywords internal
#' @noRd
NULL

#' Unit-level scores by central differences of the contributions
#'
#' @keywords internal
#' @noRd
ml_numeric_scores <- function(contrib, theta) {
  base <- contrib(theta, FALSE)$ll
  out <- matrix(0, length(base), length(theta))
  for (j in seq_along(theta)) {
    h <- 1e-5 * max(1, abs(theta[j]))
    up <- theta; up[j] <- up[j] + h
    dn <- theta; dn[j] <- dn[j] - h
    out[, j] <- (contrib(up, FALSE)$ll - contrib(dn, FALSE)$ll) / (2 * h)
  }
  out
}

ml_scores <- function(contrib, theta) {
  s <- contrib(theta, TRUE)$score
  if (is.null(s)) s <- ml_numeric_scores(contrib, theta)
  s
}

#' Hessian of the weighted log-likelihood (Jacobian of the gradient)
#'
#' @keywords internal
#' @noRd
ml_hessian <- function(contrib, theta, w) {
  p <- length(theta)
  H <- matrix(0, p, p)
  for (j in seq_len(p)) {
    h <- 1e-4 * max(1, abs(theta[j]))
    up <- theta; up[j] <- up[j] + h
    dn <- theta; dn[j] <- dn[j] - h
    H[, j] <- (colSums(w * ml_scores(contrib, up)) - colSums(w * ml_scores(contrib, dn))) / (2 * h)
  }
  (H + t(H)) / 2
}

#' Maximize sum(w * ll) by BFGS followed by Newton polishing
#'
#' @keywords internal
#' @noRd
ml_optimize <- function(theta0, contrib, w, maxit = 200, tol = 1e-10) {
  cache <- new.env(parent = emptyenv())
  eval_at <- function(theta) {
    if (is.null(cache$theta) || !identical(cache$theta, theta)) {
      cache$theta <- theta
      cache$res <- tryCatch(contrib(theta, TRUE), error = function(e) NULL)
      if (!is.null(cache$res) && is.null(cache$res$score)) {
        cache$res$score <- ml_numeric_scores(contrib, theta)
      }
    }
    cache$res
  }
  objective <- function(theta) {
    r <- eval_at(theta)
    value <- if (is.null(r)) NA_real_ else -sum(w * r$ll)
    if (is.finite(value)) value else 1e100
  }
  gradient <- function(theta) {
    r <- eval_at(theta)
    if (is.null(r)) return(rep(0, length(theta)))
    g <- -colSums(w * r$score)
    g[!is.finite(g)] <- 0
    g
  }
  opt <- stats::optim(theta0, objective, gradient, method = "BFGS",
                      control = list(maxit = max(200L, maxit), reltol = min(tol, 1e-12)))
  theta <- opt$par
  value <- objective(theta)
  for (step in seq_len(20L)) {
    g <- -gradient(theta)
    if (max(abs(g)) < 1e-8 * max(1, abs(value))) break
    H <- ml_hessian(contrib, theta, w)
    delta <- tryCatch(solve(-H, g), error = function(e) NULL)
    if (is.null(delta) || any(!is.finite(delta))) break
    accepted <- FALSE
    for (halving in 0:20) {
      candidate <- theta + delta / 2^halving
      cv <- objective(candidate)
      if (cv <= value + 1e-12 * abs(value)) {
        theta <- candidate; value <- cv; accepted <- TRUE
        break
      }
    }
    if (!accepted) break
  }
  H <- ml_hessian(contrib, theta, w)
  vcov <- tryCatch(solve(-H), error = function(e) {
    matrix(NA_real_, length(theta), length(theta))
  })
  list(theta = theta, loglik = -value, hessian = H, vcov = vcov,
       score = ml_scores(contrib, theta), iterations = opt$counts[[1]],
       convergence = opt$convergence)
}

#' Sandwich covariance from unit scores and a bread matrix
#'
#' `cluster` holds one cluster label per unit; `NULL` gives HC1 (each unit
#' its own cluster, factor n / (n - p)). Units with zero weight carry no
#' information and are not counted. `clusters_only = TRUE` uses the factor
#' G / (G - 1) alone, the convention for random-effects likelihoods whose
#' units are already individuals.
#'
#' @keywords internal
#' @noRd
ml_sandwich <- function(bread, score, w, cluster = NULL, n_obs = NULL, clusters_only = FALSE) {
  keep <- w > 0
  score <- (w * score)[keep, , drop = FALSE]
  p <- ncol(score)
  if (is.null(cluster)) {
    n <- nrow(score)
    meat <- crossprod(score) * n / (n - p)
  } else {
    cluster <- cluster[keep]
    summed <- rowsum(score, cluster)
    G <- nrow(summed)
    if (G <= 1L) stop("cluster covariance requires at least two clusters", call. = FALSE)
    n <- n_obs %||% nrow(score)
    meat <- crossprod(summed) * (G / (G - 1))
    if (!clusters_only) meat <- meat * ((n - 1) / (n - p))
  }
  bread %*% meat %*% bread
}

#' Check that weights are constant within individuals and collapse them
#'
#' @keywords internal
#' @noRd
individual_weights <- function(weights, group_size) {
  group <- rep(seq_along(group_size), group_size)
  w <- as.numeric(rowsum(weights, group)) / group_size
  if (any(abs(weights - w[group]) > 1e-10 * pmax(1, abs(weights)))) {
    stop("with model = 'random' the weights must be constant within individuals ",
         "(one weight per individual, e.g. a baseline sampling weight)", call. = FALSE)
  }
  w
}

#' Map a cluster vector (one value per observation) to individuals
#'
#' @keywords internal
#' @noRd
individual_clusters <- function(cluster, group_size) {
  group <- rep(seq_along(group_size), group_size)
  per_group <- split(cluster, group)
  out <- vapply(per_group, function(v) {
    u <- unique(as.character(v))
    if (length(u) == 1L) u else NA_character_
  }, character(1))
  if (anyNA(out)) {
    stop("for a random-effects model each cluster must contain whole individuals",
         call. = FALSE)
  }
  out
}
