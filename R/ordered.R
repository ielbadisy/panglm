#' Ordered outcome models (internal engine)
#'
#' `P(y <= j | x, u) = F(c_j - x'beta - u)`, `j = 1, ..., J - 1`, with `F`
#' the logistic or standard normal distribution function. The cutpoints are
#' estimated on an unconstrained scale (first cutpoint, then log
#' increments), so they stay ordered.
#'
#' @keywords internal
#' @noRd
NULL

ordered_link <- function(link) {
  switch(link,
    logit = list(F = stats::plogis, f = stats::dlogis, q = stats::qlogis, var = pi^2 / 3),
    probit = list(F = stats::pnorm, f = stats::dnorm, q = stats::qnorm, var = 1)
  )
}

ordered_cuts <- function(par) cumsum(c(par[1], exp(par[-1])))

ordered_cut_jacobian <- function(par) {
  m <- length(par)
  M <- matrix(0, m, m)
  M[, 1] <- 1
  if (m > 1) for (j in 2:m) M[j:m, j] <- exp(par[j])
  M
}

ordered_response <- function(y) {
  if (is.factor(y)) {
    y <- droplevels(y)
    return(list(code = as.integer(y), levels = levels(y)))
  }
  if (is.logical(y)) y <- as.numeric(y)
  if (!is.numeric(y) || any(!is.finite(y)) || any(abs(y - round(y)) > 1e-8)) {
    stop("an ordered response must be a factor, an ordered factor or integer-valued",
         call. = FALSE)
  }
  values <- sort(unique(y))
  list(code = match(y, values), levels = as.character(values))
}

ordered_pieces <- function(eta, y, cuts, lf, grad) {
  J <- length(cuts) + 1L
  a <- c(cuts, Inf)[y] - eta
  b <- c(-Inf, cuts)[y] - eta
  P <- ifelse(b > 0, lf$F(-b) - lf$F(-a), lf$F(a) - lf$F(b))
  P <- pmax(P, 1e-300)
  out <- list(logP = log(P))
  if (grad) {
    fa <- lf$f(a); fa[!is.finite(a)] <- 0
    fb <- lf$f(b); fb[!is.finite(b)] <- 0
    out$deta <- -(fa - fb) / P
    Dc <- matrix(0, length(y), J - 1L)
    i1 <- which(y <= J - 1L)
    Dc[cbind(i1, y[i1])] <- fa[i1] / P[i1]
    i2 <- which(y >= 2L)
    Dc[cbind(i2, y[i2] - 1L)] <- Dc[cbind(i2, y[i2] - 1L)] - fb[i2] / P[i2]
    out$dcut <- Dc
  }
  out
}

ordered_probabilities <- function(eta, cuts, lf) {
  Fm <- vapply(cuts, function(cj) lf$F(cj - eta), numeric(length(eta)))
  if (!is.matrix(Fm)) Fm <- matrix(Fm, nrow = length(eta))
  Fm <- cbind(0, Fm, 1)
  Fm[, -1, drop = FALSE] - Fm[, -ncol(Fm), drop = FALSE]
}

ordered_contrib_pooled <- function(X, y, J, lf) {
  k <- ncol(X)
  function(theta, score) {
    cp <- theta[k + seq_len(J - 1L)]
    pc <- ordered_pieces(as.numeric(X %*% theta[seq_len(k)]), y, ordered_cuts(cp), lf, score)
    out <- list(ll = pc$logP)
    if (score) out$score <- cbind(pc$deta * X, pc$dcut %*% ordered_cut_jacobian(cp))
    out
  }
}

ordered_contrib_random <- function(X, y, J, lf, group, gh) {
  k <- ncol(X)
  G <- max(group)
  R <- length(gh$nodes)
  log_w <- log(gh$weights)
  function(theta, score) {
    cp <- theta[k + seq_len(J - 1L)]
    cuts <- ordered_cuts(cp)
    sigma <- exp(theta[k + J])
    eta <- as.numeric(X %*% theta[seq_len(k)])
    shifts <- sqrt(2) * sigma * gh$nodes
    L <- matrix(0, G, R)
    for (q in seq_len(R)) {
      L[, q] <- log_w[q] + rowsum(ordered_pieces(eta + shifts[q], y, cuts, lf, FALSE)$logP,
                                  group)[, 1]
    }
    top <- apply(L, 1, max)
    ll <- top + log(rowSums(exp(L - top)))
    out <- list(ll = ll)
    if (score) {
      post <- exp(L - ll)
      geta <- 0; gcut <- 0; gsig <- 0
      for (q in seq_len(R)) {
        pc <- ordered_pieces(eta + shifts[q], y, cuts, lf, TRUE)
        pq <- post[group, q]
        geta <- geta + pq * pc$deta
        gcut <- gcut + pq * pc$dcut
        gsig <- gsig + pq * pc$deta * shifts[q]
      }
      out$score <- rowsum(cbind(geta * X, gcut %*% ordered_cut_jacobian(cp), gsig), group)
      dimnames(out$score) <- NULL
    }
    out
  }
}

#' Fit an ordered model; called by panglm()
#'
#' @keywords internal
#' @noRd
fit_ordered <- function(X, y, levels, link, model, group_start, group_size, weights,
                        R, maxit, tol) {
  lf <- ordered_link(link)
  J <- length(levels)
  if (J < 2L) stop("an ordered response needs at least two observed categories", call. = FALSE)
  k <- ncol(X)
  group <- rep(seq_along(group_size), group_size)
  cut_names <- paste0(levels[-J], "|", levels[-1])

  if (model == "within") return(fit_ordered_within(X, y, J, group, maxit, tol))

  w_obs <- weights %||% rep(1, length(y))
  props <- cumsum(as.numeric(tapply(w_obs, factor(y, levels = seq_len(J)), sum)))
  props <- pmin(pmax(props / props[J], 1e-4), 1 - 1e-4)[-J]
  cuts0 <- lf$q(props)
  cuts0 <- cummax(cuts0 + seq_along(cuts0) * 1e-6)
  par0 <- c(cuts0[1], log(pmax(diff(cuts0), 1e-3)))

  pooled <- ml_optimize(c(rep(0, k), par0), ordered_contrib_pooled(X, y, J, lf), w_obs,
                        maxit, tol)
  if (model == "pooling") {
    res <- pooled
    unit_weights <- w_obs
    sigma <- NULL
  } else {
    gh <- gauss_hermite_quadrature(R)
    unit_weights <- if (is.null(weights)) rep(1, length(group_size)) else
      individual_weights(weights, group_size)
    inflate <- sqrt(1 + 1 / lf$var)
    p_beta <- pooled$theta[seq_len(k)] * inflate
    p_cut <- pooled$theta[k + seq_len(J - 1L)]
    p_cut <- c(p_cut[1] * inflate, p_cut[-1] + log(inflate))
    res <- ml_optimize(c(p_beta, p_cut, 0), ordered_contrib_random(X, y, J, lf, group, gh),
                       unit_weights, maxit, tol)
    sigma <- exp(res$theta[k + J])
  }

  cp <- res$theta[k + seq_len(J - 1L)]
  p <- length(res$theta)
  G_map <- diag(1, p)
  G_map[k + seq_len(J - 1L), k + seq_len(J - 1L)] <- ordered_cut_jacobian(cp)
  if (!is.null(sigma)) G_map[p, p] <- sigma
  beta <- res$theta[seq_len(k)]
  cuts <- ordered_cuts(cp)
  coefs <- c(beta, cuts)
  names(coefs) <- c(colnames(X), cut_names)
  V_nat <- G_map %*% res$vcov %*% t(G_map)
  keep <- seq_len(k + J - 1L)
  vcov <- V_nat[keep, keep, drop = FALSE]
  dimnames(vcov) <- list(names(coefs), names(coefs))

  eta <- as.numeric(X %*% beta)
  probs <- if (is.null(sigma)) ordered_probabilities(eta, cuts, lf) else
    ordered_marginal_probabilities(eta, cuts, lf, sigma, gauss_hermite_quadrature(R))
  colnames(probs) <- levels

  n_pos <- sum(w_obs > 0)
  list(coefficients = coefs, vcov = vcov, beta_names = colnames(X), cutpoints = cuts,
       sigma = sigma, sigma_se = if (!is.null(sigma)) sqrt(V_nat[p, p]),
       loglik = res$loglik, probabilities = probs, linear.predictors = eta,
       fitted.values = NULL, ordered_levels = levels, link = link,
       theta = res$theta, theta_vcov = res$vcov, theta_score = res$score,
       theta_map = G_map, unit_weights = unit_weights, quadrature_nodes = if (!is.null(sigma)) R,
       npar = p, df.residual = n_pos - p, iterations = res$iterations,
       convergence = res$convergence)
}

ordered_marginal_probabilities <- function(eta, cuts, lf, sigma, gh) {
  out <- 0
  for (q in seq_along(gh$nodes)) {
    out <- out + gh$weights[q] * ordered_probabilities(eta + sqrt(2) * sigma * gh$nodes[q], cuts, lf)
  }
  out
}

#' Fixed-effects ordered logit by "blow-up and cluster"
#'
#' Baetschmann, Staub and Winkelmann (2015): every observation is copied
#' once per cutoff `k = 2, ..., J`, the outcome is dichotomized as
#' `1(y >= k)`, and a conditional (Chamberlain) logit is fitted to the
#' stacked data with one stratum per individual and cutoff. The copies of
#' one individual are dependent, so standard errors are clustered by
#' individual.
#'
#' @keywords internal
#' @noRd
fit_ordered_within <- function(X, y, J, group, maxit, tol) {
  idx <- integer(0); sgroup <- integer(0); dy <- numeric(0)
  n_g <- tabulate(group)
  G <- length(n_g)
  for (cut in 2:J) {
    d <- as.numeric(y >= cut)
    s <- as.numeric(rowsum(d, group))
    informative <- (s > 0 & s < n_g)[group]
    rows <- which(informative)
    idx <- c(idx, rows)
    sgroup <- c(sgroup, (cut - 2L) * G + group[rows])
    dy <- c(dy, d[rows])
  }
  if (!length(idx)) {
    stop("no individual changes category over time; the fixed-effects ordered logit ",
         "is not identified", call. = FALSE)
  }
  Xs <- X[idx, , drop = FALSE]
  sizes <- rle(sgroup)$lengths
  starts <- as.integer(c(0L, cumsum(sizes)[-length(sizes)]))
  fit <- fit_within_binomial(Xs, dy, starts, as.integer(sizes), maxit, tol)
  info <- conditional_logit_loglik_grad_cpp(fit$coefficients, Xs, dy, starts, as.integer(sizes))
  score <- as.numeric(info$score_eta) * Xs
  bread <- fit$bread
  k <- ncol(X)
  out <- list(coefficients = fit$coefficients, vcov = fit$vcov, beta_names = colnames(X),
              loglik = fit$loglik, linear.predictors = as.numeric(X %*% fit$coefficients),
              fitted.values = NULL, buc_index = idx, buc_score = score, buc_bread = bread,
              n_used_groups = length(unique(group[idx])), npar = k, df.residual = NA_real_,
              iterations = fit$iterations, convergence = fit$convergence)
  out
}

ordered_buc_vcov <- function(object, cluster_obs) {
  cl <- cluster_obs[object$buc_index]
  summed <- rowsum(object$buc_score, cl)
  G <- nrow(summed)
  if (G <= 1L) stop("cluster covariance requires at least two clusters", call. = FALSE)
  v <- object$buc_bread %*% (crossprod(summed) * G / (G - 1)) %*% object$buc_bread
  dimnames(v) <- list(names(object$coefficients), names(object$coefficients))
  v
}

#' @export
vcov.panordered <- function(object, type = NULL, cluster = NULL, ...) {
  if (is.null(type)) return(object$vcov)
  type <- match.arg(type, c("classical", "HC1", "cluster"))
  if (type == "classical") return(object$vcov_classical %||% object$vcov)
  cl_obs <- panglm_cluster_vector(object, cluster)
  if (object$model == "within") {
    if (type == "HC1") {
      stop("the fixed-effects ordered logit reports cluster-robust standard errors only",
           call. = FALSE)
    }
    return(ordered_buc_vcov(object, cl_obs))
  }
  if (object$model == "random") {
    if (type == "HC1") {
      stop("type = 'HC1' is not defined for random-effects models; use 'cluster'",
           call. = FALSE)
    }
    cl <- individual_clusters(cl_obs, object$group_size)
    v <- ml_sandwich(object$theta_vcov, object$theta_score, object$unit_weights, cl,
                     clusters_only = TRUE)
  } else {
    v <- ml_sandwich(object$theta_vcov, object$theta_score, object$unit_weights,
                     if (type == "cluster") cl_obs else NULL)
  }
  v <- object$theta_map %*% v %*% t(object$theta_map)
  keep <- seq_along(object$coefficients)
  v <- v[keep, keep, drop = FALSE]
  dimnames(v) <- list(names(object$coefficients), names(object$coefficients))
  v
}

ordered_shim <- function(object) {
  list(terms = object$terms, xlevels = object$xlevels, contrasts = object$contrasts,
       coefficients = object$coefficients[object$beta_names],
       cre_columns = object$cre_columns, index = object$index, formula = object$formula)
}

ordered_predictor <- function(object, category = NULL, cre_means = NULL) {
  lf <- ordered_link(object$link)
  k <- length(object$beta_names)
  J <- length(object$ordered_levels)
  gh <- if (object$model == "random") gauss_hermite_quadrature(object$quadrature_nodes)
  shim <- ordered_shim(object)
  function(theta, data) {
    X <- panglm_design(shim, data, cre_means)
    eta <- as.numeric(X %*% theta[seq_len(k)])
    cuts <- theta[k + seq_len(J - 1L)]
    probs <- if (is.null(gh)) ordered_probabilities(eta, cuts, lf) else
      ordered_marginal_probabilities(eta, cuts, lf, object$sigma, gh)
    if (is.null(category)) probs else probs[, category]
  }
}

#' Predict from an ordered panel model
#'
#' @param object a `"panordered"` fit (from [panglm()] with an ordered
#'   family)
#' @param newdata optional data frame; defaults to the fitted data
#' @param type `"prob"` for the matrix of category probabilities
#'   (population-averaged over the random intercept for `model = "random"`),
#'   `"class"` for the most likely category, or `"link"` for `x'beta`
#' @param ... unused
#' @return a matrix with one column per category (`"prob"`), a factor
#'   (`"class"`) or a numeric vector (`"link"`)
#' @family panglm methods
#' @export
predict.panordered <- function(object, newdata = NULL, type = c("prob", "class", "link"), ...) {
  type <- match.arg(type)
  in_sample <- is.null(newdata)
  if (object$model == "within" && type != "link") {
    stop("category probabilities are not identified for the fixed-effects ordered logit ",
         "(the individual effects are conditioned out)", call. = FALSE)
  }
  shim <- ordered_shim(object)
  if (in_sample) {
    X <- object$X[, object$beta_names, drop = FALSE]
  } else {
    X <- panglm_design(shim, newdata)
  }
  eta <- as.numeric(X %*% object$coefficients[object$beta_names])
  if (type == "link") {
    out <- eta
  } else {
    lf <- ordered_link(object$link)
    probs <- if (object$model == "random") {
      ordered_marginal_probabilities(eta, object$cutpoints, lf, object$sigma,
                                     gauss_hermite_quadrature(object$quadrature_nodes))
    } else {
      ordered_probabilities(eta, object$cutpoints, lf)
    }
    colnames(probs) <- object$ordered_levels
    out <- if (type == "prob") probs else
      factor(object$ordered_levels[max.col(probs, ties.method = "first")],
             levels = object$ordered_levels)
  }
  if (!in_sample) return(out)
  if (is.matrix(out)) out[object$inverse_order, , drop = FALSE] else out[object$inverse_order]
}

#' @export
fitted.panordered <- function(object, ...) predict.panordered(object, type = "prob")

#' @export
residuals.panordered <- function(object, ...) {
  stop("residuals are not defined for ordered outcome models", call. = FALSE)
}

#' @rdname partialeffects
#' @export
partialeffects.panordered <- function(object, variables = NULL, level = 0.95, ...) {
  if (object$model == "within") {
    stop("partial effects on category probabilities are not identified for the ",
         "fixed-effects ordered logit; use model = 'random' with cre = TRUE instead",
         call. = FALSE)
  }
  variables <- ape_variables(object$terms, object$data, object$index, variables)
  cre_columns <- object$cre_columns %||% character(0)
  cre_means <- if (length(cre_columns)) object$X[, cre_columns, drop = FALSE] else NULL
  rows <- lapply(seq_along(object$ordered_levels), function(j) {
    predictor <- ordered_predictor(object, j, cre_means)
    out <- ape_engine(predictor, object$coefficients, object$vcov, object$data,
                      variables, object$weights, level, object$xlevels)
    out <- as.data.frame(out)
    cbind(category = object$ordered_levels[j], out, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out <- out[order(match(out$term, variables), out$contrast, match(out$category, object$ordered_levels)), ]
  rownames(out) <- NULL
  class(out) <- c("panglm_ape", "data.frame")
  attr(out, "vcov_type") <- object$vcov_type
  attr(out, "part") <- "category probability"
  out
}
