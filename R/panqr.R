#' Quantile regression for panel data with individual effects
#'
#' Fixed-effects quantile regression for linear panel models: the two-step
#' estimator of Canay (2011) and the method of moments quantile regression
#' (MM-QR) of Machado and Santos Silva (2019), the estimator of the Stata
#' command `xtqreg`. A pooled quantile regression that ignores the
#' individual effects is available as a baseline. Several quantiles can be
#' estimated in one call.
#'
#' @section Estimators:
#' * `"canay"`: the individual effects are estimated from the within (mean)
#'   regression, `alpha_i = ybar_i - xbar_i'b - c` with `c = ybar - xbar'b`
#'   the overall constant (the normalization of Stata's `xtreg, fe`), and
#'   the quantile regression of the adjusted outcome `y - alpha_i` on the
#'   regressors (with an intercept) is fitted by [quantreg::rq.fit()] with
#'   `method = "br"`. The intercept is the `tau`-quantile of `c + error`.
#'   The model assumes that the individual effects are pure location shifts,
#'   and the estimator is consistent when both `N` and `T` grow.
#' * `"mmqr"`: the location-scale model
#'   `y_it = alpha_i + x_it'b + (delta_i + x_it'g) U_it`. The location
#'   parameters `b` come from the within regression of `y` on `x`; with `R`
#'   its residuals, the scale parameters `g` come from the within regression
#'   of `2 R (I(R >= 0) - mean(I(R >= 0)))`, an estimate of `|R|` that is
#'   robust to the first step; `q(tau)` is the sample `tau`-quantile of the
#'   standardized residuals `R / sigma_it`, with `sigma_it` the fitted scale
#'   including the individual effects. The quantile coefficients are
#'   `b(tau) = b + g q(tau)`. The scale variables are the regressors, as in
#'   `xtqreg`. Only slopes are reported: the intercept is absorbed in the
#'   individual effects `alpha_i(tau) = alpha_i + delta_i q(tau)`.
#' * `"pooling"`: [quantreg::rq.fit()] of `y` on the regressors, ignoring
#'   the panel structure (except for the clustered covariance).
#'
#' Time-invariant regressors are not identified with individual effects and
#' are dropped (with a warning) for `"canay"` and `"mmqr"`, as are
#' individuals observed only once (with a message).
#'
#' @section Covariance:
#' `vcov = "bootstrap"` resamples whole individuals with replacement (a
#' cluster bootstrap) `R` times, repeats every step of the estimator and
#' returns the covariance of the replicates. `vcov = "analytic"` gives:
#' * `"mmqr"`: the asymptotic covariance of Machado and Santos Silva (2019),
#'   with the matrix layout of `xtqreg` (version 1.5): it combines the
#'   influence of the location, scale and `q(tau)` estimates, with the
#'   density of the standardized error at `q(tau)` estimated by the
#'   difference quotient `2h / (q(tau + h) - q(tau - h))` (as Stata's
#'   `qreg, vce(iid)`) and the Hall and Sheather bandwidth `h`. The scale
#'   influence uses the centered `V = 2 U (I(U >= 0) - mean(I(U >= 0))) - 1`
#'   of the paper (also used by the Stata command `mmqreg`); `xtqreg` 1.5
#'   omits the `- 1`, which inflates its standard errors away from the
#'   median (by about 30% at `tau = 0.9` in simulations, where the centered
#'   version matches the Monte Carlo standard deviation).
#' * `"canay"`: a cluster-robust (by individual) sandwich built from the
#'   influence function of the two-step estimator, which corrects the
#'   quantile regression score `(tau - I(e < 0)) x` for the estimation of
#'   the individual effects and of the within slopes. The conditional
#'   density is estimated with the Gaussian kernel and bandwidth of
#'   `quantreg::summary.rq(se = "ker")`. This is a large-`N`, large-`T`
#'   approximation; Canay (2011) recommends the bootstrap. It treats the
#'   individual effects as fixed, so its intercept SE omits the sampling
#'   variation of the mean individual effect that the cluster bootstrap
#'   includes; the slope SEs of the two agree.
#' * `"pooling"`: the same kernel sandwich without the first-step terms,
#'   clustered by individual.
#'
#' The cluster sandwiches carry no finite-sample factor. For `"mmqr"` the
#' location and scale coefficients are returned with the cluster-robust
#' within covariance of the rest of the package
#' (`G / (G - 1) * (n - 1) / (n - k)`, as `plm::vcovHC(type = "sss")`).
#'
#' @param formula a formula `y ~ regressors`
#' @param data a data.frame in long format
#' @param index length-1 or length-2 character vector giving the individual
#'   (and time) columns
#' @param tau numeric vector of quantiles, each strictly between 0 and 1
#' @param model `"canay"`, `"mmqr"` or `"pooling"`
#' @param vcov `"analytic"`, `"bootstrap"` or `"none"`; the default is
#'   `"bootstrap"` for `"canay"` and `"analytic"` otherwise
#' @param R number of bootstrap replications
#' @param na.action missing-data action passed to [stats::model.frame()]
#' @return an object of class `"panqr"`, a list with components
#'   `coefficients` (a matrix with one column per quantile), `vcov` (a list
#'   of covariance matrices, one per quantile, or `NULL` for
#'   `vcov = "none"`), `vcov_type`, `tau`, `model`, `boot` (the bootstrap
#'   replicates, an array `R x k x length(tau)`, for `vcov = "bootstrap"`),
#'   `alpha` (one per individual: for `"canay"` the normalized effects
#'   subtracted from the outcome, for `"mmqr"` the location effects including
#'   the constant), for `"mmqr"` also `location`, `scale` (lists with
#'   `coefficients` and `vcov`), `q` (the quantiles of the standardized
#'   errors) and `delta` (the individual scale effects including the
#'   constant, so that `alpha_i(tau) = alpha_i + delta_i q(tau)`),
#'   `dropped` (time-invariant regressors dropped), `n_singletons`
#'   (individuals observed once, dropped for `"canay"` and `"mmqr"`), `nobs`,
#'   `n_groups`,
#'   `call` and `formula`. It has `print`, `summary`, `coef`, `vcov`, `nobs`
#'   and `confint` methods.
#' @references
#' Canay, I. A. (2011). A simple approach to quantile regression for panel
#' data. *The Econometrics Journal*, 14(3), 368-386.
#' \doi{10.1111/j.1368-423X.2011.00349.x}
#'
#' Machado, J. A. F. and Santos Silva, J. M. C. (2019). Quantiles via
#' moments. *Journal of Econometrics*, 213(1), 145-173.
#' \doi{10.1016/j.jeconom.2019.04.009}
#' @examples
#' # Simulated panel of 300 patients over 6 years. Annual health costs have
#' # a patient effect and a heavy upper tail that widens with chronic
#' # conditions: a location-scale model in which the effect of chronic
#' # conditions grows across the cost distribution.
#' set.seed(1)
#' N <- 300; TT <- 6
#' d <- data.frame(id = rep(seq_len(N), each = TT), year = rep(seq_len(TT), N))
#' alpha <- rep(rnorm(N), each = TT)
#' d$chronic <- rpois(N * TT, 1 + 0.5 * (alpha > 0))
#' d$age <- rep(sample(40:70, N, TRUE), each = TT) + d$year
#' d$cost <- 2 + alpha + 0.5 * d$chronic + 0.02 * d$age +
#'   (1 + 0.3 * d$chronic) * rexp(N * TT)
#'
#' if (requireNamespace("quantreg", quietly = TRUE)) {
#'   # median and upper tail of costs, MM-QR with analytic standard errors
#'   mm <- panqr(cost ~ chronic + age, data = d, index = c("id", "year"),
#'               tau = c(0.5, 0.9), model = "mmqr")
#'   summary(mm)
#'
#'   # Canay two-step estimator with a small cluster bootstrap
#'   cn <- panqr(cost ~ chronic + age, data = d, index = c("id", "year"),
#'               tau = 0.9, model = "canay", R = 50)
#'   coef(cn)
#' }
#' @export
panqr <- function(formula, data, index, tau = 0.5,
                  model = c("canay", "mmqr", "pooling"),
                  vcov = NULL, R = 200, na.action = stats::na.omit) {
  model <- match.arg(model)
  if (is.null(vcov)) vcov <- if (model == "canay") "bootstrap" else "analytic"
  vcov <- match.arg(vcov, c("analytic", "bootstrap", "none"))
  if (!is.numeric(tau) || !length(tau) || anyNA(tau) || any(tau <= 0 | tau >= 1)) {
    stop("'tau' must be a numeric vector of quantiles strictly between 0 and 1", call. = FALSE)
  }
  tau <- sort(unique(tau))
  if (model %in% c("canay", "pooling") && !requireNamespace("quantreg", quietly = TRUE)) {
    stop("model = '", model, "' needs the 'quantreg' package", call. = FALSE)
  }
  if (vcov == "bootstrap" && (!is.numeric(R) || length(R) != 1L || R < 2)) {
    stop("'R' must be a number of bootstrap replications >= 2", call. = FALSE)
  }
  if (missing(index)) stop("'index' is required, e.g. index = c(\"id\", \"time\")", call. = FALSE)
  if (!is.character(index) || !length(index) %in% 1:2 || anyDuplicated(index)) {
    stop("'index' must contain one or two distinct column names", call. = FALSE)
  }
  missing_index <- setdiff(index, names(data))
  if (length(missing_index)) {
    stop("index column(s) not found in data: ", paste(missing_index, collapse = ", "),
         call. = FALSE)
  }
  if (!inherits(formula, "formula") || length(formula) != 3L) {
    stop("'formula' must be a two-sided formula 'y ~ regressors'", call. = FALSE)
  }

  mf <- stats::model.frame(formula, data = data, na.action = na.action)
  model_rows <- match(rownames(mf), rownames(data))
  if (anyNA(model_rows)) stop("could not align model-frame rows with data", call. = FALSE)
  y <- stats::model.response(mf)
  if (!is.numeric(y) || any(!is.finite(y))) {
    stop("the response must be finite and numeric", call. = FALSE)
  }
  X <- stats::model.matrix(formula, data = mf)
  if (!"(Intercept)" %in% colnames(X)) {
    stop("the formula must keep the intercept (it is absorbed by the individual effects)",
         call. = FALSE)
  }
  panel_cols <- data[model_rows, index, drop = FALSE]
  if (anyNA(panel_cols)) stop("panel index columns must not contain missing values", call. = FALSE)
  panel <- build_panel_index(cbind(panel_cols, .panqr_row = seq_len(nrow(mf))), index)
  ord <- panel$data$.panqr_row
  y <- as.numeric(y[ord])
  X <- X[ord, , drop = FALSE]
  group <- rep.int(seq_len(panel$n_groups), panel$group_size)

  dropped <- character(0)
  n_singletons <- 0L
  group_id <- panel$group_id
  if (model != "pooling") {
    # individuals observed once carry no within information (their within
    # residual is zero); they are dropped
    single <- panel$group_size == 1L
    n_singletons <- sum(single)
    if (n_singletons) {
      message(n_singletons, " individual(s) observed only once dropped")
      keep <- !single[group]
      y <- y[keep]
      X <- X[keep, , drop = FALSE]
      ord <- ord[keep]
      group_id <- group_id[!single]
      group <- rep.int(seq_along(group_id), panel$group_size[!single])
    }
    if (!length(group_id)) stop("no individual is observed more than once", call. = FALSE)
    if (max(panel$group_size) < if (model == "mmqr") 3L else 2L) {
      stop("model = '", model, "' needs individuals observed in several periods",
           " (T >= ", if (model == "mmqr") 3 else 2, ")", call. = FALSE)
    }
    slopes <- X[, colnames(X) != "(Intercept)", drop = FALSE]
    if (!ncol(slopes)) stop("no regressors besides the intercept", call. = FALSE)
    Xw <- qr_within(slopes, group)
    tv <- colSums(Xw^2) > 1e-10 * pmax(1, colSums(slopes^2))
    dropped <- colnames(slopes)[!tv]
    if (length(dropped)) {
      warning("time-invariant regressor(s) dropped: ", paste(dropped, collapse = ", "),
              call. = FALSE)
      X <- X[, setdiff(colnames(X), dropped), drop = FALSE]
    }
    if (ncol(X) < 2L) stop("no time-varying regressors left", call. = FALSE)
  }

  fit <- switch(model,
    canay = panqr_canay(y, X, group, tau, se = vcov == "analytic"),
    mmqr = panqr_mmqr(y, X[, -1L, drop = FALSE], group, tau, se = vcov == "analytic"),
    pooling = panqr_pooling(y, X, group, tau, se = vcov == "analytic")
  )

  if (vcov == "bootstrap") {
    fit$boot <- panqr_bootstrap(y, X, group, tau, model, R)
    fit$vcov <- lapply(seq_along(tau), function(j) {
      v <- stats::cov(fit$boot[, , j, drop = TRUE])
      if (!is.matrix(v)) v <- matrix(v, 1, 1)
      dimnames(v) <- list(rownames(fit$coefficients), rownames(fit$coefficients))
      v
    })
  }
  # keep an explicit NULL element so that fit$vcov does not partially match
  # fit$vcov_type
  if (vcov == "none") fit["vcov"] <- list(NULL)
  if (!is.null(fit$vcov)) names(fit$vcov) <- colnames(fit$coefficients)

  fit$call <- match.call()
  fit$formula <- formula
  fit$model <- model
  fit$tau <- tau
  fit$vcov_type <- vcov
  fit$R <- if (vcov == "bootstrap") R else NULL
  fit$dropped <- dropped
  fit$index <- index
  fit$nobs <- length(y)
  fit$n_groups <- length(group_id)
  fit$n_singletons <- n_singletons
  if (!is.null(fit$alpha)) names(fit$alpha) <- as.character(group_id)
  if (!is.null(fit$delta)) names(fit$delta) <- as.character(group_id)
  fit$panel_order <- ord
  fit$na.action <- attr(mf, "na.action")
  class(fit) <- "panqr"
  fit
}

# Within (individual-demeaned) transform of the columns of M.
qr_within <- function(M, group) {
  M <- as.matrix(M)
  M - (rowsum(M, group, reorder = TRUE) / as.vector(table(group)))[group, , drop = FALSE]
}

# Within regression of y on the slope matrix Xs: slopes, overall constant
# (ybar - xbar'b), individual effects normalized as in Stata's xtreg, fe
# (alpha_i = ybar_i - xbar_i'b - cons) and the within residuals.
qr_within_fit <- function(y, Xs, group) {
  Xw <- qr_within(Xs, group)
  yw <- as.vector(qr_within(y, group))
  qx <- qr(Xw)
  if (qx$rank < ncol(Xw)) stop("the within regressors are collinear", call. = FALSE)
  b <- qr.coef(qx, yw)
  names(b) <- colnames(Xs)
  size <- as.vector(table(group))
  ybar_i <- as.vector(rowsum(y, group, reorder = TRUE)) / size
  xbar_i <- rowsum(Xs, group, reorder = TRUE) / size
  cons <- mean(y) - sum(colMeans(Xs) * b)
  alpha <- ybar_i - as.vector(xbar_i %*% b) - cons
  list(coefficients = b, cons = cons, alpha = alpha, Xw = Xw,
       residuals = as.vector(yw - Xw %*% b), xbar_i = xbar_i)
}

# Cluster-robust covariance of a within regression, panglm convention.
qr_within_vcov <- function(Xw, e, group) {
  bread <- solve(crossprod(Xw))
  s <- rowsum(Xw * e, group)
  G <- nrow(s); n <- nrow(Xw); k <- ncol(Xw)
  v <- bread %*% crossprod(s) %*% bread * (G / (G - 1)) * ((n - 1) / (n - k))
  dimnames(v) <- list(colnames(Xw), colnames(Xw))
  v
}

qr_tau_names <- function(tau) paste0("tau=", as.character(tau))

# quantreg::rq.fit(method = "br") for each tau; a k x length(tau) matrix.
qr_rq_coef <- function(X, y, tau) {
  out <- vapply(tau, function(t) quantreg::rq.fit(X, y, tau = t, method = "br")$coefficients,
                numeric(ncol(X)))
  out <- matrix(out, nrow = ncol(X), dimnames = list(colnames(X), qr_tau_names(tau)))
  out
}

# Kernel density weights of quantreg::summary.rq(se = "ker").
qr_kernel_density <- function(e, tau) {
  n <- length(e)
  h <- qr_bandwidth_hs(tau, n)
  while ((tau - h < 0) || (tau + h > 1)) h <- h / 2
  h <- (stats::qnorm(tau + h) - stats::qnorm(tau - h)) *
    min(sqrt(stats::var(e)), (stats::quantile(e, 0.75) - stats::quantile(e, 0.25)) / 1.34)
  stats::dnorm(e / h) / h
}

# Hall and Sheather (1988) bandwidth, as quantreg::bandwidth.rq(hs = TRUE).
qr_bandwidth_hs <- function(tau, n, alpha = 0.05) {
  x0 <- stats::qnorm(tau)
  f0 <- stats::dnorm(x0)
  n^(-1 / 3) * stats::qnorm(1 - alpha / 2)^(2 / 3) * ((1.5 * f0^2) / (2 * x0^2 + 1))^(1 / 3)
}

# Cluster sandwich (sum f x x')^{-1} (sum_g s_g s_g') (sum f x x')^{-1}.
qr_sandwich <- function(X, f, score, group) {
  Hinv <- solve(crossprod(X * sqrt(f)))
  s <- rowsum(score, group)
  v <- Hinv %*% crossprod(s) %*% Hinv
  dimnames(v) <- list(colnames(X), colnames(X))
  v
}

panqr_pooling <- function(y, X, group, tau, se) {
  coef <- qr_rq_coef(X, y, tau)
  out <- list(coefficients = coef)
  if (se) {
    out$vcov <- lapply(seq_along(tau), function(j) {
      e <- as.vector(y - X %*% coef[, j])
      f <- qr_kernel_density(e, tau[j])
      qr_sandwich(X, f, X * (tau[j] - (e < 0)), group)
    })
  }
  out
}

panqr_canay <- function(y, X, group, tau, se) {
  Xs <- X[, -1L, drop = FALSE]
  w <- qr_within_fit(y, Xs, group)
  yadj <- y - w$alpha[group]
  coef <- qr_rq_coef(X, yadj, tau)
  out <- list(coefficients = coef, alpha = w$alpha, within = w$coefficients)
  if (se) {
    size <- as.vector(table(group))
    xbar_c <- sweep(w$xbar_i, 2, colMeans(Xs))[group, , drop = FALSE]
    Qinv <- solve(crossprod(w$Xw))
    u <- w$residuals
    out$vcov <- lapply(seq_along(tau), function(j) {
      e <- as.vector(yadj - X %*% coef[, j])
      f <- qr_kernel_density(e, tau[j])
      fX <- X * f
      m_i <- (rowsum(fX, group, reorder = TRUE) / size)[group, , drop = FALSE]
      m_c <- sweep(m_i, 2, colMeans(fX))
      H <- crossprod(fX, xbar_c)
      score <- X * (tau[j] - (e < 0)) - m_c * u + (w$Xw * u) %*% t(H %*% Qinv)
      qr_sandwich(X, f, score, group)
    })
  }
  out
}

# Sample tau-quantile, the smallest value with empirical cdf >= tau (the
# minimizer of the check function; it is unique unless n * tau is an integer,
# in which case the lower of the two minimizing order statistics is taken).
qr_sample_quantile <- function(x, tau) {
  x <- sort(x)
  x[max(1L, ceiling(length(x) * tau - 1e-10))]
}

panqr_mmqr <- function(y, Xs, group, tau, se) {
  loc <- qr_within_fit(y, Xs, group)
  u <- loc$residuals
  up <- (u >= 0) - mean(u >= 0)
  au <- 2 * u * up
  sc <- qr_within_fit(au, Xs, group)
  g <- sc$coefficients
  size <- as.vector(table(group))
  au_i <- as.vector(rowsum(au, group, reorder = TRUE)) / size
  sigma <- au_i[group] + as.vector(loc$Xw %*% g)
  if (any(sigma <= 0)) {
    warning(sprintf("%.1f%% of the fitted values of the scale function are not positive",
                    100 * mean(sigma <= 0)), call. = FALSE)
  }
  us <- u / sigma
  q <- vapply(tau, function(t) qr_sample_quantile(us, t), numeric(1))
  names(q) <- qr_tau_names(tau)
  b <- loc$coefficients
  coef <- outer(b, rep(1, length(tau))) + outer(g, q)
  dimnames(coef) <- list(names(b), qr_tau_names(tau))
  out <- list(coefficients = coef, q = q, alpha = loc$alpha + loc$cons,
              delta = sc$alpha + sc$cons,
              location = list(coefficients = b,
                              vcov = qr_within_vcov(loc$Xw, u, group)),
              scale = list(coefficients = g,
                           vcov = qr_within_vcov(sc$Xw, sc$residuals, group)))
  if (se) {
    # scale of individual i at the overall regressor means
    si <- au_i[group] - as.vector(sweep(loc$xbar_i, 2, colMeans(Xs))[group, , drop = FALSE] %*% g)
    out$vcov <- lapply(seq_along(tau), function(j) mmqr_vcov(loc$Xw, us, sigma, si, g, q[j], tau[j]))
  }
  out
}

# Asymptotic covariance of the MM-QR slopes, transcribed from xtqreg 1.5.
# Xw: within regressors; us: standardized residuals; sigma: fitted scale;
# si: scale evaluated at the overall regressor means; g: scale slopes.
# centered = TRUE uses V = av - 1 as in the paper; FALSE reproduces xtqreg.
mmqr_vcov <- function(Xw, us, sigma, si, g, q, tau, centered = TRUE) {
  n <- length(us)
  k <- ncol(Xw)
  Qinv <- solve(crossprod(Xw))
  Pxx <- crossprod(Xw * abs(sigma))
  Px <- colSums(Xw * sigma^2)
  vp <- (us >= 0) - mean(us >= 0)
  av <- 2 * us * vp - centered
  # density of the standardized error at q(tau): Hall-Sheather difference
  # quotient, as Stata's qreg, vce(iid) on a constant
  h <- qr_bandwidth_hs(tau, n)
  while ((tau - h <= 0) || (tau + h >= 1)) h <- h / 2
  fq <- 2 * h / (qr_sample_quantile(us, tau + h) - qr_sample_quantile(us, tau - h))
  w <- (tau - (us <= q)) / fq - us - q * av
  O <- matrix(0, 2 * k + 1, 2 * k + 1)
  i1 <- seq_len(k); i2 <- k + seq_len(k); i3 <- 2 * k + 1
  O[i1, i1] <- mean(us^2) * Pxx
  O[i1, i2] <- O[i2, i1] <- mean(us * av) * Pxx
  O[i2, i2] <- mean(av^2) * Pxx
  O[i3, i1] <- O[i1, i3] <- mean(us * w) * Px
  O[i3, i2] <- O[i2, i3] <- mean(av * w) * Px
  O[i3, i3] <- sum(w^2) * mean(si^2)
  XI <- cbind(Qinv, q * Qinv, g / sum(si))
  v <- XI %*% O %*% t(XI)
  dimnames(v) <- list(colnames(Xw), colnames(Xw))
  v
}

# Cluster bootstrap: resample individuals with replacement, refit.
panqr_bootstrap <- function(y, X, group, tau, model, R) {
  rows <- split(seq_along(y), group)
  G <- length(rows)
  k <- if (model == "mmqr") ncol(X) - 1L else ncol(X)
  out <- array(NA_real_, c(R, k, length(tau)))
  for (r in seq_len(R)) {
    pick <- sample.int(G, G, replace = TRUE)
    idx <- unlist(rows[pick], use.names = FALSE)
    g_new <- rep.int(seq_len(G), lengths(rows[pick]))
    est <- tryCatch(
      switch(model,
        canay = panqr_canay(y[idx], X[idx, , drop = FALSE], g_new, tau, FALSE)$coefficients,
        mmqr = suppressWarnings(panqr_mmqr(y[idx], X[idx, -1L, drop = FALSE], g_new,
                                           tau, FALSE)$coefficients),
        pooling = panqr_pooling(y[idx], X[idx, , drop = FALSE], g_new, tau, FALSE)$coefficients
      ),
      error = function(e) NULL)
    if (!is.null(est)) out[r, , ] <- est
  }
  ok <- stats::complete.cases(matrix(out, R))
  if (!all(ok)) {
    warning(sum(!ok), " bootstrap replication(s) failed and were discarded", call. = FALSE)
    out <- out[ok, , , drop = FALSE]
  }
  out
}

#' Covariance of a panqr fit
#'
#' @param object a `"panqr"` fit
#' @param tau `NULL` or one of the fitted quantiles
#' @param ... unused
#' @return the covariance matrix of the coefficients at quantile `tau`; when
#'   several quantiles were fitted and `tau` is `NULL`, a list of matrices,
#'   one per quantile
#' @examples
#' set.seed(2)
#' d <- data.frame(id = rep(1:100, each = 5), t = rep(1:5, 100))
#' d$x <- rnorm(500) + rep(rnorm(100), each = 5)
#' d$y <- d$x + rep(rnorm(100), each = 5) + (1 + 0.2 * abs(d$x)) * rnorm(500)
#' fit <- panqr(y ~ x, data = d, index = c("id", "t"), tau = c(0.25, 0.75),
#'              model = "mmqr")
#' vcov(fit, tau = 0.75)
#' @export
vcov.panqr <- function(object, tau = NULL, ...) {
  if (is.null(object$vcov)) stop("the fit has no covariance (vcov = 'none')", call. = FALSE)
  if (is.null(tau)) {
    if (length(object$tau) == 1L) return(object$vcov[[1]])
    return(object$vcov)
  }
  object$vcov[[qr_tau_position(object, tau)]]
}

qr_tau_position <- function(object, tau) {
  j <- which(abs(object$tau - tau) < 1e-10)
  if (length(tau) != 1L || !length(j)) {
    stop("'tau' must be one of the fitted quantiles: ", paste(object$tau, collapse = ", "),
         call. = FALSE)
  }
  j
}

#' @export
coef.panqr <- function(object, ...) {
  cf <- object$coefficients
  if (ncol(cf) == 1L) stats::setNames(cf[, 1], rownames(cf)) else cf
}

#' @export
nobs.panqr <- function(object, ...) object$nobs

#' @export
confint.panqr <- function(object, parm, level = 0.95, tau = NULL, ...) {
  if (missing(parm)) parm <- rownames(object$coefficients)
  one <- function(j) {
    cf <- object$coefficients[, j]
    se <- sqrt(diag(object$vcov[[j]]))[parm]
    crit <- stats::qnorm(1 - (1 - level) / 2)
    out <- cbind(cf[parm] - crit * se, cf[parm] + crit * se)
    colnames(out) <- sprintf("%.1f %%", c((1 - level) / 2, 1 - (1 - level) / 2) * 100)
    rownames(out) <- parm
    out
  }
  if (is.null(object$vcov)) stop("the fit has no covariance (vcov = 'none')", call. = FALSE)
  if (!is.null(tau)) return(one(qr_tau_position(object, tau)))
  if (length(object$tau) == 1L) return(one(1L))
  stats::setNames(lapply(seq_along(object$tau), one), colnames(object$coefficients))
}

qr_label <- function(x) {
  switch(x$model,
    canay = "Canay (2011) two-step fixed-effects quantile regression",
    mmqr = "MM-QR fixed-effects quantile regression (Machado and Santos Silva 2019)",
    pooling = "pooled quantile regression")
}

#' @export
print.panqr <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("panqr(", qr_label(x), ")\n\n", sep = "")
  print(x$coefficients, digits = digits)
  cat("\nN =", x$nobs, " groups =", x$n_groups, "\n")
  invisible(x)
}

#' Summarize a panqr fit
#'
#' @param object a `"panqr"` fit
#' @param ... unused
#' @return an object of class `"summary.panqr"`: a list with one
#'   coefficient table per quantile (`coefficients`: estimate, standard
#'   error, z value and p-value), the location and scale tables and the
#'   quantiles `q` of the standardized errors for MM-QR, and the fit
#'   information (`call`, `model`, `tau`, `nobs`, `n_groups`, `vcov_type`,
#'   `R`, `dropped`).
#' @examples
#' set.seed(2)
#' d <- data.frame(id = rep(1:100, each = 5), t = rep(1:5, 100))
#' d$x <- rnorm(500) + rep(rnorm(100), each = 5)
#' d$y <- d$x + rep(rnorm(100), each = 5) + (1 + 0.2 * abs(d$x)) * rnorm(500)
#' summary(panqr(y ~ x, data = d, index = c("id", "t"), tau = 0.9, model = "mmqr"))
#' @export
summary.panqr <- function(object, ...) {
  tab <- function(cf, v) {
    se <- if (is.null(v)) rep(NA_real_, length(cf)) else sqrt(diag(v))
    z <- cf / se
    cbind(Estimate = cf, `Std. Error` = se, `z value` = z,
          `Pr(>|z|)` = 2 * stats::pnorm(-abs(z)))
  }
  tables <- lapply(seq_along(object$tau), function(j)
    tab(object$coefficients[, j], object$vcov[[j]]))
  names(tables) <- colnames(object$coefficients)
  out <- list(call = object$call, label = qr_label(object), coefficients = tables,
              model = object$model, tau = object$tau, nobs = object$nobs,
              n_groups = object$n_groups, vcov_type = object$vcov_type, R = object$R,
              dropped = object$dropped, q = object$q)
  if (object$model == "mmqr") {
    out$location <- tab(object$location$coefficients, object$location$vcov)
    out$scale <- tab(object$scale$coefficients, object$scale$vcov)
  }
  class(out) <- "summary.panqr"
  out
}

#' @export
print.summary.panqr <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nCall:\n")
  print(x$call)
  cat("\nEstimator:", x$label, "\n")
  cat("N =", x$nobs, " groups =", x$n_groups, " vcov:", x$vcov_type,
      if (identical(x$vcov_type, "bootstrap")) paste0("(cluster, R = ", x$R, ")") else "", "\n")
  if (length(x$dropped)) cat("Dropped (time invariant):", paste(x$dropped, collapse = ", "), "\n")
  if (!is.null(x$location)) {
    cat("\nLocation (within regression, cluster-robust):\n")
    stats::printCoefmat(x$location, digits = digits, signif.stars = TRUE, has.Pvalue = TRUE)
    cat("\nScale (within regression, cluster-robust):\n")
    stats::printCoefmat(x$scale, digits = digits, signif.stars = TRUE, has.Pvalue = TRUE)
  }
  for (j in seq_along(x$tau)) {
    cat("\nQuantile tau = ", format(x$tau[j]),
        if (!is.null(x$q)) paste0(" (q(tau) = ", format(x$q[j], digits = digits), ")") else "",
        ":\n", sep = "")
    stats::printCoefmat(x$coefficients[[j]], digits = digits, signif.stars = TRUE,
                        has.Pvalue = TRUE, na.print = "")
  }
  cat("\n")
  invisible(x)
}
