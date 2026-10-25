#' Sample-selection correction for linear panel models
#'
#' Two-step selection correction of Wooldridge (1995) for linear panel
#' models in which the outcome is observed only for a selected subsample
#' (attrition in a cohort, survival, response to a wave), when selection
#' may depend on the idiosyncratic error of the outcome equation.
#'
#' @section Estimator:
#' Step 1. For each period `t` a probit of the selection indicator `s_it`
#' on the selection regressors `z_it` and, with `mundlak = TRUE`, the
#' individual means `zbar_i` of the time-varying selection regressors (the
#' Mundlak version of Chamberlain's device) is fitted on all observations
#' of that period (by iteratively reweighted least squares, as
#' [stats::glm()], with the convergence tolerance `1e-12`). The inverse Mills ratio is
#' `lambda_it = dnorm(xb_it) / pnorm(xb_it)` with `xb_it` the probit
#' index. The means are taken over all observations of the individual in
#' the selection sample, selected or not. A period in which every
#' observation (or none) is selected has no probit; its inverse Mills
#' ratio is set to zero and it gets no correction term (for instance the
#' first wave of a cohort, before any attrition).
#'
#' Step 2. The outcome equation is fitted on the selected observations,
#' adding the inverse Mills ratio interacted with period dummies,
#' `lambda_it * d_t` (`lambda = "period"`, Wooldridge's specification, one
#' coefficient per period with a probit) or a single term `lambda_it`
#' (`lambda = "common"`). With `model = "within"` the individual effects
#' are removed by the within transformation on the selected sample; with
#' `model = "pooling"` the outcome equation is estimated by pooled OLS and,
#' when `mundlak = TRUE`, includes the same Mundlak means `zbar_i` to
#' account for the correlation between the individual effect and the
#' regressors. Period effects are not added automatically; include them in
#' `outcome` (e.g. `+ factor(wave)`) when needed.
#'
#' The selection regressors should contain the outcome regressors plus at
#' least one excluded variable (an instrument that affects selection but
#' not the outcome); without one the correction is identified only by the
#' nonlinearity of the inverse Mills ratio, and a warning is given.
#'
#' @section Inference:
#' The coefficient standard errors default to a panel (cluster) bootstrap:
#' whole individuals are resampled with replacement and both steps,
#' including the Mundlak means and the per-period probits, are re-estimated
#' in each replicate. The covariance is the sample covariance of the `boot`
#' replicate estimates; replicates that fail are dropped and counted.
#'
#' The cluster-robust covariance `vcov(fit, type = "cluster")` treats the
#' inverse Mills ratios as known regressors. It ignores the estimation of
#' the first step and is therefore not valid for the coefficients when
#' there is selection bias. It uses the factor
#' `G / (G - 1) * (n - 1) / (n - k)`, as in `plm::vcovHC(type = "sss")`.
#'
#' The test of selection bias is the joint Wald test that all inverse
#' Mills ratio coefficients are zero. Under this null hypothesis the first
#' step does not affect the asymptotic distribution of the second-step
#' estimator, so the cluster-robust Wald statistic is valid (Wooldridge
#' 2010, section 19.9); it is chi-squared with as many degrees of freedom
#' as correction terms. The same statistic computed with the bootstrap
#' covariance is reported as well.
#'
#' @param selection a two-sided formula `s ~ z1 + x1 + ...` for the
#'   selection equation; `s` is a 0/1 (or logical) indicator equal to one
#'   when the outcome is observed
#' @param outcome a two-sided formula `y ~ x1 + ...` for the outcome
#'   equation, fitted on the observations with `s == 1`
#' @param data a data.frame in long format, with one row per individual and
#'   period whether or not the outcome is observed
#' @param index length-2 character vector giving the individual and period
#'   columns
#' @param model `"within"` (fixed effects in the outcome equation) or
#'   `"pooling"` (pooled OLS, with Mundlak means when `mundlak = TRUE`)
#' @param mundlak logical; add the individual means of the time-varying
#'   selection regressors to the probits (and to the pooled outcome
#'   equation)
#' @param lambda `"period"` (inverse Mills ratio interacted with period
#'   dummies) or `"common"` (a single inverse Mills ratio term)
#' @param boot number of panel bootstrap replicates; `0` skips the
#'   bootstrap and makes the cluster-robust covariance the default
#' @param seed optional random seed for the bootstrap; the global random
#'   number state is restored afterwards
#' @param na.action missing-data action passed to [stats::model.frame()]
#'   for each formula
#' @return an object of class `"panselect"`, a list with components
#'   `coefficients`, `vcov` (the default covariance), `vcov_type`
#'   (`"boot"` or `"cluster"`), `vcov_boot`, `vcov_cluster`,
#'   `vcov_classical`, `boot` (the replicate estimates, the number of
#'   requested and failed replicates), `probit` (a list with the probit
#'   coefficients of each period, `NULL` for a period without probit),
#'   `lambda_terms` (names of the correction terms), `selection_test` (the
#'   cluster-robust and bootstrap Wald tests), `imr` (the inverse Mills
#'   ratio of every observation of the selection sample, in the row order
#'   of the data), `residuals`, `nobs` (selected observations used),
#'   `nobs_selection`, `n_groups`, `n_groups_selection`, `df.residual`,
#'   `periods`, `model`, `mundlak`, `lambda`, `call`. It has `print`,
#'   `summary`, `coef`, `vcov`, `nobs` and `confint` methods.
#' @references
#' Wooldridge, J. M. (1995). Selection corrections for panel data models
#' under conditional mean independence assumptions. *Journal of
#' Econometrics*, 68, 115-132.
#'
#' Wooldridge, J. M. (2010). *Econometric Analysis of Cross Section and
#' Panel Data*, 2nd edition, section 19.9. MIT Press.
#' @family panel models
#' @seealso [summary.panselect()], [vcov.panselect()] and [predict.panselect()]
#'   for the fitted model, [tidy-methods] for `tidy()` and `glance()`.
#' @examples
#' # Simulated cohort of 300 patients followed over 4 waves. A quality of
#' # life score is observed only for patients who answer the questionnaire;
#' # non-response depends on the same health shock that lowers the score,
#' # and on a reminder letter that does not affect the score itself.
#' set.seed(1)
#' N <- 300; TT <- 4
#' d <- data.frame(id = rep(seq_len(N), each = TT), wave = rep(seq_len(TT), N))
#' alpha <- rep(rnorm(N), each = TT)
#' d$treat <- rnorm(N * TT) + 0.5 * alpha
#' d$reminder <- rnorm(N * TT)
#' v <- rnorm(N * TT)
#' u <- 0.7 * v + rnorm(N * TT, sd = 0.7)
#' d$respond <- as.numeric(0.3 + 0.5 * d$treat + d$reminder + 0.3 * alpha + v > 0)
#' d$qol <- ifelse(d$respond == 1, 1 + 0.5 * d$treat + alpha + u, NA)
#'
#' fit <- panselect(respond ~ treat + reminder, qol ~ treat, data = d,
#'                  index = c("id", "wave"), model = "within", boot = 49, seed = 2)
#' summary(fit)
#' @export
panselect <- function(selection, outcome, data, index,
                      model = c("within", "pooling"),
                      mundlak = TRUE,
                      lambda = c("period", "common"),
                      boot = 199L, seed = NULL,
                      na.action = stats::na.omit) {
  model <- match.arg(model)
  lambda <- match.arg(lambda)
  if (!is.logical(mundlak) || length(mundlak) != 1L || is.na(mundlak)) {
    stop("'mundlak' must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.numeric(boot) || length(boot) != 1L || is.na(boot) || boot < 0 ||
      boot != round(boot) || boot == 1) {
    stop("'boot' must be 0 or an integer of at least 2", call. = FALSE)
  }
  for (f in list(selection, outcome)) {
    if (!inherits(f, "formula") || length(f) != 3L) {
      stop("'selection' and 'outcome' must be two-sided formulas", call. = FALSE)
    }
  }
  if (missing(index)) stop("'index' is required, e.g. index = c(\"id\", \"time\")", call. = FALSE)
  if (!is.character(index) || length(index) != 2L || anyDuplicated(index)) {
    stop("'index' must contain two distinct column names, c(id, period)", call. = FALSE)
  }
  missing_index <- setdiff(index, names(data))
  if (length(missing_index)) {
    stop("index column(s) not found in data: ", paste(missing_index, collapse = ", "),
         call. = FALSE)
  }
  data <- as.data.frame(data)
  rownames(data) <- NULL

  # selection sample: every row with complete selection variables
  mf_s <- stats::model.frame(selection, data = data, na.action = na.action)
  rows_s <- match(rownames(mf_s), rownames(data))
  if (anyNA(data[rows_s, index])) {
    stop("panel index columns must not contain missing values", call. = FALSE)
  }
  s <- stats::model.response(mf_s)
  if (is.logical(s)) s <- as.numeric(s)
  if (!is.numeric(s) || !all(s %in% c(0, 1))) {
    stop("the selection response must be a 0/1 or logical indicator", call. = FALSE)
  }
  Z <- stats::model.matrix(selection, data = mf_s)
  if (!"(Intercept)" %in% colnames(Z)) {
    stop("the selection equation must include an intercept", call. = FALSE)
  }

  # outcome sample: selected rows with complete outcome variables
  sel_rows <- rows_s[s == 1]
  if (!length(sel_rows)) stop("no observation is selected (s == 1)", call. = FALSE)
  mf_o <- stats::model.frame(outcome, data = data[sel_rows, , drop = FALSE],
                             na.action = na.action)
  rows_o <- match(rownames(mf_o), rownames(data))
  y_o <- stats::model.response(mf_o)
  if (!is.numeric(y_o) || any(!is.finite(y_o))) {
    stop("the outcome must be finite and numeric", call. = FALSE)
  }
  X_o <- stats::model.matrix(outcome, data = mf_o)

  # panel order of the selection sample
  panel <- build_panel_index(cbind(data[rows_s, index, drop = FALSE],
                                   .ps_row = seq_along(rows_s)), index)
  ord <- panel$data$.ps_row
  Z <- Z[ord, , drop = FALSE]
  s <- as.numeric(s[ord])
  data_row <- rows_s[ord]
  gid <- rep.int(seq_len(panel$n_groups), panel$group_size)
  tt <- panel$data[[index[2]]]
  periods <- sort(unique(tt))

  # outcome values aligned with the selection sample (NA when not used)
  pos <- match(data_row, rows_o)
  use <- !is.na(pos)
  if (!any(use)) stop("no selected observation has a complete outcome", call. = FALSE)
  X <- matrix(NA_real_, length(s), ncol(X_o), dimnames = list(NULL, colnames(X_o)))
  X[use, ] <- X_o[pos[use], , drop = FALSE]
  y <- rep(NA_real_, length(s))
  y[use] <- as.numeric(y_o[pos[use]])

  # time-varying selection regressors (Mundlak means) on the full sample
  zv <- setdiff(colnames(Z), "(Intercept)")
  varying <- zv[vapply(zv, function(v) {
    m <- iv_between(Z[, v], gid, panel$group_size)
    any(abs(Z[, v] - m) > 1e-10 * max(1, abs(Z[, v])))
  }, logical(1))]
  mean_cols <- if (mundlak) varying else character(0)

  xo <- setdiff(colnames(X_o), "(Intercept)")
  if (!length(setdiff(zv, xo))) {
    warning("the selection equation has no regressor excluded from the outcome ",
            "equation; the correction is identified only by the nonlinearity of ",
            "the inverse Mills ratio", call. = FALSE)
  }
  if (model == "within" && length(periods) < 2L) {
    stop("model = 'within' requires at least two periods", call. = FALSE)
  }

  spec <- list(model = model, lambda = lambda, mean_cols = mean_cols, periods = periods)
  fit <- ps_core(Z, s, gid, tt, y, X, use, spec)

  # cluster-robust covariance (first step treated as known)
  fit$vcov_cluster <- ps_cluster_vcov(fit$X_transformed, fit$residuals_transformed,
                                      fit$bread, fit$cluster_id)
  dimnames(fit$vcov_cluster) <- dimnames(fit$vcov_classical)

  # panel bootstrap over individuals
  boot_res <- NULL
  if (boot > 0) {
    boot_res <- ps_bootstrap(Z, s, gid, tt, y, X, use, spec, panel$n_groups,
                             as.integer(boot), seed, names(fit$coefficients))
    keep <- stats::complete.cases(boot_res$draws)
    if (sum(keep) < 2L) stop("fewer than two bootstrap replicates succeeded", call. = FALSE)
    fit$vcov_boot <- stats::cov(boot_res$draws[keep, , drop = FALSE])
  } else {
    fit$vcov_boot <- NULL
  }
  fit$boot <- boot_res
  fit$vcov_type <- if (boot > 0) "boot" else "cluster"
  fit$vcov <- if (boot > 0) fit$vcov_boot else fit$vcov_cluster

  fit$selection_test <- list(
    cluster = ps_wald(fit$coefficients, fit$vcov_cluster, fit$lambda_terms),
    boot = if (boot > 0) ps_wald(fit$coefficients, fit$vcov_boot, fit$lambda_terms) else NULL
  )
  imr <- rep(NA_real_, nrow(data))
  imr[data_row] <- fit$imr
  fit$imr <- imr
  resid <- rep(NA_real_, nrow(data))
  resid[data_row[use]] <- fit$residuals
  fit$residuals <- resid[!is.na(resid)]
  fit$call <- match.call()
  fit$selection <- selection
  fit$outcome <- outcome
  fit$model <- model
  fit$mundlak <- mundlak
  fit$lambda <- lambda
  fit$index <- index
  fit$periods <- periods
  fit$mean_cols <- mean_cols
  fit$nobs_selection <- length(s)
  fit$n_groups_selection <- panel$n_groups
  fit$na.action <- list(selection = attr(mf_s, "na.action"), outcome = attr(mf_o, "na.action"))
  class(fit) <- "panselect"
  fit
}

#' Both estimation steps on panel-ordered arrays
#'
#' @keywords internal
#' @noRd
ps_core <- function(Z, s, gid, tt, y, X, use, spec, boot = FALSE) {
  G <- max(gid)
  size <- tabulate(gid, G)
  M <- NULL
  if (length(spec$mean_cols)) {
    M <- iv_between(Z[, spec$mean_cols, drop = FALSE], gid, size)
    colnames(M) <- paste0("mean_", spec$mean_cols)
  }
  D <- cbind(Z, M)

  # step 1: one probit per period
  imr <- rep(0, length(s))
  probit <- stats::setNames(vector("list", length(spec$periods)), as.character(spec$periods))
  has_probit <- logical(length(spec$periods))
  for (j in seq_along(spec$periods)) {
    r <- which(tt == spec$periods[j])
    if (!length(r) || all(s[r] == 1) || all(s[r] == 0)) next
    Dt <- iv_independent(D[r, , drop = FALSE])
    pf <- if (boot) {
      suppressWarnings(ps_probit(Dt, s[r]))
    } else {
      ps_probit(Dt, s[r])
    }
    if (!pf$converged) {
      if (boot) stop("probit did not converge")
      warning("the probit for period ", spec$periods[j], " did not converge", call. = FALSE)
    }
    b <- stats::setNames(rep(NA_real_, ncol(D)), colnames(D))
    b[colnames(Dt)] <- pf$coefficients
    probit[[j]] <- b
    xb <- as.numeric(Dt %*% pf$coefficients)
    imr[r] <- exp(stats::dnorm(xb, log = TRUE) - stats::pnorm(xb, log.p = TRUE))
    has_probit[j] <- TRUE
  }
  if (!any(has_probit)) {
    stop("every period has all observations selected (or none): no selection to correct",
         call. = FALSE)
  }

  # step 2: outcome equation on the selected sample
  L <- if (spec$lambda == "common") {
    matrix(imr[use], ncol = 1L, dimnames = list(NULL, "lambda"))
  } else {
    pp <- spec$periods[has_probit]
    out <- vapply(pp, function(p) imr[use] * (tt[use] == p), numeric(sum(use)))
    out <- matrix(out, nrow = sum(use))
    colnames(out) <- paste0("lambda:", pp)
    out
  }
  W <- X[use, , drop = FALSE]
  if (spec$model == "pooling" && !is.null(M)) W <- cbind(W, M[use, , drop = FALSE])
  W <- cbind(W, L)
  yy <- y[use]
  g_o <- gid[use]
  g_o <- match(g_o, unique(g_o))
  N_o <- max(g_o)
  if (spec$model == "within") {
    W <- W[, colnames(W) != "(Intercept)", drop = FALSE]
    size_o <- tabulate(g_o, N_o)
    W <- W - iv_between(W, g_o, size_o)
    yy <- yy - as.numeric(iv_between(yy, g_o, size_o))
  }
  W <- iv_independent(W)
  lambda_terms <- intersect(colnames(L), colnames(W))
  if (!length(lambda_terms)) {
    stop("the inverse Mills ratio terms are collinear with the outcome regressors",
         call. = FALSE)
  }
  qw <- qr(W)
  b <- stats::setNames(as.numeric(qr.coef(qw, yy)), colnames(W))
  e <- as.numeric(yy - W %*% b)
  if (boot) return(b)
  n <- length(yy); k <- ncol(W)
  df <- n - k - if (spec$model == "within") N_o else 0L
  bread <- chol2inv(qr.R(qw))
  bread <- bread[order(qw$pivot), order(qw$pivot), drop = FALSE]
  dimnames(bread) <- list(colnames(W), colnames(W))
  list(coefficients = b, residuals = e, residuals_transformed = e,
       X_transformed = W, bread = bread,
       vcov_classical = sum(e^2) / df * bread,
       df.residual = df, nobs = n, n_groups = N_o, cluster_id = g_o,
       lambda_terms = lambda_terms, probit = probit,
       probit_periods = spec$periods[has_probit], imr = imr)
}

#' Probit by IRLS with a tight convergence criterion
#'
#' @keywords internal
#' @noRd
ps_probit <- function(X, s) {
  stats::glm.fit(X, s, family = stats::binomial(link = "probit"),
                 control = stats::glm.control(epsilon = 1e-12, maxit = 100))
}

#' Cluster-robust sandwich with the sss small-sample factor
#'
#' @keywords internal
#' @noRd
ps_cluster_vcov <- function(W, e, bread, cluster) {
  score <- rowsum(W * e, cluster)
  G <- nrow(score); n <- nrow(W); k <- ncol(W)
  if (G < 2L) stop("cluster covariance requires at least two clusters", call. = FALSE)
  bread %*% crossprod(score) %*% bread * (G / (G - 1)) * ((n - 1) / (n - k))
}

#' Panel bootstrap: resample whole individuals and refit both steps
#'
#' @keywords internal
#' @noRd
ps_bootstrap <- function(Z, s, gid, tt, y, X, use, spec, N, B, seed, terms) {
  if (!is.null(seed)) {
    old <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      get(".Random.seed", envir = globalenv(), inherits = FALSE)
    } else NULL
    on.exit({
      if (is.null(old)) {
        rm(".Random.seed", envir = globalenv())
      } else {
        assign(".Random.seed", old, envir = globalenv())
      }
    }, add = TRUE)
    set.seed(seed)
  }
  rows_by_group <- split(seq_along(gid), gid)
  draws <- matrix(NA_real_, B, length(terms), dimnames = list(NULL, terms))
  failed <- 0L
  for (b in seq_len(B)) {
    pick <- sample.int(N, N, replace = TRUE)
    r <- unlist(rows_by_group[pick], use.names = FALSE)
    g <- rep.int(seq_len(N), lengths(rows_by_group)[pick])
    est <- tryCatch(
      ps_core(Z[r, , drop = FALSE], s[r], g, tt[r], y[r], X[r, , drop = FALSE],
              use[r], spec, boot = TRUE),
      error = function(err) NULL)
    if (is.null(est) || !all(terms %in% names(est))) {
      failed <- failed + 1L
      next
    }
    draws[b, ] <- est[terms]
  }
  list(draws = draws, R = B, failed = failed)
}

#' Joint Wald test that the inverse Mills ratio coefficients are zero
#'
#' @keywords internal
#' @noRd
ps_wald <- function(coefs, V, terms) {
  b <- coefs[terms]
  stat <- as.numeric(t(b) %*% solve(V[terms, terms, drop = FALSE], b))
  df <- length(terms)
  list(statistic = stat, df = df, p.value = stats::pchisq(stat, df, lower.tail = FALSE))
}

#' Covariance of a panselect fit
#'
#' @param object a `"panselect"` fit
#' @param type `NULL` (the default covariance of the fit), `"boot"` (panel
#'   bootstrap), `"cluster"` (cluster-robust by individual, ignoring the
#'   estimation of the first step) or `"classical"` (homoskedastic, also
#'   ignoring the first step)
#' @param ... unused
#' @return the covariance matrix of the outcome-equation coefficients
#' @family panselect methods
#' @examples
#' set.seed(3)
#' d <- data.frame(id = rep(1:100, each = 3), t = rep(1:3, 100))
#' d$x <- rnorm(300); d$z <- rnorm(300); v <- rnorm(300)
#' d$s <- as.numeric(0.5 + d$x + d$z + v > 0)
#' d$y <- ifelse(d$s == 1, 1 + d$x + 0.5 * v + rnorm(300), NA)
#' fit <- panselect(s ~ x + z, y ~ x, data = d, index = c("id", "t"), boot = 0)
#' vcov(fit, type = "cluster")
#' @export
vcov.panselect <- function(object, type = NULL, ...) {
  if (is.null(type)) return(object$vcov)
  type <- match.arg(type, c("boot", "cluster", "classical"))
  switch(type,
    boot = {
      if (is.null(object$vcov_boot)) {
        stop("no bootstrap covariance: the model was fitted with boot = 0", call. = FALSE)
      }
      object$vcov_boot
    },
    cluster = object$vcov_cluster,
    classical = object$vcov_classical)
}

#' @export
coef.panselect <- function(object, ...) object$coefficients

#' @export
nobs.panselect <- function(object, ...) object$nobs

#' @export
confint.panselect <- function(object, parm, level = 0.95, ...) {
  confint.panglm(object, parm, level = level, ...)
}

ps_label <- function(x) {
  paste0(if (x$model == "within") "fixed effects" else "pooled OLS",
         if (x$mundlak) " with Mundlak means" else "",
         ", inverse Mills ratio ",
         if (x$lambda == "period") "by period" else "common to all periods")
}

#' @export
print.panselect <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("panselect(Wooldridge selection correction, ", ps_label(x), ")\n\n", sep = "")
  print(x$coefficients, digits = digits)
  cat("\nSelected N =", x$nobs, " groups =", x$n_groups,
      " (selection sample N =", x$nobs_selection, ")\n")
  invisible(x)
}

#' Summarize a panselect fit
#'
#' @param object a `"panselect"` fit
#' @param ... unused
#' @return an object of class `"summary.panselect"`: a list with the
#'   coefficient table (`coefficients`: estimate, standard error, z value
#'   and p-value, using the default covariance of the fit), the selection
#'   test (`selection_test`), the probit coefficients (`probit`) and the
#'   fit information.
#' @family panselect methods
#' @examples
#' set.seed(3)
#' d <- data.frame(id = rep(1:100, each = 3), t = rep(1:3, 100))
#' d$x <- rnorm(300); d$z <- rnorm(300); v <- rnorm(300)
#' d$s <- as.numeric(0.5 + d$x + d$z + v > 0)
#' d$y <- ifelse(d$s == 1, 1 + d$x + 0.5 * v + rnorm(300), NA)
#' summary(panselect(s ~ x + z, y ~ x, data = d, index = c("id", "t"),
#'                   boot = 29, seed = 1))
#' @export
summary.panselect <- function(object, ...) {
  coefs <- object$coefficients
  se <- sqrt(diag(object$vcov))
  z <- coefs / se
  table <- cbind(Estimate = coefs, `Std. Error` = se, `z value` = z,
                 `Pr(>|z|)` = 2 * stats::pnorm(-abs(z)))
  out <- list(call = object$call, label = ps_label(object), coefficients = table,
              selection_test = object$selection_test, probit = object$probit,
              probit_periods = object$probit_periods, lambda_terms = object$lambda_terms,
              nobs = object$nobs, n_groups = object$n_groups,
              nobs_selection = object$nobs_selection,
              n_groups_selection = object$n_groups_selection,
              df.residual = object$df.residual, vcov_type = object$vcov_type,
              boot = if (is.null(object$boot)) NULL else
                c(R = object$boot$R, failed = object$boot$failed))
  class(out) <- "summary.panselect"
  out
}

#' @export
print.summary.panselect <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nCall:\n")
  print(x$call)
  cat("\nWooldridge (1995) panel selection correction:", x$label, "\n")
  cat("Selection sample: N =", x$nobs_selection, " groups =", x$n_groups_selection, "\n")
  cat("Outcome sample:   N =", x$nobs, " groups =", x$n_groups, "\n")
  cat("Probits fitted for periods:", paste(x$probit_periods, collapse = ", "), "\n")
  cat("\nOutcome equation (vcov: ",
      if (x$vcov_type == "boot") {
        paste0("panel bootstrap, ", x$boot[["R"]] - x$boot[["failed"]], " replicates")
      } else {
        "cluster-robust, ignores the first-step estimation"
      }, "):\n", sep = "")
  stats::printCoefmat(x$coefficients, digits = digits, signif.stars = TRUE, has.Pvalue = TRUE)
  if (!is.null(x$boot) && x$boot[["failed"]] > 0) {
    cat("Note:", x$boot[["failed"]], "bootstrap replicate(s) failed and were dropped\n")
  }
  st <- x$selection_test$cluster
  cat("\nTest of selection bias (H0: all inverse Mills ratio coefficients are zero):\n")
  cat("  cluster-robust Wald chi2(", st$df, ") = ", format(st$statistic, digits = digits),
      ", p = ", format.pval(st$p.value, digits = digits), "\n", sep = "")
  sb <- x$selection_test$boot
  if (!is.null(sb)) {
    cat("  bootstrap Wald chi2(", sb$df, ") = ", format(sb$statistic, digits = digits),
        ", p = ", format.pval(sb$p.value, digits = digits), "\n", sep = "")
  }
  cat("Residual df:", x$df.residual, "\n\n")
  invisible(x)
}
