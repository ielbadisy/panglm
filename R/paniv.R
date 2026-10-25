#' Instrumental-variable estimators for linear panel models
#'
#' Two-stage least squares for linear panel models with endogenous
#' regressors: fixed-effects 2SLS, random-effects generalized 2SLS
#' (Balestra and Varadharajan-Krishnakumar) or error-components 2SLS
#' (Baltagi), pooled 2SLS, and the Hausman and Taylor (1981) estimator,
#' which also identifies the effects of time-invariant regressors.
#'
#' @section Formula:
#' The formula has two parts separated by `|`. For `model = "pooling"`,
#' `"within"` and `"random"` the second part lists the complete instrument
#' set, as in [plm::plm()]: exogenous regressors are repeated on both sides
#' and excluded instruments appear only on the right, so that
#' `y ~ x1 + x2 | z1 + x2` treats `x1` as endogenous, instrumented by `z1`.
#' Regressors absent from the second part are the endogenous ones.
#'
#' For `model = "hausmantaylor"` the conventions of [plm::pht()] are used:
#' the first part lists all regressors and the second part lists those that
#' are exogenous (uncorrelated with the individual effect), e.g.
#' `y ~ x1 + x2 + z1 + z2 | x1 + z1` makes `x2` and `z2` endogenous.
#' Whether a regressor is time varying or time invariant is detected from
#' the data. Identification requires at least as many exogenous time-varying
#' regressors as endogenous time-invariant ones.
#'
#' @section Estimators:
#' * `"within"`: 2SLS on within-demeaned regressors and instruments
#'   (one-way or, with `effect = "twoways"`, two-way demeaning). Time
#'   invariant columns are dropped. The classical covariance uses the
#'   degrees of freedom `n - k - N` (`N` individual effects), as in `plm`.
#' * `"random"`: the variance components are estimated in the Swamy-Arora
#'   way from within and between 2SLS residuals (the unbalanced-panel
#'   degrees of freedom of `plm` are used, which reduce to the usual ones
#'   for balanced panels). The data are then quasi-demeaned with
#'   `theta_i = 1 - sqrt(s2_nu / (s2_nu + T_i s2_mu))`. With
#'   `inst.method = "bvk"` (the default, G2SLS) the instruments are
#'   quasi-demeaned as well; with `inst.method = "baltagi"` (EC2SLS) the
#'   instrument matrix is the within and between transforms of the
#'   instruments.
#' * `"pooling"`: ordinary 2SLS on the stacked data.
#' * `"hausmantaylor"`: the variance components follow [plm::pht()]:
#'   `s2_nu` is the within residual variance, and `s2_1` comes from a 2SLS
#'   regression of the estimated individual effects on the time-invariant
#'   regressors, instrumented by the exogenous time-invariant and
#'   time-varying regressors. The GLS-transformed regression is then
#'   estimated by 2SLS with instruments: within transforms of all
#'   time-varying regressors, individual means of the exogenous time-varying
#'   regressors, the exogenous time-invariant regressors and a constant.
#'   For unbalanced panels `s2_mu = (s2_1 - s2_nu) / T_h` with `T_h` the
#'   harmonic mean of the `T_i`, as in [plm::pht()]. A negative `s2_mu` is
#'   set to zero (with a warning), which reduces the estimator to pooled
#'   2SLS on the same instruments.
#'
#' In every case the classical covariance is `s2 inverse(Xhat'Xhat)` with `s2`
#' the mean square of the structural residuals of the transformed
#' regression. The cluster covariance is the sandwich
#' `B (sum_g u_g u_g') B` with `u_g = Xhat_g' e_g` summed within the panel
#' individual (or the supplied cluster) and the factor
#' `G / (G - 1) * (n - 1) / (n - k)`, the convention of the rest of the
#' package and of `plm::vcovHC(type = "sss")`.
#'
#' @section Diagnostics:
#' For each endogenous regressor, the first-stage regression of its
#' transformed values on the transformed instruments is compared with the
#' regression on the transformed included exogenous regressors only. The
#' result reports the partial R-squared of the excluded instruments and the
#' F statistic for their joint significance (with the absorbed fixed
#' effects counted in the residual degrees of freedom); for
#' `vcov = "cluster"` a cluster-robust Wald F is reported as well. These
#' first-stage statistics are not computed for `"hausmantaylor"`, whose
#' instruments are built from the regressors themselves. When the model is
#' overidentified, the Sargan statistic `n e'P_W e / e'e` (classical) or
#' the cluster-robust Hansen J statistic of the two-step efficient GMM
#' estimator (cluster) tests the validity of the overidentifying
#' restrictions; it is chi-squared with `L - k` degrees of freedom under the
#' null.
#'
#' @param formula a two-part formula `y ~ regressors | instruments`; see
#'   the Formula section.
#'   `lag(x, k)`, `lead(x, k)` and `diff(x)` are taken within each
#'   individual on the grid of sorted distinct time values (needs
#'   `index = c(id, time)`); the first waves of each individual get `NA`
#'   and leave the sample, and `lag(x, 0:2)` gives several lags
#' @param data a data.frame in long format
#' @param index length-1 or length-2 character vector giving the individual
#'   (and time) columns
#' @param model one of `"within"`, `"random"`, `"pooling"` or
#'   `"hausmantaylor"`
#' @param effect `"individual"`, or `"twoways"` (only for
#'   `model = "within"`)
#' @param inst.method instrument set for `model = "random"`: `"bvk"`
#'   (G2SLS, the default, as in `plm`) or `"baltagi"` (EC2SLS)
#' @param vcov `"classical"` or `"cluster"` (cluster-robust by individual);
#'   see [vcov.paniv()] to recompute after fitting
#' @param na.action missing-data action passed to [stats::model.frame()]
#' @return an object of class `"paniv"`, a list with components
#'   `coefficients`, `vcov`, `vcov_type`, `residuals` (structural residuals
#'   `y - X b` on the original scale, in the panel order), `model`,
#'   `effect`, `inst.method`, `endogenous` (names of the endogenous
#'   regressors), `instruments` (names of the excluded instruments, or for
#'   Hausman-Taylor the four variable groups in `varlist`), `first_stage`
#'   (a data frame with one row per endogenous regressor, or `NULL`),
#'   `overid` (the Sargan or Hansen test, or `NULL` when exactly
#'   identified), `ercomp` (variance components and `theta` for the random
#'   and Hausman-Taylor estimators), `nobs`, `n_groups`, `df.residual`,
#'   `call` and `formula`. Methods: `print`, `summary`, `coef`, `vcov`,
#'   `confint`, `predict`, `nobs`, `tidy` and `glance`.
#' @references
#' Baltagi, B. H. (1981). Simultaneous equations with error components.
#' *Journal of Econometrics*, 17, 189-200.
#'
#' Balestra, P. and Varadharajan-Krishnakumar, J. (1987). Full information
#' estimations of a system of simultaneous equations with error component
#' structure. *Econometric Theory*, 3, 223-246.
#'
#' Hausman, J. A. and Taylor, W. E. (1981). Panel data and unobservable
#' individual effects. *Econometrica*, 49, 1377-1398.
#' @family panel models
#' @seealso [summary.paniv()], [vcov.paniv()] and [predict.paniv()] for the
#'   fitted model, [tidy-methods] for `tidy()` and `glance()`.
#' @section Validation:
#' Checked in the package's validation suite on GitHub: the within, random and
#' pooled 2SLS estimators against `plm::plm()` (same `inst.method`),
#' `AER::ivreg()` and `fixest::feols()` on `plm::Crime`, balanced and
#' unbalanced; the Hausman-Taylor estimator against `plm::pht()` and
#' `plm::plm(random.method = "ht")` on `plm::Wages` (Baltagi 2013, Table 7.4).
#' @examples
#' # Simulated panel of 500 patients over 4 waves. Private insurance is
#' # endogenous for doctor visits (a time-varying health shock drives both)
#' # and is instrumented by an employer offer of insurance.
#' set.seed(1)
#' N <- 500; TT <- 4
#' d <- data.frame(id = rep(seq_len(N), each = TT), wave = rep(seq_len(TT), N))
#' alpha <- rep(rnorm(N), each = TT)       # unobserved frailty
#' s0 <- rep(rnorm(N), each = TT)          # family background
#' need <- rnorm(N * TT)                   # time-varying health shock
#' d$school <- round(12 + 2 * s0 + alpha)  # correlated with the frailty
#' d$female <- rep(rbinom(N, 1, 0.5), each = TT)
#' d$offer <- rbinom(N * TT, 1, plogis(s0))
#' d$income <- rnorm(N * TT) + 0.5 * alpha
#' d$insured <- as.numeric(d$offer + need + rnorm(N * TT) > 0.5)
#' d$visits <- 1 + 0.6 * d$insured + 0.2 * d$income + alpha + need + rnorm(N * TT)
#'
#' # the within estimator is biased upwards; FE-2SLS is not (true value 0.6)
#' coef(panglm(visits ~ insured + income, data = d, index = c("id", "wave"),
#'             model = "within"))
#' fe2sls <- paniv(visits ~ insured + income | offer + income, data = d,
#'                 index = c("id", "wave"), model = "within", vcov = "cluster")
#' summary(fe2sls)
#'
#' # Hausman-Taylor: the effect of schooling (time invariant, correlated
#' # with the frailty) on health costs, identified by the exogenous
#' # time-varying employer offer; income is endogenous (true values 0.4,
#' # 0.2, 0.3, 0.1)
#' d$cost <- 2 + 0.4 * d$offer + 0.2 * d$income + 0.3 * d$female +
#'   0.1 * d$school + alpha + rnorm(N * TT)
#' ht <- paniv(cost ~ offer + income + female + school | offer + female,
#'             data = d, index = c("id", "wave"), model = "hausmantaylor")
#' summary(ht)
#' @export
paniv <- function(formula, data, index,
                  model = c("within", "random", "pooling", "hausmantaylor"),
                  effect = c("individual", "twoways"),
                  inst.method = c("bvk", "baltagi"),
                  vcov = c("classical", "cluster"),
                  na.action = stats::na.omit) {
  model <- match.arg(model)
  effect <- match.arg(effect)
  inst.method <- match.arg(inst.method)
  vcov <- match.arg(vcov)
  if (effect == "twoways" && model != "within") {
    stop("effect = 'twoways' is only available for model = 'within'", call. = FALSE)
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
  if (effect == "twoways" && length(index) < 2L) {
    stop("effect = 'twoways' requires a two-element index = c(id, time)", call. = FALSE)
  }

  formula <- panel_guard_formula(formula)
  parts <- iv_split_formula(formula)
  mf <- panel_model_frame(parts$combined, data, index, na.action = na.action)
  model_rows <- match(rownames(mf), rownames(data))
  if (anyNA(model_rows)) stop("could not align model-frame rows with data", call. = FALSE)
  y <- stats::model.response(mf)
  if (!is.numeric(y) || any(!is.finite(y))) {
    stop("the response must be finite and numeric", call. = FALSE)
  }
  X <- stats::model.matrix(parts$regressors, data = mf)
  Z <- stats::model.matrix(parts$instruments, data = mf)

  panel_cols <- data[model_rows, index, drop = FALSE]
  if (anyNA(panel_cols)) stop("panel index columns must not contain missing values", call. = FALSE)
  panel <- build_panel_index(cbind(panel_cols, .paniv_row = seq_len(nrow(mf))), index)
  ord <- panel$data$.paniv_row
  y <- as.numeric(y[ord])
  X <- X[ord, , drop = FALSE]
  Z <- Z[ord, , drop = FALSE]
  group <- rep.int(seq_len(panel$n_groups), panel$group_size)
  time_id <- if (length(index) > 1L) panel$data[[index[2]]] else NULL

  fit <- if (model == "hausmantaylor") {
    paniv_hausmantaylor(y, X, Z, group, panel)
  } else {
    paniv_standard(y, X, Z, group, panel, model, effect, inst.method, time_id)
  }

  fit$call <- match.call()
  fit$formula <- formula
  fit$model <- model
  fit$effect <- effect
  fit$inst.method <- if (model == "random") inst.method else NULL
  fit$index <- index
  fit$nobs <- length(y)
  fit$n_groups <- panel$n_groups
  fit$cluster_id <- rep(panel$group_id, panel$group_size)
  fit$time_id <- time_id
  fit$xb <- y - fit$residuals
  fit$panel_order <- ord
  fit$na.action <- attr(mf, "na.action")
  fit$vcov_classical <- fit$vcov
  fit$vcov_type <- "classical"
  class(fit) <- "paniv"

  if (vcov == "cluster") {
    fit$vcov <- vcov.paniv(fit, type = "cluster")
    fit$vcov_type <- "cluster"
  }
  fit$first_stage <- if (is.null(fit$stage)) NULL else iv_first_stage(fit, vcov)
  fit$overid <- iv_overid(fit, vcov)
  fit
}

#' Split a two-part IV formula `y ~ x | z` by hand (no Formula dependency)
#'
#' @keywords internal
#' @noRd
iv_split_formula <- function(formula) {
  if (!inherits(formula, "formula") || length(formula) != 3L) {
    stop("'formula' must be a two-sided formula 'y ~ regressors | instruments'",
         call. = FALSE)
  }
  env <- environment(formula) %||% parent.frame()
  rhs <- formula[[3L]]
  if (!is.call(rhs) || !identical(rhs[[1L]], as.name("|"))) {
    stop("'formula' must have two parts, 'y ~ regressors | instruments'", call. = FALSE)
  }
  first <- rhs[[2L]]
  second <- rhs[[3L]]
  if ((is.call(first) && identical(first[[1L]], as.name("|"))) ||
      (is.call(second) && identical(second[[1L]], as.name("|")))) {
    stop("'formula' must have exactly two parts separated by a single '|'", call. = FALSE)
  }
  lhs <- formula[[2L]]
  list(
    regressors = stats::as.formula(call("~", first), env = env),
    instruments = stats::as.formula(call("~", second), env = env),
    combined = stats::as.formula(call("~", lhs, call("+", first, second)), env = env)
  )
}

#' Group means expanded to the observation level
#'
#' @keywords internal
#' @noRd
iv_between <- function(M, group, size) {
  M <- as.matrix(M)
  rowsum(M, group, reorder = FALSE)[group, , drop = FALSE] / size[group]
}

#' Keep the linearly independent columns of a matrix (in their order)
#'
#' @keywords internal
#' @noRd
iv_independent <- function(M, tol = 1e-7) {
  if (!ncol(M)) return(M)
  scale <- sqrt(colSums(M^2))
  keep <- scale > tol * max(1, sqrt(nrow(M)))
  M <- M[, keep, drop = FALSE]
  if (!ncol(M)) return(M)
  dec <- qr(M, tol = tol)
  M[, sort(dec$pivot[seq_len(dec$rank)]), drop = FALSE]
}

#' Core 2SLS on already transformed data
#'
#' @keywords internal
#' @noRd
iv_2sls <- function(y, X, W) {
  if (ncol(W) < ncol(X)) {
    stop("insufficient number of instruments: ", ncol(W), " instrument(s) for ",
         ncol(X), " regressor(s) after the panel transformation", call. = FALSE)
  }
  qW <- qr(W)
  Xhat <- qr.fitted(qW, X)
  if (!is.matrix(Xhat)) Xhat <- matrix(Xhat, ncol = ncol(X))
  colnames(Xhat) <- colnames(X)
  XtX <- crossprod(Xhat)
  if (qr(Xhat)$rank < ncol(X)) {
    stop("the model is not identified: the instruments do not span the regressors ",
         "(check that each endogenous regressor has a relevant excluded instrument)",
         call. = FALSE)
  }
  bread <- solve(XtX)
  b <- as.numeric(bread %*% crossprod(Xhat, y))
  names(b) <- colnames(X)
  e <- as.numeric(y - X %*% b)
  dimnames(bread) <- list(colnames(X), colnames(X))
  list(coefficients = b, residuals = e, Xhat = Xhat, bread = bread, W = W, X = X, y = y)
}

#' Transform a matrix for the within / two-way within estimator
#'
#' @keywords internal
#' @noRd
iv_within <- function(M, group, panel, effect, time_id) {
  if (effect == "twoways") {
    return(demean_twoway(M, rep(0, nrow(M)), panel$group_start, panel$group_size,
                         time_id)$X)
  }
  M - iv_between(M, group, panel$group_size)
}

#' FE-2SLS, pooled 2SLS and random-effects G2SLS / EC2SLS
#'
#' @keywords internal
#' @noRd
paniv_standard <- function(y, X, Z, group, panel, model, effect, inst.method, time_id) {
  Ti <- panel$group_size
  n <- length(y)
  N <- panel$n_groups
  endogenous <- setdiff(colnames(X), colnames(Z))
  excluded <- setdiff(colnames(Z), colnames(X))
  if (!length(endogenous)) {
    warning("every regressor appears among the instruments: no endogenous regressor, ",
            "the fit reduces to least squares", call. = FALSE)
  }
  ercomp <- NULL

  if (model == "pooling") {
    Xt <- iv_drop_columns(X, X)
    Wt <- iv_independent(Z)
    yt <- y
    absorbed <- 0
    transform <- function(M) M
    gls <- NULL
  } else if (model == "within") {
    Xw <- iv_within(X, group, panel, effect, time_id)
    Xt <- iv_drop_columns(Xw, Xw, drop_intercept = TRUE)
    Wt <- iv_independent(iv_within(Z, group, panel, effect, time_id))
    yt <- as.numeric(iv_within(matrix(y), group, panel, effect, time_id))
    n_time <- if (effect == "twoways") length(unique(time_id)) else 0L
    absorbed <- if (effect == "twoways") N + n_time - 1 else N
    transform <- function(M) iv_within(M, group, panel, effect, time_id)
  } else {
    # Swamy-Arora components from within and between 2SLS residuals
    Xw_all <- iv_within(X, group, panel, "individual", NULL)
    Xw <- iv_drop_columns(Xw_all, Xw_all, drop_intercept = TRUE)
    Ww <- iv_independent(iv_within(Z, group, panel, "individual", NULL))
    yw <- y - as.numeric(iv_between(matrix(y), group, Ti))
    fw <- iv_2sls(yw, Xw, Ww)
    Xb_all <- iv_between(X, group, Ti)
    Xb <- iv_independent(Xb_all)
    Wb <- iv_independent(iv_between(Z, group, Ti))
    yb <- as.numeric(iv_between(matrix(y), group, Ti))
    fb <- iv_2sls(yb, Xb, Wb)
    quad_w <- sum(fw$residuals^2)
    quad_b <- sum(fb$residuals^2)
    # degrees of freedom of plm's Swamy-Arora estimator: for balanced panels
    # N - K_b - 1 and T (N - K_b - 1), with K_b the number of identified
    # between slopes; for unbalanced panels the trace-based quadratic-form
    # expectation of plm (dfcor = 3), with K the number of pooled slopes
    m_w <- n - N - ncol(Xw)
    if (length(unique(Ti)) == 1L) {
      m_id_nu <- N - (ncol(Xb) - ("(Intercept)" %in% colnames(Xb))) - 1
      m_id_eta <- Ti[1L] * m_id_nu
    } else {
      K <- ncol(X) - ("(Intercept)" %in% colnames(X))
      m_id_nu <- N - K - 1
      cols <- colnames(Xb)
      m_id_eta <- n - sum(diag(solve(crossprod(Xb),
                                     crossprod(Xb * Ti[group], X[, cols, drop = FALSE]))))
    }
    s2_nu <- quad_w / m_w
    s2_mu <- max(0, (quad_b - m_id_nu * s2_nu) / m_id_eta)
    theta_i <- 1 - sqrt(s2_nu / (s2_nu + Ti * s2_mu))
    theta <- theta_i[group]
    ercomp <- list(sigma2 = c(idios = s2_nu, id = s2_mu), theta = theta_i)
    transform <- function(M) M - theta * iv_between(M, group, Ti)
    Xt <- iv_drop_columns(transform(X), X)
    yt <- as.numeric(transform(matrix(y)))
    if (inst.method == "bvk") {
      Wt <- iv_independent(transform(Z))
    } else {
      Wt <- iv_independent(iv_ec2sls_instruments(Z, group, panel))
    }
    absorbed <- 0
  }

  res <- iv_2sls(yt, Xt, Wt)
  k <- ncol(Xt)
  df_resid <- n - k - absorbed
  s2 <- sum(res$residuals^2) / df_resid
  V <- res$bread * s2
  coefs <- res$coefficients
  structural <- y - as.numeric(X[, names(coefs), drop = FALSE] %*% coefs)

  # first-stage ingredients: restricted instrument set = included exogenous
  included <- intersect(colnames(Xt), colnames(Z))
  stage <- NULL
  endo_kept <- intersect(endogenous, colnames(Xt))
  if (length(endo_kept)) {
    restricted <- if (model == "random" && inst.method == "baltagi") {
      iv_ec2sls_instruments(Z[, included, drop = FALSE], group, panel)
    } else {
      transform(Z[, included, drop = FALSE])
    }
    stage <- list(endogenous = endo_kept, W = Wt,
                  W_restricted = iv_independent(as.matrix(restricted)),
                  absorbed = absorbed)
  }

  list(coefficients = coefs, vcov = V, residuals = structural,
       residuals_transformed = res$residuals, Xhat = res$Xhat, X_transformed = Xt,
       y_transformed = yt, W = Wt, bread = res$bread, sigma2 = s2,
       df.residual = df_resid, absorbed = absorbed, endogenous = endogenous,
       instruments = excluded, ercomp = ercomp, stage = stage)
}

#' EC2SLS instrument matrix: within and between transforms of Z
#'
#' @keywords internal
#' @noRd
iv_ec2sls_instruments <- function(Z, group, panel) {
  Wq <- iv_within(Z, group, panel, "individual", NULL)
  Wp <- iv_between(Z, group, panel$group_size)
  colnames(Wq) <- paste0("within_", colnames(Z))
  colnames(Wp) <- paste0("mean_", colnames(Z))
  cbind(Wq, Wp)
}

#' Drop unidentified columns (zero after transformation or collinear)
#'
#' @keywords internal
#' @noRd
iv_drop_columns <- function(M, reference, drop_intercept = FALSE) {
  if (drop_intercept) M <- M[, setdiff(colnames(M), "(Intercept)"), drop = FALSE]
  kept <- iv_independent(M)
  dropped <- setdiff(colnames(M), colnames(kept))
  if (length(dropped)) {
    warning("dropping unidentified column(s): ", paste(dropped, collapse = ", "),
            call. = FALSE)
  }
  if (!ncol(kept)) stop("the model contains no estimable covariates", call. = FALSE)
  kept
}

#' Hausman-Taylor (1981) estimator, variance components as in plm::pht()
#'
#' @keywords internal
#' @noRd
paniv_hausmantaylor <- function(y, X, Z, group, panel) {
  Ti <- panel$group_size
  n <- length(y)
  N <- panel$n_groups
  if (!"(Intercept)" %in% colnames(X)) {
    stop("model = 'hausmantaylor' requires an intercept", call. = FALSE)
  }
  exogenous <- setdiff(colnames(Z), "(Intercept)")
  extra <- setdiff(exogenous, colnames(X))
  if (length(extra)) {
    stop("for model = 'hausmantaylor' the second part of the formula lists the exogenous ",
         "regressors and must be a subset of the first part; not regressors: ",
         paste(extra, collapse = ", "), call. = FALSE)
  }
  regressors <- setdiff(colnames(X), "(Intercept)")
  Xw <- X[, regressors, drop = FALSE] - iv_between(X[, regressors, drop = FALSE], group, Ti)
  scale <- pmax(1, sqrt(colSums(X[, regressors, drop = FALSE]^2)))
  invariant <- regressors[sqrt(colSums(Xw^2)) <= 1e-8 * scale]
  varying <- setdiff(regressors, invariant)
  x1 <- intersect(varying, exogenous)
  x2 <- setdiff(varying, exogenous)
  z1 <- intersect(invariant, exogenous)
  z2 <- setdiff(invariant, exogenous)
  if (!length(varying)) {
    stop("model = 'hausmantaylor' needs at least one time-varying regressor", call. = FALSE)
  }
  if (length(z2) > length(x1)) {
    stop("the Hausman-Taylor model is not identified: ", length(z2),
         " endogenous time-invariant regressor(s) but only ", length(x1),
         " exogenous time-varying regressor(s)", call. = FALSE)
  }

  # step 1: within estimator of the time-varying slopes
  Xv <- Xw[, varying, drop = FALSE]
  ybar <- as.numeric(iv_between(matrix(y), group, Ti))
  yw <- y - ybar
  within <- iv_2sls(yw, Xv, Xv)
  s2_nu <- sum(within$residuals^2) / (n - N)

  # step 2: individual effects on the time-invariant regressors by 2SLS
  alpha <- ybar - as.numeric(iv_between(X[, varying, drop = FALSE], group, Ti) %*%
                               within$coefficients)
  XC <- X[, c("(Intercept)", invariant), drop = FALSE]
  WC <- X[, c("(Intercept)", z1, x1), drop = FALSE]
  zo <- iv_2sls(alpha, XC, WC)
  s2_one <- sum(zo$residuals^2) / N
  T_h <- N / sum(1 / Ti)
  s2_mu <- (s2_one - s2_nu) / T_h
  if (s2_mu < 0) {
    warning("negative estimate of the individual-effect variance set to zero",
            call. = FALSE)
    s2_mu <- 0
  }
  theta_i <- 1 - sqrt(s2_nu / (s2_nu + Ti * s2_mu))
  theta <- theta_i[group]

  # step 3: 2SLS on the GLS-transformed data
  yt <- y - theta * ybar
  Xt <- X - theta * iv_between(X, group, Ti)
  Wv <- Xv
  colnames(Wv) <- paste0("within_", varying)
  Wm <- iv_between(X[, x1, drop = FALSE], group, Ti)
  colnames(Wm) <- paste0("mean_", x1)
  W <- cbind(`(Intercept)` = 1, Wv, X[, z1, drop = FALSE], Wm)
  W <- iv_independent(W)
  res <- iv_2sls(yt, Xt, W)
  k <- ncol(Xt)
  df_resid <- n - k
  s2 <- sum(res$residuals^2) / df_resid
  coefs <- res$coefficients

  list(coefficients = coefs, vcov = res$bread * s2,
       residuals = y - as.numeric(X %*% coefs),
       residuals_transformed = res$residuals, Xhat = res$Xhat, X_transformed = Xt,
       y_transformed = yt, W = W, bread = res$bread, sigma2 = s2,
       df.residual = df_resid, absorbed = 0, endogenous = c(x2, z2),
       instruments = character(0),
       varlist = list(xv = x1, nv = x2, xc = z1, nc = z2),
       ercomp = list(sigma2 = c(idios = s2_nu, id = s2_mu, one = s2_one), theta = theta_i),
       stage = NULL)
}

#' First-stage partial R-squared and F per endogenous regressor
#'
#' @keywords internal
#' @noRd
iv_first_stage <- function(fit, vcov_type) {
  st <- fit$stage
  W <- st$W
  Wr <- st$W_restricted
  n <- nrow(W)
  q <- ncol(W) - ncol(Wr)
  df_u <- n - ncol(W) - st$absorbed
  rows <- lapply(st$endogenous, function(v) {
    x <- fit$X_transformed[, v]
    qu <- qr(W)
    ru <- qr.resid(qu, x)
    rr <- if (ncol(Wr)) qr.resid(qr(Wr), x) else x
    rss_u <- sum(ru^2); rss_r <- sum(rr^2)
    partial <- (rss_r - rss_u) / rss_r
    f <- ((rss_r - rss_u) / q) / (rss_u / df_u)
    out <- data.frame(regressor = v, partial.R2 = partial, F = f, df1 = q, df2 = df_u,
                      p.value = stats::pf(f, q, df_u, lower.tail = FALSE),
                      stringsAsFactors = FALSE)
    if (vcov_type == "cluster" && q > 0) {
      # cluster-robust Wald test of the excluded instruments in the first stage
      Wx <- W[, setdiff(colnames(W), colnames(Wr)), drop = FALSE]
      Wx_r <- if (ncol(Wr)) qr.resid(qr(Wr), Wx) else Wx
      bread <- solve(crossprod(Wx_r))
      g <- Wx_r * ru
      score <- rowsum(g, fit$cluster_id)
      G <- nrow(score)
      V <- bread %*% crossprod(score) %*% bread * (G / (G - 1)) * ((n - 1) / (n - ncol(W)))
      pi_hat <- bread %*% crossprod(Wx_r, x)
      fr <- as.numeric(t(pi_hat) %*% solve(V, pi_hat)) / q
      out$F.robust <- fr
      out$p.robust <- stats::pf(fr, q, G - 1, lower.tail = FALSE)
    }
    out
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Sargan (classical) or Hansen J (cluster) overidentification test
#'
#' @keywords internal
#' @noRd
iv_overid <- function(fit, vcov_type) {
  W <- fit$W
  X <- fit$X_transformed
  y <- fit$y_transformed
  df <- ncol(W) - ncol(X)
  if (df < 1L) return(NULL)
  e <- fit$residuals_transformed
  if (vcov_type == "classical") {
    stat <- length(e) * sum(qr.fitted(qr(W), e)^2) / sum(e^2)
    name <- "Sargan"
  } else {
    h <- rowsum(W * e, fit$cluster_id)
    S <- crossprod(h)
    Si <- tryCatch(solve(S), error = function(err) NULL)
    if (is.null(Si)) return(NULL)
    WX <- crossprod(W, X); Wy <- crossprod(W, y)
    b2 <- solve(t(WX) %*% Si %*% WX, t(WX) %*% Si %*% Wy)
    g <- crossprod(W, y - X %*% b2)
    stat <- as.numeric(t(g) %*% Si %*% g)
    name <- "Hansen J"
  }
  list(method = name, statistic = stat, df = df,
       p.value = stats::pchisq(stat, df, lower.tail = FALSE))
}

#' Covariance of a paniv fit
#'
#' @param object a `"paniv"` fit
#' @param type `NULL` (the covariance stored on the fit), `"classical"` or
#'   `"cluster"`
#' @param cluster optional clustering vector with one value per observation
#'   used in the fit (in the row order of the data); defaults to the panel
#'   individual
#' @param ... unused
#' @return the covariance matrix of the coefficients
#' @family paniv methods
#' @examples
#' set.seed(2)
#' d <- data.frame(id = rep(1:50, each = 4), t = rep(1:4, 50))
#' d$z <- rnorm(200); d$u <- rnorm(200)
#' d$x <- d$z + d$u + rnorm(200); d$y <- 1 + d$x + d$u + rnorm(200)
#' fit <- paniv(y ~ x | z, data = d, index = c("id", "t"), model = "pooling")
#' vcov(fit, type = "cluster")
#' @export
vcov.paniv <- function(object, type = NULL, cluster = NULL, ...) {
  if (is.null(type)) return(object$vcov)
  type <- match.arg(type, c("classical", "cluster"))
  if (type == "classical") return(object$vcov_classical)
  cl <- panglm_cluster_vector(object, cluster)
  score <- object$Xhat * object$residuals_transformed
  score_sum <- rowsum(score, cl)
  G <- nrow(score_sum); n <- nrow(score); k <- ncol(score)
  if (G < 2L) stop("cluster covariance requires at least two clusters", call. = FALSE)
  v <- object$bread %*% crossprod(score_sum) %*% object$bread
  v <- v * (G / (G - 1)) * ((n - 1) / (n - k))
  dimnames(v) <- dimnames(object$vcov_classical)
  v
}

#' @export
coef.paniv <- function(object, ...) object$coefficients

#' @export
nobs.paniv <- function(object, ...) object$nobs

#' @export
confint.paniv <- function(object, parm, level = 0.95, ...) {
  confint.panglm(object, parm, level = level, ...)
}

#' @export
print.paniv <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("paniv(", iv_label(x), ")\n\n", sep = "")
  print(x$coefficients, digits = digits)
  cat("\nN =", x$nobs, " groups =", x$n_groups, "\n")
  invisible(x)
}

iv_label <- function(x) {
  switch(x$model,
    within = if (identical(x$effect, "twoways")) "two-way fixed-effects 2SLS" else "fixed-effects 2SLS",
    pooling = "pooled 2SLS",
    random = if (identical(x$inst.method, "baltagi")) "random-effects EC2SLS (Baltagi)" else
      "random-effects G2SLS (Balestra-Varadharajan-Krishnakumar)",
    hausmantaylor = "Hausman-Taylor"
  )
}

#' Summarize a paniv fit
#'
#' @param object a `"paniv"` fit
#' @param ... unused
#' @return an object of class `"summary.paniv"`: a list with the coefficient
#'   table (`coefficients`: estimate, standard error, z value and p-value),
#'   `first_stage`, `overid`, `ercomp`, `endogenous`, `instruments`,
#'   `varlist` (Hausman-Taylor) and the fit information (`call`, `model`,
#'   `nobs`, `n_groups`, `df.residual`, `vcov_type`).
#' @family paniv methods
#' @examples
#' set.seed(2)
#' d <- data.frame(id = rep(1:50, each = 4), t = rep(1:4, 50))
#' d$z <- rnorm(200); d$u <- rnorm(200)
#' d$x <- d$z + d$u + rnorm(200); d$y <- 1 + d$x + d$u + rnorm(200)
#' summary(paniv(y ~ x | z, data = d, index = c("id", "t"), model = "within"))
#' @export
summary.paniv <- function(object, ...) {
  coefs <- object$coefficients
  se <- sqrt(diag(object$vcov))
  z <- coefs / se
  table <- cbind(Estimate = coefs, `Std. Error` = se, `z value` = z,
                 `Pr(>|z|)` = 2 * stats::pnorm(-abs(z)))
  out <- list(call = object$call, label = iv_label(object), coefficients = table,
              first_stage = object$first_stage, overid = object$overid,
              ercomp = object$ercomp, endogenous = object$endogenous,
              instruments = object$instruments, varlist = object$varlist,
              model = object$model, nobs = object$nobs, n_groups = object$n_groups,
              df.residual = object$df.residual, vcov_type = object$vcov_type)
  class(out) <- "summary.paniv"
  out
}

#' @export
print.summary.paniv <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nCall:\n")
  print(x$call)
  cat("\nEstimator:", x$label, "\n")
  cat("N =", x$nobs, " groups =", x$n_groups, " vcov:", x$vcov_type, "\n")
  if (!is.null(x$varlist)) {
    show <- function(v) if (length(v)) paste(v, collapse = ", ") else "(none)"
    cat("Exogenous time-varying:  ", show(x$varlist$xv), "\n")
    cat("Endogenous time-varying: ", show(x$varlist$nv), "\n")
    cat("Exogenous time-invariant:", show(x$varlist$xc), "\n")
    cat("Endogenous time-invariant:", show(x$varlist$nc), "\n")
  } else {
    cat("Endogenous:", if (length(x$endogenous)) paste(x$endogenous, collapse = ", ") else "(none)",
        "\nExcluded instruments:",
        if (length(x$instruments)) paste(x$instruments, collapse = ", ") else "(none)", "\n")
  }
  cat("\nCoefficients:\n")
  stats::printCoefmat(x$coefficients, digits = digits, signif.stars = TRUE, has.Pvalue = TRUE)
  if (!is.null(x$ercomp)) {
    cat("\nVariance components: idiosyncratic = ",
        format(x$ercomp$sigma2[["idios"]], digits = digits), ", individual = ",
        format(x$ercomp$sigma2[["id"]], digits = digits), "\n", sep = "")
    th <- unique(round(x$ercomp$theta, 10))
    cat("theta:", if (length(th) == 1L) format(th, digits = digits) else
      paste0(format(min(th), digits = digits), " to ", format(max(th), digits = digits)), "\n")
  }
  if (!is.null(x$first_stage)) {
    cat("\nFirst stage (excluded instruments):\n")
    fs <- x$first_stage
    rownames(fs) <- fs$regressor
    fs$regressor <- NULL
    print(fs, digits = digits)
  }
  if (!is.null(x$overid)) {
    cat("\n", x$overid$method, " overidentification test: chi2(", x$overid$df, ") = ",
        format(x$overid$statistic, digits = digits), ", p = ",
        format.pval(x$overid$p.value, digits = digits), "\n", sep = "")
  }
  cat("Residual df:", x$df.residual, "\n\n")
  invisible(x)
}
