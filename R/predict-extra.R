# predict() methods for pangmm(), paniv() and panselect() fits.

# ---- paniv ------------------------------------------------------------------

#' Predictions from a panel instrumental-variables fit
#'
#' Linear predictions from a [paniv()] fit. For `model = "pooling"`,
#' `"random"` and `"hausmantaylor"` the prediction is the linear index
#' `x'b` (with the intercept), the conditional mean given the regressors
#' that does not involve the individual effect; this is what `predict()`
#' returns for `AER::ivreg()` and `plm::plm()` fits of these models. For
#' `model = "within"` the estimated individual effects
#' `a_i = mean_t(y_it - x_it'b)` are added (and, for `effect = "twoways"`,
#' the time effects, from the two-way decomposition of `y - x'b`), as in
#' [predict.panglm()].
#'
#' Without `newdata` the in-sample fitted values are returned in the row
#' order of the fitting data (with `NA` for the rows dropped by
#' `na.action`). With `newdata`, the regressors of the first part of the
#' formula are evaluated on `newdata`; factor regressors must carry the
#' levels used in the fit. For within fits `newdata` must contain the
#' panel index columns.
#'
#' @param object a `"paniv"` fit
#' @param newdata optional data frame
#' @param allow.new.levels for within fits, if `TRUE` an individual (or
#'   time period) of `newdata` that is not in the fit gets an effect of zero,
#'   so that its prediction is the linear index `x'b`; if `FALSE` (default)
#'   it is an error
#' @param ... unused
#' @return a numeric vector of predictions, one per row of `newdata`
#'   (`NA` where a regressor is missing)
#' @family paniv methods
#' @section Validation:
#' Checked against `predict()` of `plm::plm()` instrumental-variable fits and
#' `AER::ivreg()` in the package's validation suite on GitHub.
#' @examples
#' set.seed(1)
#' d <- data.frame(id = rep(1:50, each = 4), t = rep(1:4, 50))
#' d$z <- rnorm(200); d$u <- rnorm(200)
#' d$x <- d$z + d$u + rnorm(200)
#' d$y <- 1 + d$x + d$u + rep(rnorm(50), each = 4) + rnorm(200)
#' fit <- paniv(y ~ x | z, data = d, index = c("id", "t"), model = "within")
#' head(predict(fit))
#' predict(fit, newdata = data.frame(id = 1:2, t = 5, x = 0))
#' @export
predict.paniv <- function(object, newdata = NULL, allow.new.levels = FALSE, ...) {
  b <- object$coefficients
  within <- object$model == "within"
  if (within) {
    fe <- object$residuals - object$residuals_transformed
    ids <- object$cluster_id
  }
  if (is.null(newdata)) {
    if (is.null(object$xb)) stop("the fit does not store its linear index; refit it", call. = FALSE)
    values <- if (within) object$xb + fe else object$xb
    values_mf <- numeric(length(values))
    values_mf[object$panel_order] <- values
    return(stats::napredict(object$na.action, values_mf))
  }
  newdata <- as.data.frame(newdata)
  X <- iv_newdata_design(object, newdata)
  eta <- as.numeric(X %*% b)
  if (!within) return(eta)

  idn <- object$index[1]
  if (!idn %in% names(newdata)) {
    stop("newdata must contain the panel index column '", idn, "'", call. = FALSE)
  }
  if (object$effect == "twoways") {
    if (length(object$index) < 2L || !object$index[2] %in% names(newdata)) {
      stop("newdata must contain the time index column", call. = FALSE)
    }
    eff <- iv_twoway_effects(fe, ids, object$time_id)
    a <- unname(eff$individual[match(as.character(newdata[[idn]]), names(eff$individual))])
    g <- unname(eff$time[match(as.character(newdata[[object$index[2]]]), names(eff$time))])
    if ((anyNA(a) || anyNA(g)) && !allow.new.levels) {
      stop("newdata contains individuals or periods not present in the fitted model",
           call. = FALSE)
    }
    a[is.na(a)] <- 0
    g[is.na(g)] <- 0
    return(eta + as.numeric(a) + as.numeric(g))
  }
  alpha <- tapply(fe, as.character(ids), mean)
  a <- unname(alpha[match(as.character(newdata[[idn]]), names(alpha))])
  if (anyNA(a) && !allow.new.levels) {
    stop("newdata contains individual levels not present in the fitted model", call. = FALSE)
  }
  a[is.na(a)] <- 0
  eta + as.numeric(a)
}

# Design matrix of the regressors on newdata, columns matching coef()
iv_newdata_design <- function(object, newdata) {
  parts <- iv_split_formula(object$formula)
  mf <- stats::model.frame(parts$regressors, data = newdata, na.action = stats::na.pass)
  X <- stats::model.matrix(parts$regressors, data = mf)
  miss <- setdiff(names(object$coefficients), colnames(X))
  if (length(miss)) {
    stop("newdata does not produce the regressor column(s) ", paste(miss, collapse = ", "),
         " (check the factor levels)", call. = FALSE)
  }
  X[, names(object$coefficients), drop = FALSE]
}

# Individual and time effects from the in-sample two-way effect a_i + g_t,
# by alternating projections (exact, as fe lies in the span of the dummies)
iv_twoway_effects <- function(fe, id, time) {
  if (is.null(time)) stop("the fit does not store its time index", call. = FALSE)
  id <- as.character(id)
  time <- as.character(time)
  a <- stats::setNames(numeric(length(unique(id))), sort(unique(id)))
  g <- stats::setNames(numeric(length(unique(time))), sort(unique(time)))
  for (it in seq_len(10000L)) {
    g_new <- tapply(fe - a[id], time, mean)[names(g)]
    a_new <- tapply(fe - g_new[time], id, mean)[names(a)]
    delta <- max(abs(a_new - a), abs(g_new - g))
    a <- a_new
    g <- g_new
    if (delta < 1e-13 * max(1, max(abs(fe)))) break
  }
  # normalise the time effects to mean zero (predictions are unaffected)
  shift <- mean(g)
  list(individual = a + shift, time = g - shift)
}

# ---- panselect --------------------------------------------------------------

#' Predictions from a panel sample-selection fit
#'
#' Predictions from a [panselect()] fit on `newdata`, a long data frame with
#' the panel index columns and the variables of both equations (the
#' outcome and the selection indicator themselves are not needed).
#'
#' * `type = "unconditional"`: the linear index of the outcome equation,
#'   `x'b`, the mean of the outcome in the whole population (selected or
#'   not) given the regressors. With `model = "pooling"` it includes the
#'   intercept and, with `mundlak = TRUE`, the Mundlak means; with
#'   `model = "within"` it excludes the individual effect (and the
#'   intercept, which the within transformation removes).
#' * `type = "conditional"`: the mean in the selected population,
#'   `x'b + sum_t g_t lambda_it`, adding the inverse Mills ratio terms with
#'   their estimated coefficients (one per period with
#'   `lambda = "period"`, a single one with `lambda = "common"`).
#' * `type = "selection"`: the selection probability `pnorm(z'd_t)` of the
#'   first-step probit of each row's period.
#'
#' The inverse Mills ratios and probabilities use the first-step probit of
#' the row's period; the Mundlak means are the means of the time-varying
#' selection regressors over the rows of each individual in `newdata`, so
#' `newdata` should hold every period of an individual (as in the fit).
#' A period without a probit in the fit (every observation selected) has an
#' inverse Mills ratio of zero and a missing selection probability; a
#' period that is not in the fit gives a missing prediction for the
#' conditional and selection types.
#'
#' @param object a `"panselect"` fit
#' @param newdata a data frame in long format (see Details); the fit does
#'   not store its data, so `newdata` is required (pass the fitting data
#'   for in-sample predictions)
#' @param type `"unconditional"`, `"conditional"` or `"selection"`
#' @param ... unused
#' @return a numeric vector with one prediction per row of `newdata`
#' @family panselect methods
#' @section Validation:
#' On one period, checked against `predict()` of
#' `sampleSelection::selection(method = "2step")` in the package's validation
#' suite on GitHub.
#' @examples
#' set.seed(1)
#' d <- data.frame(id = rep(1:150, each = 3), t = rep(1:3, 150))
#' d$x <- rnorm(450); d$z <- rnorm(450); v <- rnorm(450)
#' d$s <- as.numeric(0.5 + d$x + d$z + v > 0)
#' d$y <- ifelse(d$s == 1, 1 + d$x + 0.5 * v + rnorm(450), NA)
#' fit <- panselect(s ~ x + z, y ~ x, data = d, index = c("id", "t"),
#'                  model = "pooling", boot = 0)
#' head(predict(fit, newdata = d, type = "conditional"))
#' head(predict(fit, newdata = d, type = "selection"))
#' @export
predict.panselect <- function(object, newdata,
                              type = c("unconditional", "conditional", "selection"), ...) {
  type <- match.arg(type)
  if (missing(newdata) || is.null(newdata)) {
    stop("panselect fits do not store the data: supply 'newdata' (for in-sample ",
         "predictions, the data used in the fit)", call. = FALSE)
  }
  newdata <- as.data.frame(newdata)
  missing_index <- setdiff(object$index, names(newdata))
  if (length(missing_index)) {
    stop("newdata must contain the panel index column(s) ",
         paste(missing_index, collapse = ", "), call. = FALSE)
  }
  b <- object$coefficients
  n <- nrow(newdata)
  id <- newdata[[object$index[1]]]
  tm <- newdata[[object$index[2]]]

  # selection regressors with the Mundlak means over each individual's rows
  needs_sel <- type != "unconditional" || (object$model == "pooling" && length(object$mean_cols))
  D <- NULL
  if (needs_sel) {
    sf <- stats::delete.response(stats::terms(object$selection))
    mf <- stats::model.frame(sf, data = newdata, na.action = stats::na.pass)
    Z <- stats::model.matrix(sf, data = mf)
    D <- Z
    if (length(object$mean_cols)) {
      miss <- setdiff(object$mean_cols, colnames(Z))
      if (length(miss)) {
        stop("newdata does not produce the selection column(s) ", paste(miss, collapse = ", "),
             call. = FALSE)
      }
      M <- vapply(object$mean_cols, function(v) stats::ave(Z[, v], id, FUN = mean), numeric(n))
      M <- matrix(M, nrow = n, dimnames = list(NULL, paste0("mean_", object$mean_cols)))
      D <- cbind(Z, M)
    }
  }

  if (type == "selection" || type == "conditional") {
    xb <- rep(NA_real_, n)
    has <- rep(FALSE, n)
    known <- as.character(tm) %in% names(object$probit)
    for (p in names(object$probit)) {
      d <- object$probit[[p]]
      r <- which(as.character(tm) == p)
      if (!length(r) || is.null(d)) next
      d <- d[!is.na(d)]
      miss <- setdiff(names(d), colnames(D))
      if (length(miss)) {
        stop("newdata does not produce the selection column(s) ", paste(miss, collapse = ", "),
             call. = FALSE)
      }
      xb[r] <- as.numeric(D[r, names(d), drop = FALSE] %*% d)
      has[r] <- TRUE
    }
    if (type == "selection") return(stats::pnorm(xb))
    imr <- ifelse(has, exp(stats::dnorm(xb, log = TRUE) - stats::pnorm(xb, log.p = TRUE)), 0)
    imr[!known] <- NA_real_
  }

  # outcome equation index
  of <- stats::delete.response(stats::terms(object$outcome))
  mf <- stats::model.frame(of, data = newdata, na.action = stats::na.pass)
  X <- stats::model.matrix(of, data = mf)
  if (object$model == "pooling" && length(object$mean_cols)) {
    X <- cbind(X, D[, paste0("mean_", object$mean_cols), drop = FALSE])
  }
  xterms <- setdiff(names(b), object$lambda_terms)
  miss <- setdiff(xterms, colnames(X))
  if (length(miss)) {
    stop("newdata does not produce the outcome column(s) ", paste(miss, collapse = ", "),
         " (check the factor levels)", call. = FALSE)
  }
  eta <- as.numeric(X[, xterms, drop = FALSE] %*% b[xterms])
  if (type == "unconditional") return(eta)

  g <- if (object$lambda == "common") {
    rep(unname(b["lambda"]), n)
  } else {
    unname(b[paste0("lambda:", as.character(tm))])
  }
  g[is.na(g)] <- 0
  eta + g * imr
}

# ---- pangmm -----------------------------------------------------------------

#' Predictions from a dynamic panel GMM fit
#'
#' Linear predictions from a [pangmm()] fit on `newdata`, a long data frame
#' with the panel index columns and the variables of the first part of the
#' formula. Lags are taken within individuals over the sorted periods of
#' `newdata`, as in the fit.
#'
#' * `type = "level"`: `x_it'b + d_t (+ c)`, the regressors (with their
#'   lags, including the lagged dependent variable) times their
#'   coefficients, plus the time effect and, for system GMM, the constant.
#'   The individual effect is not identified by GMM and is excluded, so for
#'   difference GMM (no constant) this is the level up to an individual
#'   intercept. On the fitting data these are the fitted values of the
#'   level equations of system GMM.
#' * `type = "difference"`: the prediction of the first difference
#'   `y_it - y_i,t-1`, the difference of the level predictions, in which
#'   the individual effect cancels. On the fitting data these are the
#'   fitted values of the differenced equations.
#'
#' Periods that precede the first time effect of the fit have a time effect
#' of zero (the reference); a later period that is not in the fit gives a
#' missing prediction, with a warning. A prediction is missing whenever a
#' required lag is unavailable.
#'
#' @param object a `"pangmm"` fit
#' @param newdata a data frame in long format; the fit does not store its
#'   data, so `newdata` is required (pass the fitting data for in-sample
#'   predictions)
#' @param type `"level"` or `"difference"`
#' @param ... unused
#' @return a numeric vector with one prediction per row of `newdata`
#' @family pangmm methods
#' @section Validation:
#' Checked against the transformed equations stored by `plm::pgmm()` in the
#' package's validation suite on GitHub.
#' @examples
#' set.seed(3)
#' N <- 100; TT <- 6
#' d <- data.frame(id = rep(seq_len(N), each = TT), t = rep(seq_len(TT), N))
#' d$x <- rnorm(N * TT)
#' d$y <- ave(rep(rnorm(N), each = TT) + 0.3 * d$x + rnorm(N * TT), d$id,
#'            FUN = function(e) stats::filter(e, 0.4, method = "recursive"))
#' fit <- pangmm(y ~ lag(y, 1) + x | lag(y, 2:99), data = d, index = c("id", "t"),
#'               transformation = "ld")
#' head(predict(fit, newdata = d))
#' head(predict(fit, newdata = d, type = "difference"))
#' @export
predict.pangmm <- function(object, newdata, type = c("level", "difference"), ...) {
  type <- match.arg(type)
  if (missing(newdata) || is.null(newdata)) {
    stop("pangmm fits do not store the data: supply 'newdata' (for in-sample ",
         "predictions, the data used in the fit)", call. = FALSE)
  }
  newdata <- as.data.frame(newdata)
  index <- object$index
  missing_index <- setdiff(index, names(newdata))
  if (length(missing_index)) {
    stop("newdata must contain the panel index column(s) ",
         paste(missing_index, collapse = ", "), call. = FALSE)
  }
  formula <- object$formula
  env <- environment(formula) %||% parent.frame()
  parts <- gmm_split_formula(formula)
  main <- gmm_terms(parts$main, merge = FALSE)
  id <- newdata[[index[1]]]
  tm <- newdata[[index[2]]]
  if (anyNA(id) || anyNA(tm)) stop("panel index columns must not contain missing values",
                                   call. = FALSE)
  uid <- sort(unique(id))
  utm <- sort(unique(tm))
  ii <- match(id, uid)
  tt <- match(tm, utm)
  if (anyDuplicated(cbind(ii, tt))) {
    stop("each individual-time index pair must identify a unique observation", call. = FALSE)
  }
  N <- length(uid)
  TT <- length(utm)
  grid <- function(e) {
    v <- eval(e, newdata, env)
    if (is.logical(v)) v <- as.numeric(v)
    if (!is.numeric(v) || length(v) != nrow(newdata)) {
      stop("'", gmm_deparse(e), "' must evaluate to a numeric vector with one value per row",
           call. = FALSE)
    }
    M <- matrix(NA_real_, N, TT)
    M[cbind(ii, tt)] <- as.numeric(v)
    M
  }
  lagged <- function(M, k) {
    if (k == 0) return(M)
    out <- matrix(NA_real_, N, TT)
    if (k < TT) out[, (k + 1):TT] <- M[, seq_len(TT - k), drop = FALSE]
    out
  }

  b <- object$coefficients
  eta <- matrix(0, N, TT)
  for (term in main) {
    lab <- gmm_deparse(term$expr)
    M <- grid(term$expr)
    for (k in term$lags) {
      nm <- gmm_lagname(lab, k)
      if (nm %in% names(b)) eta <- eta + b[[nm]] * lagged(M, k)
    }
  }

  # time effects: zero before the first dummy period of the fit
  dummies <- setdiff(object$timedummies, "(Intercept)")
  delta <- stats::setNames(numeric(TT), as.character(utm))
  if (length(dummies)) {
    pos <- match(as.character(utm), dummies)
    delta[!is.na(pos)] <- b[dummies[pos[!is.na(pos)]]]
    unmatched <- which(is.na(pos))
    if (length(unmatched)) {
      first <- dummies[1L]
      num <- suppressWarnings(as.numeric(c(first, as.character(utm[unmatched]))))
      earlier <- if (!anyNA(num)) num[-1L] < num[1L] else as.character(utm[unmatched]) < first
      if (any(!earlier)) {
        warning("period(s) ", paste(utm[unmatched][!earlier], collapse = ", "),
                " of newdata have no time effect in the fit: predictions are NA",
                call. = FALSE)
        delta[unmatched[!earlier]] <- NA_real_
      }
    }
  }
  const <- if ("(Intercept)" %in% names(b)) b[["(Intercept)"]] else 0
  level <- eta + rep(delta, each = N) + const
  out <- if (type == "level") {
    level
  } else {
    D <- matrix(NA_real_, N, TT)
    if (TT > 1L) D[, -1] <- level[, -1, drop = FALSE] - level[, -TT, drop = FALSE]
    D
  }
  as.numeric(out[cbind(ii, tt)])
}
