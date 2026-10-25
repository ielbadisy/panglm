#' Two-part model for panel cost data
#'
#' Fits the two-part model of Duan et al. (1983), the standard model for
#' health care expenditures with a mass of zeros: a binary model for
#' `P(y > 0)` on all observations, and a model for `E[y | y > 0]` on the
#' positive observations. The unconditional mean is
#' `E[y] = P(y > 0) * E[y | y > 0]`.
#'
#' The positive part is either a gamma GLM with a log link
#' (`positive = "gamma"`, no retransformation needed) or a linear model for
#' `log(y)` (`positive = "lognormal"`). For the log model,
#' `E[y | y > 0] = exp(x'beta) * phi` needs a retransformation factor
#' `phi`: Duan's smearing estimate, the (weighted) mean of the exponentiated
#' log-scale residuals (`retransform = "smearing"`), or the normal-theory
#' factor `exp(sigma^2 / 2)` (`retransform = "normal"`; with
#' `model = "random"`, `sigma^2` is the sum of the individual and
#' idiosyncratic variance components). Both assume homoskedastic log-scale
#' errors.
#'
#' With `model = "random"` each part has its own random intercept and the
#' parts are estimated separately; the combined mean multiplies the two
#' population-averaged means, which assumes the two random intercepts are
#' independent. With `correlated = TRUE` (random effects, lognormal positive
#' part) the two random intercepts are instead bivariate normal with
#' correlation `rho` and both parts are estimated jointly by maximum
#' likelihood (Van Ourti 2004): individuals with a high propensity to use
#' care may also have high spending when they do. The positive-part random
#' intercept is integrated out in closed form and the binary-part one by
#' adaptive Gauss-Hermite quadrature. Predictions then integrate over the
#' joint distribution (normal-theory retransformation). A likelihood-ratio
#' test of `rho = 0` compares `logLik()` with the separate fit.
#' `cre = TRUE` adds the individual means of the time-varying
#' regressors (computed over all observations of the individual, including
#' zeros) to both parts, the Mundlak correlated random effects device.
#'
#' @param formula model formula for the positive part, e.g.
#'   `cost ~ age + treatment`; also used for the binary part unless
#'   `binary_formula` is given
#' @param data a data.frame
#' @param index length-1 or length-2 character vector of panel index columns
#' @param binary_formula optional right-hand-side formula for the binary
#'   part, e.g. `~ age + treatment + smoker`
#' @param model `"pooling"` or `"random"`, used for both parts
#' @param link link of the binary part, `"probit"` or `"logit"`
#' @param positive `"gamma"` (gamma GLM, log link) or `"lognormal"`
#'   (linear model for `log(y)`)
#' @param retransform for `positive = "lognormal"`: `"smearing"` or
#'   `"normal"`
#' @param correlated logical; with `model = "random"` and
#'   `positive = "lognormal"`, let the random intercepts of the two parts be
#'   correlated and fit them jointly
#' @param cre logical; add individual means of the time-varying regressors
#'   to both parts
#' @param weights optional observation weights (`model = "pooling"` only)
#' @param vcov covariance type passed to [panglm()] for both parts
#' @param R,maxit,tol passed to [panglm()]
#' @return an object of class `"twopart"` with components `binary` and
#'   `positive` (the two [panglm()] fits), `retransformation` (the factor
#'   `phi`, 1 for the gamma part) and the settings used. Methods: `print`,
#'   `summary`, `coef`, `vcov`, `predict`, `fitted`, `logLik`, `nobs` and
#'   [partialeffects()].
#' @references Duan, N., Manning, W. G., Morris, C. N. and Newhouse, J. P.
#'   (1983). A comparison of alternative models for the demand for medical
#'   care. *Journal of Business and Economic Statistics*, 1, 115-126.
#'   \doi{10.1080/07350015.1983.10509330}
#'
#'   Duan, N. (1983). Smearing estimate: a nonparametric retransformation
#'   method. *Journal of the American Statistical Association*, 78,
#'   605-610. \doi{10.1080/01621459.1983.10478017}
#'
#'   Van Ourti, T. (2004). Measuring horizontal inequity in Belgian health
#'   care using a Gaussian random effects two part count data model.
#'   *Health Economics*, 13, 705-724. \doi{10.1002/hec.869}
#' @examples
#' data(copd)
#' set.seed(1)
#' copd$cost <- ifelse(copd$exacerbations > 0,
#'                     exp(6 + 0.1 * copd$crp + rnorm(nrow(copd), 0, 0.5)), 0)
#' tp <- twopart(cost ~ treatment + age + crp, data = copd,
#'               index = c("id", "visit"), positive = "gamma")
#' summary(tp)
#' head(predict(tp))
#' partialeffects(tp)
#' @export
twopart <- function(formula, data, index, binary_formula = NULL,
                    model = c("pooling", "random"),
                    link = c("probit", "logit"),
                    positive = c("gamma", "lognormal"),
                    retransform = c("smearing", "normal"),
                    correlated = FALSE, cre = FALSE, weights = NULL,
                    vcov = c("classical", "HC1", "cluster"),
                    R = 21, maxit = 100, tol = 1e-10) {
  model <- match.arg(model)
  link <- match.arg(link)
  positive <- match.arg(positive)
  retransform <- match.arg(retransform)
  vcov <- match.arg(vcov)
  if (missing(index)) stop("'index' is required, e.g. index = c(\"id\", \"time\")", call. = FALSE)
  missing_index <- setdiff(index, names(data))
  if (length(missing_index)) {
    stop("index column(s) not found in data: ", paste(missing_index, collapse = ", "),
         call. = FALSE)
  }
  if (correlated && !(model == "random" && positive == "lognormal")) {
    stop("correlated = TRUE needs model = 'random' and positive = 'lognormal'", call. = FALSE)
  }
  if (correlated && vcov == "HC1") {
    stop("vcov = 'HC1' is not available for random effects: use 'cluster'", call. = FALSE)
  }
  if (!is.null(weights) && (!is.numeric(weights) || length(weights) != nrow(data))) {
    stop("'weights' must be a numeric vector with nrow(data) values", call. = FALSE)
  }

  formula <- stats::as.formula(formula)
  rhs_binary <- if (is.null(binary_formula)) {
    stats::delete.response(stats::terms(formula))
  } else {
    stats::delete.response(stats::terms(stats::as.formula(binary_formula)))
  }
  binary_formula <- stats::reformulate(attr(rhs_binary, "term.labels"),
                                       intercept = attr(rhs_binary, "intercept") == 1)

  vars <- unique(c(all.vars(formula), all.vars(binary_formula), index))
  vars <- intersect(vars, names(data))
  keep <- stats::complete.cases(data[, vars, drop = FALSE])
  d <- data[keep, vars, drop = FALSE]
  w <- if (is.null(weights)) NULL else weights[keep]

  y <- stats::model.response(stats::model.frame(formula, d))
  if (!is.numeric(y) || any(!is.finite(y)) || any(y < 0)) {
    stop("the response must be finite, numeric and nonnegative", call. = FALSE)
  }
  if (all(y > 0) || all(y == 0)) {
    stop("a two-part model needs both zero and positive outcomes", call. = FALSE)
  }

  rhs_positive <- stats::reformulate(attr(stats::terms(formula), "term.labels"),
                                     intercept = attr(stats::terms(formula), "intercept") == 1)
  cre_columns <- character(0)
  cre_spec <- NULL
  if (cre) {
    cre_spec <- list(binary = twopart_cre_spec(binary_formula, d, index[1]),
                     positive = twopart_cre_spec(rhs_positive, d, index[1]))
    d <- twopart_add_means(d, cre_spec, index[1])
    cre_columns <- unique(c(cre_spec$binary$names, cre_spec$positive$names))
    binary_formula <- twopart_extend(binary_formula, cre_spec$binary$names)
    rhs_positive <- twopart_extend(rhs_positive, cre_spec$positive$names)
  }

  d$.twopart_any <- as.numeric(y > 0)
  fit_binary <- panglm(
    stats::update(binary_formula, .twopart_any ~ .), data = d, index = index,
    model = model, family = binomialfamily(link), weights = w, vcov = vcov,
    R = R, maxit = maxit, tol = tol
  )

  is_pos <- y > 0
  dp <- d[is_pos, , drop = FALSE]
  wp <- if (is.null(w)) NULL else w[is_pos]
  dp$.twopart_y <- if (positive == "gamma") y[is_pos] else log(y[is_pos])
  fit_positive <- panglm(
    stats::update(rhs_positive, .twopart_y ~ .), data = dp, index = index,
    model = model, family = if (positive == "gamma") gammafamily() else gaussianfamily(),
    weights = wp, vcov = vcov, R = R, maxit = maxit, tol = tol
  )

  phi <- 1
  if (positive == "lognormal") {
    phi <- if (retransform == "smearing") {
      resid <- fit_positive$y - as.numeric(fit_positive$X %*% fit_positive$coefficients)
      pw <- fit_positive$weights %||% rep(1, length(resid))
      sum(pw * exp(resid)) / sum(pw)
    } else if (model == "random") {
      exp((fit_positive$sigma_v2 + fit_positive$sigma_mu2) / 2)
    } else {
      exp(fit_positive$dispersion / 2)
    }
  }

  joint <- NULL
  if (correlated) {
    joint <- fit_twopart_cre(fit_binary, fit_positive, index, link, vcov, R, maxit, tol)
    retransform <- "normal"
    kb <- joint$kb; kp <- joint$kp
    ib <- seq_len(kb); ip <- kb + seq_len(kp)
    fit_binary$coefficients[] <- joint$theta[ib]
    fit_positive$coefficients[] <- joint$theta[ip]
    fit_binary$vcov[] <- joint$vcov[ib, ib]
    fit_positive$vcov[] <- joint$vcov[ip, ip]
    fit_binary$sigma <- joint$sigma_binary
    fit_positive$sigma_mu2 <- joint$sigma_positive^2
    fit_positive$sigma_v2 <- joint$sigma_e^2
    fit_binary$vcov_type <- fit_positive$vcov_type <- vcov
    phi <- NA_real_
  }

  out <- list(binary = fit_binary, positive = fit_positive, retransformation = phi,
              correlated = correlated, joint = joint,
              positive_type = positive, retransform = if (positive == "lognormal") retransform,
              link = link, model = model, formula = formula, binary_formula = binary_formula,
              index = index, data = d[, setdiff(names(d), ".twopart_any"), drop = FALSE],
              weights = w, y = y, cre = cre, cre_columns = cre_columns, cre_spec = cre_spec,
              nobs = length(y), n_positive = sum(is_pos), call = match.call())
  class(out) <- "twopart"
  out
}

twopart_cre_spec <- function(rhs, data, id_name) {
  terms_object <- stats::terms(rhs)
  X <- stats::model.matrix(terms_object, data)
  X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  xlev <- stats::.getXlevels(terms_object, stats::model.frame(terms_object, data))
  if (!ncol(X)) return(list(terms = terms_object, xlevels = xlev, columns = character(0),
                            names = character(0)))
  id <- data[[id_name]]
  varying <- vapply(seq_len(ncol(X)), function(j) {
    any(abs(X[, j] - stats::ave(X[, j], id)) > 1e-8 * max(1, abs(X[, j])))
  }, logical(1))
  columns <- colnames(X)[varying]
  list(terms = terms_object, xlevels = xlev, columns = columns,
       names = if (length(columns)) paste0("mean_", make.names(columns)) else character(0))
}

twopart_add_means <- function(data, cre_spec, id_name) {
  id <- data[[id_name]]
  for (part in cre_spec) {
    if (!length(part$columns)) next
    X <- stats::model.matrix(part$terms,
                             stats::model.frame(part$terms, data, xlev = part$xlevels,
                                                na.action = stats::na.pass))
    for (j in seq_along(part$columns)) {
      data[[part$names[j]]] <- stats::ave(X[, part$columns[j]], id)
    }
  }
  data
}

twopart_extend <- function(rhs, names) {
  if (!length(names)) return(rhs)
  terms_object <- stats::terms(rhs)
  stats::reformulate(c(attr(terms_object, "term.labels"), names),
                     intercept = attr(terms_object, "intercept") == 1)
}

twopart_predictor <- function(object, part = c("total", "binary", "positive")) {
  part <- match.arg(part)
  k_binary <- length(object$binary$coefficients)
  index_binary <- seq_len(k_binary)
  index_positive <- k_binary + seq_along(object$positive$coefficients)
  if (isTRUE(object$correlated)) {
    return(function(theta, data) {
      a <- as.numeric(panglm_design(object$binary, data) %*% theta[index_binary])
      b <- as.numeric(panglm_design(object$positive, data) %*% theta[index_positive])
      m <- twopart_cre_means(a, b, object$joint, object$link)
      switch(part, total = m$total, binary = m$prob, positive = m$total / m$prob)
    })
  }
  binary_mean <- panglm_mean_function(object$binary)
  positive_mean <- if (object$positive_type == "gamma") {
    panglm_mean_function(object$positive)
  } else {
    phi <- object$retransformation
    function(eta) phi * exp(eta)
  }
  function(theta, data) {
    p <- if (part != "positive") {
      binary_mean(as.numeric(panglm_design(object$binary, data) %*% theta[index_binary]))
    } else 1
    m <- if (part != "binary") {
      positive_mean(as.numeric(panglm_design(object$positive, data) %*% theta[index_positive]))
    } else 1
    p * m
  }
}

#' Predict from a two-part model
#'
#' @param object a `"twopart"` fit
#' @param newdata optional data frame; defaults to the (complete-case) data
#'   used for fitting. For `cre = TRUE` fits, the individual means are
#'   recomputed within `newdata` unless its `mean_*` columns are supplied.
#' @param type `"response"` for `E[y] = P(y > 0) E[y | y > 0]`, `"prob"`
#'   for `P(y > 0)`, or `"positive"` for `E[y | y > 0]` (retransformed for
#'   the lognormal part)
#' @param ... unused
#' @return a numeric vector with one prediction per row of `newdata`
#' @export
predict.twopart <- function(object, newdata = NULL,
                            type = c("response", "prob", "positive"), ...) {
  type <- match.arg(type)
  if (is.null(newdata)) {
    newdata <- object$data
  } else if (object$cre && !all(object$cre_columns %in% names(newdata))) {
    if (!object$index[1] %in% names(newdata)) {
      stop("newdata must contain panel index column '", object$index[1], "'", call. = FALSE)
    }
    newdata <- twopart_add_means(newdata, object$cre_spec, object$index[1])
  }
  part <- switch(type, response = "total", prob = "binary", positive = "positive")
  twopart_predictor(object, part)(stats::coef(object), newdata)
}

#' @export
fitted.twopart <- function(object, ...) predict.twopart(object)

#' @export
coef.twopart <- function(object, ...) {
  b <- object$binary$coefficients
  p <- object$positive$coefficients
  stats::setNames(c(b, p), c(paste0("binary:", names(b)), paste0("positive:", names(p))))
}

#' @export
vcov.twopart <- function(object, ...) {
  if (isTRUE(object$correlated)) {
    v <- object$joint$vcov
    dimnames(v) <- list(names(stats::coef(object)), names(stats::coef(object)))
    return(v)
  }
  vb <- object$binary$vcov
  vp <- object$positive$vcov
  kb <- nrow(vb); kp <- nrow(vp)
  out <- matrix(0, kb + kp, kb + kp)
  out[seq_len(kb), seq_len(kb)] <- vb
  out[kb + seq_len(kp), kb + seq_len(kp)] <- vp
  dimnames(out) <- list(names(stats::coef(object)), names(stats::coef(object)))
  out
}

#' @export
logLik.twopart <- function(object, ...) {
  if (isTRUE(object$correlated)) {
    return(structure(object$joint$loglik, df = object$joint$npar, nobs = object$nobs,
                     class = "logLik"))
  }
  lb <- stats::logLik(object$binary)
  lp <- stats::logLik(object$positive)
  value <- as.numeric(lb) + as.numeric(lp)
  if (object$positive_type == "lognormal") {
    # Jacobian of the log transformation, so the likelihood is on the cost scale
    pw <- object$positive$weights %||% rep(1, length(object$positive$y))
    value <- value - sum(pw * object$positive$y)
  }
  attr(value, "df") <- attr(lb, "df") + attr(lp, "df")
  attr(value, "nobs") <- object$nobs
  class(value) <- "logLik"
  value
}

#' @export
nobs.twopart <- function(object, ...) object$nobs

#' @export
print.twopart <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nTwo-part model (", x$model, if (isTRUE(x$correlated)) ", correlated", "): ",
      x$link, " for P(y > 0); ",
      if (x$positive_type == "gamma") "gamma (log link)" else "log-linear",
      " for E[y | y > 0]\n", sep = "")
  cat("N =", x$nobs, " positive =", x$n_positive, " groups =", x$binary$n_groups, "\n\n")
  print(round(stats::coef(x), digits))
  cat("\n")
  invisible(x)
}

#' @export
summary.twopart <- function(object, ...) {
  part_test <- function(fit) {
    cre_wald_test(list(coefficients = fit$coefficients, vcov = fit$vcov,
                       cre_columns = intersect(object$cre_columns, names(fit$coefficients))))
  }
  out <- list(object = object, binary = summary(object$binary),
              positive = summary(object$positive),
              cre_tests = list(binary = part_test(object$binary),
                               positive = part_test(object$positive)))
  class(out) <- "summary.twopart"
  out
}

#' @export
print.summary.twopart <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  object <- x$object
  cat("\nTwo-part model (", object$model, if (isTRUE(object$correlated)) ", correlated",
      ")\n", sep = "")
  cat("N =", object$nobs, " positive =", object$n_positive,
      " groups =", object$binary$n_groups, "\n")
  cat("\nPart 1: P(y > 0), binomial", object$link, "\n")
  stats::printCoefmat(x$binary$coefficients, digits = digits, has.Pvalue = TRUE)
  cat("\nPart 2: E[y | y > 0],",
      if (object$positive_type == "gamma") "gamma GLM with log link" else "linear model for log(y)",
      "\n")
  stats::printCoefmat(x$positive$coefficients, digits = digits, has.Pvalue = TRUE)
  if (isTRUE(object$correlated)) {
    vp <- object$joint$variance_params
    se <- sqrt(diag(object$joint$variance_vcov))
    cat("\nCorrelated random intercepts (joint maximum likelihood)\n")
    print(round(cbind(Estimate = vp, `Std. Error` = se), digits))
  } else if (object$positive_type == "lognormal") {
    cat("\nRetransformation (", object$retransform, "): ",
        format(object$retransformation, digits = digits), "\n", sep = "")
  }
  for (part in c("binary", "positive")) {
    test <- x$cre_tests[[part]]
    if (!is.null(test)) {
      cat("Mundlak test, ", part, " part: chi2(", test$df, ") = ",
          format(test$statistic, digits = digits), ", p = ",
          format.pval(test$p.value, digits = digits), "\n", sep = "")
    }
  }
  cat("Log-likelihood:", format(as.numeric(stats::logLik(object)), digits = digits), "\n\n")
  invisible(x)
}
