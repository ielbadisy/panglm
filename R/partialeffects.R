#' Average partial effects
#'
#' Average partial effects (APE, also called average marginal effects) of
#' each regressor on the expected outcome, with delta-method standard
#' errors. For a continuous regressor the APE is the sample average of the
#' derivative of the mean with respect to that regressor; for a binary
#' (0/1 or logical) regressor and for each non-reference level of a factor
#' it is the sample average of the discrete change in the mean, all other
#' regressors held at their observed values. Derivatives and changes are
#' taken through the model formula, so interactions, polynomial terms and
#' factor codings are handled automatically.
#'
#' The mean used is the one returned by `predict(type = "response")`:
#' the GLM mean for pooled models; the mean including the estimated
#' individual (and time) effects for fixed-effects models; and the
#' population-averaged mean (integrating over the random intercept, and
#' the random slopes when `random` has slope variables) for
#' random-effects models, which is the average structural function of
#' Wooldridge (2005). For correlated random effects fits (`cre = TRUE`) the
#' individual means are held fixed, so the APE of a time-varying regressor
#' is the within-individual effect. For two-part models the mean is
#' `P(y > 0) * E[y | y > 0]` (or one of its parts, see `part`).
#'
#' Standard errors use the covariance stored on the fit (so a fit made with
#' `vcov = "cluster"` gives cluster-robust APE standard errors) and treat
#' variance components, dispersion and fixed-effect intercepts as fixed.
#' Exact conditional logit fits are not supported, because their
#' individual intercepts are conditioned out and the mean is not identified.
#'
#' @param object a `"panglm"` or `"twopart"` fit
#' @param variables character vector of variable names (as they appear in
#'   the data); defaults to every regressor in the model formula(s)
#' @param level confidence level for the Wald intervals
#' @param part for two-part models: `"total"` (the unconditional mean),
#'   `"binary"` (the probability of a positive outcome) or `"positive"` (the
#'   mean conditional on a positive outcome)
#' @param ... unused
#' @return a data frame of class `"panglm_ape"` with columns `term`,
#'   `contrast`, `estimate`, `std.error`, `statistic`, `p.value`,
#'   `conf.low` and `conf.high`, one row per continuous or binary regressor
#'   and one row per non-reference factor level.
#' @references Wooldridge, J. M. (2005). Unobserved heterogeneity and
#'   estimation of average partial effects. In *Identification and
#'   Inference for Econometric Models*, Cambridge University Press.
#' @family panglm methods
#' @seealso [effectcurve()] for the effect of a spline term across its range.
#' @section Validation:
#' Numerical accuracy is checked in the package tests and in the validation
#' suite; see the validation article on the package website.
#' @examples
#' data(copd)
#' fit <- panglm(hospitalized ~ treatment + age + crp, data = copd,
#'               index = c("id", "visit"), model = "random",
#'               family = "probit", cre = TRUE)
#' partialeffects(fit)
#' @export
partialeffects <- function(object, ...) UseMethod("partialeffects")

#' @rdname partialeffects
#' @export
partialeffects.panglm <- function(object, variables = NULL, level = 0.95, ...) {
  if (is.null(object$data)) {
    stop("this fit does not store its data; re-fit it with the current panglm", call. = FALSE)
  }
  predictor <- panglm_predictor(object)
  variables <- ape_variables(object$terms, object$data, object$index, variables)
  out <- ape_engine(predictor, object$coefficients, object$vcov, object$data,
                    variables, object$weights, level, object$xlevels)
  attr(out, "vcov_type") <- object$vcov_type
  out
}

#' @rdname partialeffects
#' @export
partialeffects.twopart <- function(object, variables = NULL, level = 0.95,
                                   part = c("total", "binary", "positive"), ...) {
  part <- match.arg(part)
  predictor <- twopart_predictor(object, part)
  terms_union <- stats::terms(stats::as.formula(paste(
    "~", paste(c(attr(stats::terms(object$formula), "term.labels"),
                 attr(stats::terms(object$binary_formula), "term.labels")), collapse = " + ")
  )))
  variables <- ape_variables(terms_union, object$data, object$index, variables)
  variables <- setdiff(variables, object$cre_columns)
  xlev <- c(object$binary$xlevels, object$positive$xlevels)
  xlev <- xlev[!duplicated(names(xlev))]
  out <- ape_engine(predictor, stats::coef(object), stats::vcov(object), object$data,
                    variables, object$weights, level, xlev)
  attr(out, "vcov_type") <- object$binary$vcov_type
  attr(out, "part") <- part
  out
}

#' @export
print.panglm_ape <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  part <- attr(x, "part")
  cat("\nAverage partial effects",
      if (!is.null(part)) paste0(" (", part, " mean)") else "", "\n", sep = "")
  cat("Delta-method standard errors, vcov:", attr(x, "vcov_type") %||% "classical", "\n\n")
  table <- as.data.frame(x)
  rownames(table) <- NULL
  numeric_columns <- vapply(table, is.numeric, logical(1))
  table[numeric_columns] <- lapply(table[numeric_columns], signif, digits = digits)
  print(table, row.names = FALSE)
  cat("\n")
  invisible(x)
}

panglm_mean_function <- function(object) {
  family <- object$family
  if (family$family == "tobit") {
    s <- tobit_total_sd(object)
    return(function(eta) tobit_censored_mean(eta, s, family$left, family$right))
  }
  if (object$model != "random") return(function(eta) linkinv_r(eta, family$link_id))
  switch(family$family,
    gaussian = function(eta) eta,
    poisson = function(eta) exp(eta),
    negbin = {
      multiplier <- object$marginal_mean_multiplier
      function(eta) multiplier * exp(eta)
    },
    binomial = ,
    beta = {
      gh <- marginal_quadrature(object$quadrature_nodes %||% 21L)
      function(eta) random_binomial_marginal_mean(eta, object$sigma, family, gh)
    },
    gamma = {
      sigma <- object$sigma
      function(eta) exp(eta + sigma^2 / 2)
    }
  )
}

#' Mean prediction as a function of (coefficients, data)
#'
#' For the fitted data only (fixed-effect shifts and held CRE means are
#' aligned with `object$data`).
#'
#' @keywords internal
#' @noRd
panglm_predictor <- function(object, hold_cre = TRUE) {
  if (object$model == "within" && object$family$family == "binomial") {
    stop("partial effects are not identified for exact conditional logit fits; ",
         "use model = 'random' with cre = TRUE instead", call. = FALSE)
  }
  shift <- 0
  if (object$model == "within") {
    shift <- as.numeric(object$linear.predictors - object$X %*% object$coefficients)
  }
  cre_columns <- object$cre_columns %||% character(0)
  cre_means <- if (hold_cre && length(cre_columns)) object$X[, cre_columns, drop = FALSE] else NULL
  if (!is.null(object$varcomp)) {
    # random slopes: the population-averaged mean depends on z_it
    return(function(theta, data) {
      X <- panglm_design(object, data, cre_means)
      slopes_marginal_mean(object$family, as.numeric(X %*% theta),
                           slopes_z(object$random, data), object$varcomp)
    })
  }
  mean_fun <- panglm_mean_function(object)
  function(theta, data) {
    X <- panglm_design(object, data, cre_means)
    mean_fun(as.numeric(X %*% theta) + shift)
  }
}

ape_variables <- function(terms_object, data, index, variables) {
  available <- setdiff(intersect(all.vars(stats::delete.response(terms_object)), names(data)), index)
  if (is.null(variables)) return(available)
  unknown <- setdiff(variables, available)
  if (length(unknown)) {
    stop("variable(s) not among the model regressors: ", paste(unknown, collapse = ", "),
         call. = FALSE)
  }
  variables
}

ape_engine <- function(predictor, theta, V, data, variables, weights, level, xlevels) {
  base <- predictor(theta, data)
  w <- weights %||% rep(1, nrow(data))
  usable <- is.finite(base) & w > 0
  w <- w[usable] / sum(w[usable])
  wmean <- function(values) sum(w * values[usable])
  crit <- stats::qnorm(1 - (1 - level) / 2)

  # each contrast is a function of theta: average of (pred(d1) - pred(d0)) / denom
  contrast_rows <- list()
  for (v in variables) {
    x <- data[[v]]
    discrete <- is.factor(x) || is.character(x) || is.logical(x) ||
      (is.numeric(x) && all(x %in% c(0, 1)))
    if (discrete) {
      if (is.logical(x)) {
        settings <- list(list(label = "TRUE - FALSE", d1 = TRUE, d0 = FALSE))
      } else if (is.numeric(x)) {
        settings <- list(list(label = "1 - 0", d1 = 1, d0 = 0))
      } else {
        lev <- xlevels[[v]] %||% levels(factor(x))
        settings <- lapply(lev[-1], function(l) {
          list(label = paste(l, "-", lev[1]), d1 = l, d0 = lev[1])
        })
      }
      for (s in settings) {
        d1 <- data; d0 <- data
        d1[[v]] <- ape_set(x, s$d1)
        d0[[v]] <- ape_set(x, s$d0)
        contrast_rows[[length(contrast_rows) + 1L]] <-
          list(term = v, contrast = s$label, d1 = d1, d0 = d0, denom = 1)
      }
    } else {
      step <- 1e-4 * stats::sd(x)
      if (!is.finite(step) || step == 0) next
      d1 <- data; d0 <- data
      d1[[v]] <- x + step
      d0[[v]] <- x - step
      contrast_rows[[length(contrast_rows) + 1L]] <-
        list(term = v, contrast = "dY/dX", d1 = d1, d0 = d0, denom = 2 * step)
    }
  }

  results <- lapply(contrast_rows, function(r) {
    f <- function(th) wmean(predictor(th, r$d1) - predictor(th, r$d0)) / r$denom
    estimate <- f(theta)
    gradient <- vapply(seq_along(theta), function(j) {
      d <- 1e-4 * max(1, abs(theta[j]))
      up <- theta; up[j] <- up[j] + d
      dn <- theta; dn[j] <- dn[j] - d
      (f(up) - f(dn)) / (2 * d)
    }, numeric(1))
    if (estimate == 0 && all(gradient == 0)) return(NULL)
    se <- sqrt(as.numeric(t(gradient) %*% V %*% gradient))
    data.frame(term = r$term, contrast = r$contrast, estimate = estimate,
               std.error = se, statistic = estimate / se,
               p.value = 2 * stats::pnorm(-abs(estimate / se)),
               conf.low = estimate - crit * se, conf.high = estimate + crit * se,
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, results)
  if (is.null(out)) stop("no regressor with an identified partial effect", call. = FALSE)
  rownames(out) <- NULL
  class(out) <- c("panglm_ape", "data.frame")
  out
}

ape_set <- function(x, value) {
  n <- length(x)
  if (is.factor(x)) return(factor(rep(value, n), levels = levels(x)))
  rep(value, n)
}
