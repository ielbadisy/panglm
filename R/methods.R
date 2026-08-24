#' @export
coef.panglm <- function(object, ...) object$coefficients

#' Extract fitted values from a panglm model
#'
#' Values are returned in the row order of the data supplied to [panglm()],
#' even though estimation internally sorts observations by panel index.
#' Fixed-effects binomial values are conditional inclusion probabilities,
#' conditional on each individual's observed number of successes.
#'
#' @param object a fitted `"panglm"` object
#' @param ... unused
#' @export
fitted.panglm <- function(object, ...) {
  if (is.null(object$fitted.values)) return(NULL)
  object$fitted.values[object$inverse_order %||% seq_along(object$fitted.values)]
}

#' @export
residuals.panglm <- function(object, ...) {
  if (is.null(object$fitted.values) || is.null(object$y)) {
    stop("residuals are not available for this fitted model", call. = FALSE)
  }
  values <- object$y - object$fitted.values
  values[object$inverse_order %||% seq_along(values)]
}

#' @export
logLik.panglm <- function(object, ...) {
  val <- object$loglik
  attr(val, "df") <- length(object$coefficients)
  attr(val, "nobs") <- object$nobs
  class(val) <- "logLik"
  val
}

#' @export
nobs.panglm <- function(object, ...) object$nobs

#' Confidence intervals for panglm coefficients
#'
#' Wald-type intervals `estimate +/- z * se`, using whichever `vcov()` the
#' object currently carries (classical by default; HC1/cluster if the fit
#' or a prior [vcov.panglm()] call set it).
#'
#' @param object a `"panglm"` fit
#' @param parm which parameters (names or indices); defaults to all
#' @param level confidence level
#' @param ... unused
#' @examples
#' data(copd)
#' fit <- panglm(fev1 ~ crp, data = copd, index = c("id", "visit"),
#'               model = "within", family = "gaussian")
#' confint(fit)
#' @export
confint.panglm <- function(object, parm, level = 0.95, ...) {
  coefs <- object$coefficients
  if (missing(parm)) parm <- names(coefs)
  se <- sqrt(diag(object$vcov))[parm]
  crit <- stats::qnorm(1 - (1 - level) / 2)
  lo <- coefs[parm] - crit * se
  hi <- coefs[parm] + crit * se
  out <- cbind(lo, hi)
  pct <- c((1 - level) / 2, 1 - (1 - level) / 2) * 100
  colnames(out) <- sprintf("%.1f %%", pct)
  rownames(out) <- parm
  out
}

#' Predict from a fitted panglm model
#'
#' Pooled predictions use the ordinary GLM mean. Random-effects response
#' predictions are population averages over the fitted latent-effect
#' distribution. Fixed-effects predictions include the estimated individual
#' and, where applicable, time effects. New fixed-effect levels produce an
#' error unless `allow.new.levels = TRUE`, in which case their contribution
#' is set to zero.
#'
#' Exact conditional binomial fits provide in-sample conditional inclusion
#' probabilities. They cannot provide predictions for arbitrary `newdata`
#' because the individual intercept is conditioned out and the probability
#' depends on the complete response set and its success total.
#'
#' @param object a fitted `"panglm"` object
#' @param newdata optional data frame; when omitted, return in-sample values
#' @param type either `"link"` or `"response"`
#' @param allow.new.levels if `TRUE`, assign zero to unseen fixed-effect
#'   levels; otherwise unseen levels are an error
#' @param ... unused
#' @return a numeric vector of predictions
#' @examples
#' data(copd)
#' fit <- panglm(fev1 ~ crp, data = copd, index = c("id", "visit"),
#'               model = "random", family = "gaussian")
#' predict(fit, newdata = copd[1:4, ], type = "response")
#' @export
predict.panglm <- function(object, newdata = NULL,
                           type = c("link", "response"),
                           allow.new.levels = FALSE, ...) {
  type <- match.arg(type)
  if (is.null(newdata)) {
    values <- if (type == "response") object$fitted.values else object$linear.predictors
    if (is.null(values)) stop("in-sample predictions are unavailable for this fit", call. = FALSE)
    return(values[object$inverse_order %||% seq_along(values)])
  }

  if (object$model == "within" && object$family$family == "binomial") {
    stop("new-data prediction is not defined for exact conditional binomial fits; ",
         "use fitted() for probabilities conditional on each observed group total",
         call. = FALSE)
  }

  terms_object <- stats::delete.response(object$terms %||% stats::terms(object$formula))
  model_frame <- stats::model.frame(
    terms_object, data = newdata, xlev = object$xlevels, na.action = stats::na.pass
  )
  X <- stats::model.matrix(
    terms_object, data = model_frame, contrasts.arg = object$contrasts
  )
  coefficient_names <- names(object$coefficients)
  missing_columns <- setdiff(coefficient_names, colnames(X))
  if (length(missing_columns)) {
    stop("newdata does not produce required model column(s): ",
         paste(missing_columns, collapse = ", "), call. = FALSE)
  }
  eta <- as.numeric(X[, coefficient_names, drop = FALSE] %*% object$coefficients)

  if (object$model == "within") {
    id_name <- object$index[1]
    if (!id_name %in% names(newdata)) {
      stop("newdata must contain panel index column '", id_name, "'", call. = FALSE)
    }
    id_match <- match(as.character(newdata[[id_name]]), names(object$individual_effects))
    individual <- unname(object$individual_effects[id_match])
    if (anyNA(individual) && !allow.new.levels) {
      stop("newdata contains individual levels not present in the fitted model", call. = FALSE)
    }
    individual[is.na(individual)] <- 0
    eta <- eta + individual

    if (object$effect == "twoways") {
      time_name <- object$index[2]
      if (is.null(time_name) || !time_name %in% names(newdata)) {
        stop("newdata must contain the fitted time index column", call. = FALSE)
      }
      time_match <- match(as.character(newdata[[time_name]]), names(object$time_effects))
      time_effect <- unname(object$time_effects[time_match])
      if (anyNA(time_effect) && !allow.new.levels) {
        stop("newdata contains time levels not present in the fitted model", call. = FALSE)
      }
      time_effect[is.na(time_effect)] <- 0
      eta <- eta + time_effect
    }
  }

  if (type == "link") return(eta)
  if (object$model != "random") return(linkinv_r(eta, object$family$link_id))

  switch(object$family$family,
    gaussian = eta,
    poisson = exp(eta),
    negbin = {
      multiplier <- object$marginal_mean_multiplier
      if (!is.finite(multiplier)) {
        stop("the fitted beta-negative-binomial model has no finite marginal mean ",
             "because parameter 'a' is not greater than one", call. = FALSE)
      }
      multiplier * exp(eta)
    },
    binomial = random_binomial_marginal_mean(
      eta, object$sigma, object$family, object$quadrature_nodes %||% 21L
    )
  )
}

`%||%` <- function(x, y) if (is.null(x)) y else x
