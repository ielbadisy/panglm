#' @export
coef.panglm <- function(object, ...) object$coefficients

#' @export
fitted.panglm <- function(object, ...) object$fitted.values

#' @export
residuals.panglm <- function(object, ...) {
  if (is.null(object$fitted.values) || is.null(object$y)) {
    stop("residuals are not available for this model (no stored fitted values); ",
         "available for model = 'pooling' or 'within'/'random' with family = 'gaussian'",
         call. = FALSE)
  }
  object$y - object$fitted.values
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

#' @export
predict.panglm <- function(object, newdata = NULL, type = c("link", "response"), ...) {
  type <- match.arg(type)
  if (is.null(newdata)) {
    stop("predict.panglm currently requires 'newdata'", call. = FALSE)
  }
  X <- stats::model.matrix(stats::delete.response(stats::terms(object$formula)), data = newdata)
  common <- intersect(colnames(X), names(object$coefficients))
  eta <- as.numeric(X[, common, drop = FALSE] %*% object$coefficients[common])
  if (type == "link") return(eta)
  linkinv_r(eta, object$family$link_id)
}
