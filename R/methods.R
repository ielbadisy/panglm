#' @export
coef.panglm <- function(object, ...) object$coefficients

#' @export
vcov.panglm <- function(object, ...) object$vcov

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

  link_id <- object$family$link_id
  switch(link_id + 1L,
    eta,                              # 1: identity
    exp(eta),                         # 2: log
    1 / (1 + exp(-eta)),              # 3: logit
    stats::pnorm(eta)                 # 4: probit
  )
}
