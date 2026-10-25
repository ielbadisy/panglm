#' Family specifications for panglm
#'
#' Lightweight family objects used by [panglm()], analogous in spirit to
#' [stats::family] objects but keyed to the integer codes the C++ backend
#' expects (see `src/families.h`).
#'
#' @param link canonical link name
#' @return a list with `family`, `link`, `family_id`, `link_id`
#' @seealso [panglm()], whose `family` argument takes these objects, and
#'   [panglm-deprecated] for the former names.
#' @examples
#' gaussianfamily()
#' poissonfamily()
#' binomialfamily("probit")
#' negbinfamily()
#' gammafamily()
#' orderedfamily("probit")
#' betafamily()
#' tobitfamily(left = 0)
#'
#' data(copd)
#' fit <- panglm(hospitalized ~ crp, data = copd, index = c("id", "visit"),
#'               model = "within", family = binomialfamily("logit"))
#' coef(fit)
#' @name panglm-families
NULL

#' @rdname panglm-families
#' @export
gaussianfamily <- function(link = "identity") {
  link <- match.arg(link, "identity")
  list(family = "gaussian", link = link, family_id = 0L, link_id = 0L)
}

#' @rdname panglm-families
#' @export
poissonfamily <- function(link = "log") {
  link <- match.arg(link, "log")
  list(family = "poisson", link = link, family_id = 1L, link_id = 1L)
}

#' @rdname panglm-families
#' @export
binomialfamily <- function(link = c("logit", "probit")) {
  link <- match.arg(link)
  link_id <- if (link == "logit") 2L else 3L
  list(family = "binomial", link = link, family_id = 2L, link_id = link_id)
}

#' @rdname panglm-families
#' @export
negbinfamily <- function(link = "log") {
  link <- match.arg(link, "log")
  list(family = "negbin", link = link, family_id = 1L, link_id = 1L)
}

#' @details `gammafamily()` is the gamma GLM with a log link, the usual
#'   choice for skewed positive outcomes such as health care costs: the
#'   variance is proportional to the squared mean and the coefficients act
#'   multiplicatively on the mean, with no retransformation needed.
#' @rdname panglm-families
#' @export
gammafamily <- function(link = "log") {
  link <- match.arg(link, "log")
  list(family = "gamma", link = link, family_id = 3L, link_id = 1L)
}

#' @details `orderedfamily()` is the ordered logit or ordered probit model
#'   for ordinal outcomes such as self-assessed health:
#'   `P(y <= j) = F(c_j - x'beta - u_i)`. The response may be a factor, an
#'   ordered factor or integer-valued; categories follow the factor levels
#'   or the sorted values. The intercept is absorbed in the cutpoints `c_j`,
#'   which are reported after the slopes. `model = "pooling"` and
#'   `model = "random"` (Gaussian random intercept by Gauss-Hermite
#'   quadrature) are available for both links, `model = "within"` for the
#'   logit link only (the "blow-up and cluster" fixed-effects ordered logit
#'   of Baetschmann, Staub and Winkelmann 2015). The fit has class
#'   `c("panordered", "panglm")`; see [predict.panordered()].
#' @references Baetschmann, G., Staub, K. E. and Winkelmann, R. (2015).
#'   Consistent estimation of the fixed effects ordered logit model.
#'   *Journal of the Royal Statistical Society A*, 178, 685-703.
#'   \doi{10.1111/rssa.12090}
#' @rdname panglm-families
#' @export
orderedfamily <- function(link = c("logit", "probit")) {
  link <- match.arg(link)
  list(family = "ordered", link = link, family_id = NA_integer_,
       link_id = if (link == "logit") 2L else 3L)
}

#' @details `betafamily()` is beta regression for continuous outcomes in
#'   the open interval (0, 1), such as utility index scores rescaled to
#'   (0, 1), proportions of days in hospital or adherence fractions:
#'   `y ~ Beta(mu * phi, (1 - mu) * phi)` with `g(mu) = x'beta`, so `E(y) = mu`
#'   and `Var(y) = mu (1 - mu) / (1 + phi)`. The precision `phi` is
#'   reported as `fit$precision`. `model = "pooling"` and `model = "random"`
#'   (Gaussian random intercept on the link scale, adaptive Gauss-Hermite
#'   quadrature) are available; there is no fixed-effects estimator.
#'   Outcomes equal to 0 or 1 must be transformed first, e.g. with
#'   `(y * (n - 1) + 0.5) / n` (Smithson and Verkuilen 2006).
#' @references Ferrari, S. and Cribari-Neto, F. (2004). Beta regression for
#'   modelling rates and proportions. *Journal of Applied Statistics*, 31,
#'   799-815. \doi{10.1080/0266476042000214501}
#'
#'   Smithson, M. and Verkuilen, J. (2006). A better lemon squeezer?
#'   Maximum-likelihood regression with beta-distributed dependent
#'   variables. *Psychological Methods*, 11, 54-71.
#'   \doi{10.1037/1082-989X.11.1.54}
#' @rdname panglm-families
#' @export
betafamily <- function(link = c("logit", "probit")) {
  link <- match.arg(link)
  list(family = "beta", link = link, family_id = NA_integer_,
       link_id = if (link == "logit") 2L else 3L)
}

#' @keywords internal
#' @noRd
binomial_response_to_numeric <- function(y) {
  if (is.factor(y)) y <- as.numeric(y) - 1
  else if (is.character(y)) y <- as.numeric(factor(y)) - 1
  else if (is.logical(y)) y <- as.numeric(y)
  else y <- as.numeric(y)
  if (length(unique(y)) != 2) stop("the response must have exactly 2 distinct values for family = 'binomial'", call. = FALSE)
  if (min(y) != 0) y <- y - min(y)
  y / max(y)
}

resolve_family <- function(family) {
  if (is.character(family)) {
    family <- switch(family,
      gaussian = gaussianfamily(),
      poisson  = poissonfamily(),
      binomial = binomialfamily("logit"),
      probit   = binomialfamily("probit"),
      negbin   = negbinfamily(),
      gamma    = gammafamily(),
      Gamma    = gammafamily(),
      ologit   = orderedfamily("logit"),
      oprobit  = orderedfamily("probit"),
      beta     = betafamily(),
      tobit    = tobitfamily(),
      stop("unknown family '", family, "'", call. = FALSE)
    )
  }
  if (!is.list(family) || is.null(family$family)) {
    stop("family must be a character string or a panglm family spec ",
         "(see ?`panglm-families`)", call. = FALSE)
  }
  family
}
