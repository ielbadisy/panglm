#' Family specifications for panglm
#'
#' Lightweight family objects used by [panglm()], analogous in spirit to
#' [stats::family] objects but keyed to the integer codes the C++ backend
#' expects (see `src/families.h`).
#'
#' @param link canonical link name
#' @return a list with `family`, `link`, `family_id`, `link_id`
#' @name panglm-families
NULL

#' @rdname panglm-families
#' @export
gaussian_family <- function(link = "identity") {
  link <- match.arg(link, "identity")
  list(family = "gaussian", link = link, family_id = 0L, link_id = 0L)
}

#' @rdname panglm-families
#' @export
poisson_family <- function(link = "log") {
  link <- match.arg(link, "log")
  list(family = "poisson", link = link, family_id = 1L, link_id = 1L)
}

#' @rdname panglm-families
#' @export
binomial_family <- function(link = c("logit", "probit")) {
  link <- match.arg(link)
  link_id <- if (link == "logit") 2L else 3L
  list(family = "binomial", link = link, family_id = 2L, link_id = link_id)
}

#' @rdname panglm-families
#' @export
negbin_family <- function(link = "log") {
  link <- match.arg(link, "log")
  list(family = "negbin", link = link, family_id = 1L, link_id = 1L)
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
      gaussian = gaussian_family(),
      poisson  = poisson_family(),
      binomial = binomial_family("logit"),
      probit   = binomial_family("probit"),
      negbin   = negbin_family(),
      stop("unknown family '", family, "'", call. = FALSE)
    )
  }
  if (!is.list(family) || is.null(family$family_id)) {
    stop("family must be a character string or a panglm family spec ",
         "(see ?`panglm-families`)", call. = FALSE)
  }
  family
}
