# tidy() and glance() methods for the panglm model classes other than
# "panglm" (see R/broom.R for that one), so every fit works with broom and
# modelsummary table-building workflows. Inference is Wald (normal) based,
# as in the print and summary methods of each class.

# Coefficient table from estimates and standard errors. The statistic and
# p-value are NA when the standard error is zero or missing (a reference
# category or a term without a variance estimate).
tidy_wald <- function(term, estimate, se, conf.int, conf.level, extra = NULL) {
  tidy_check_level(conf.level, "conf.level")
  estimate <- as.numeric(estimate)
  se <- as.numeric(se)
  ok <- is.finite(se) & se > 0
  zval <- ifelse(ok, estimate / se, NA_real_)
  out <- data.frame(term = as.character(term), estimate = estimate, std.error = se,
                    statistic = zval, p.value = 2 * stats::pnorm(-abs(zval)),
                    row.names = NULL, stringsAsFactors = FALSE)
  if (conf.int) {
    crit <- stats::qnorm(1 - (1 - conf.level) / 2)
    out$conf.low <- estimate - crit * se
    out$conf.high <- estimate + crit * se
  }
  if (!is.null(extra)) out <- cbind(out, extra, stringsAsFactors = FALSE)
  out
}

tidy_check_level <- function(level, name = "level") {
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1) {
    stop("'", name, "' must be a single number strictly between 0 and 1", call. = FALSE)
  }
  invisible(level)
}

tidy_vcov_se <- function(V, nms) {
  if (is.null(V)) return(rep(NA_real_, length(nms)))
  se <- sqrt(pmax(diag(V), 0))
  se[is.na(diag(V))] <- NA_real_
  as.numeric(se)
}

glance_ll <- function(object) {
  ll <- tryCatch(stats::logLik(object), error = function(e) NULL)
  if (is.null(ll) || !is.finite(as.numeric(ll))) {
    return(data.frame(logLik = NA_real_, AIC = NA_real_, BIC = NA_real_, df.logLik = NA_real_))
  }
  df <- attr(ll, "df")
  n <- attr(ll, "nobs") %||% stats::nobs(object)
  data.frame(logLik = as.numeric(ll), AIC = -2 * as.numeric(ll) + 2 * df,
             BIC = -2 * as.numeric(ll) + log(n) * df, df.logLik = df)
}

or_na <- function(x, type = NA_real_) if (is.null(x) || !length(x)) type else x

#' Tidy and glance methods for the panglm model classes
#'
#' Implement the `generics::tidy()` and `generics::glance()` generics for the
#' fits of [twopart()], [pangmm()], [paniv()], [panqr()], [panselect()],
#' [panfmm()], [attgt()], [etwfe()] and [eventstudy()], the hurdle model, and for the test
#' objects of the Hausman, joint Wald and dispersion tests, so that they work
#' with `broom` and `modelsummary`.
#'
#' `tidy()` returns one row per coefficient with columns `term`, `estimate`,
#' `std.error`, `statistic` (the Wald z), `p.value` and, with
#' `conf.int = TRUE`, the Wald interval `conf.low`, `conf.high`. The
#' estimates and standard errors are those of `coef()` and `vcov()`, so the
#' interval equals `confint()`. Extra columns identify the rows where a fit
#' has several blocks of coefficients: `component` (twopart: `"binary"` or
#' `"positive"`; panfmm: `"class"` or `"prior"`), `class` (panfmm), `tau`
#' (panqr: one row per term and quantile), `group`, `time`, `event`
#' (attgt, etwfe, eventstudy).
#'
#' For [attgt()], `type = "attgt"` (the default) gives one row per group-time
#' effect `ATT(g, t)`, named as in `coef()`; the other types are the
#' aggregations of [summary.panglm_attgt()]. For [etwfe()], `type =
#' "simple"` (the default) gives the overall ATT, `"group"`, `"calendar"` and
#' `"event"` the aggregations of [summary.panglm_etwfe()], and
#' `"coefficients"` the coefficients of the underlying regression. For
#' [eventstudy()] the reference period is a row with estimate and standard
#' error zero and a missing statistic. A statistic is missing whenever its
#' standard error is zero or missing.
#'
#' For the test objects, `tidy()` and `glance()` both return a one-row data
#' frame in the style of `broom`'s `htest` methods: `statistic`, `p.value`,
#' `parameter` (the degrees of freedom; for an F test also `den.df`),
#' `distribution` and `method`.
#'
#' `glance()` returns one row with the estimator, the variance type, the
#' numbers of observations and individuals, the log-likelihood, AIC and BIC
#' where the model is fitted by maximum likelihood (twopart, panfmm), and
#' the specification tests reported by `summary()` (for example the Sargan,
#' Hansen and Arellano-Bond tests of [pangmm()], the pre-trend test of
#' [attgt()] and [eventstudy()]).
#'
#' @param x a fitted model or test object
#' @param conf.int if `TRUE`, add `conf.low`/`conf.high` columns
#' @param conf.level confidence level for the interval
#' @param type for [attgt()] and [etwfe()] fits, which effects to tidy (see
#'   Details)
#' @param ... further arguments passed to the `summary()` aggregation
#'   (attgt, etwfe), otherwise unused
#' @return a `data.frame`: one row per coefficient (or effect) for `tidy()` of
#'   a fit, one row for `glance()` and for the test objects (see Details).
#' @seealso [tidy.panglm()] and [glance.panglm()] for [panglm()] fits.
#' @section Validation:
#' Numerical accuracy is checked in the package tests and in the validation
#' suite; see the validation article on the package website.
#' @examples
#' data(copd)
#' if (requireNamespace("generics", quietly = TRUE)) {
#'   set.seed(1)
#'   d <- data.frame(id = rep(1:100, each = 5), t = rep(1:5, 100), x = rnorm(500))
#'   d$y <- d$x + rep(rnorm(100), each = 5) + rnorm(500)
#'   fit <- panqr(y ~ x, data = d, index = c("id", "t"), tau = c(0.25, 0.75))
#'   generics::tidy(fit, conf.int = TRUE)
#'   generics::glance(fit)
#' }
#' @name tidy-methods
NULL

# ---- twopart ----------------------------------------------------------------

#' @rdname tidy-methods
#' @exportS3Method generics::tidy
tidy.twopart <- function(x, conf.int = FALSE, conf.level = 0.95, ...) {
  b <- stats::coef(x)
  comp <- sub(":.*$", "", names(b))
  tidy_wald(names(b), b, tidy_vcov_se(stats::vcov(x), names(b)), conf.int, conf.level,
            extra = data.frame(component = comp, stringsAsFactors = FALSE))
}

#' @rdname tidy-methods
#' @exportS3Method generics::glance
glance.twopart <- function(x, ...) {
  cbind(data.frame(model = x$model, link = x$link, positive = x$positive_type,
                   correlated = isTRUE(x$correlated), stringsAsFactors = FALSE),
        glance_ll(x),
        data.frame(nobs = x$nobs, n.positive = x$n_positive,
                   n.groups = or_na(x$binary$n_groups, NA_integer_)))
}

# ---- pangmm -----------------------------------------------------------------

#' @rdname tidy-methods
#' @exportS3Method generics::tidy
tidy.pangmm <- function(x, conf.int = FALSE, conf.level = 0.95, ...) {
  b <- stats::coef(x)
  tidy_wald(names(b), b, tidy_vcov_se(stats::vcov(x), names(b)), conf.int, conf.level)
}

#' @rdname tidy-methods
#' @exportS3Method generics::glance
glance.pangmm <- function(x, ...) {
  ar <- x$artest
  arp <- function(k) {
    if (is.null(ar) || !any(ar$order == k)) NA_real_ else ar$p.value[ar$order == k][1L]
  }
  data.frame(
    model = x$model,
    transformation = x$transformation,
    effect = x$effect,
    vcov.type = x$vcov_type,
    nobs = x$nobs,
    n.groups = x$n_groups,
    n.periods = x$n_periods,
    n.instruments = x$n_instruments,
    sargan = or_na(x$sargan$statistic), sargan.p.value = or_na(x$sargan$p.value),
    hansen = or_na(x$hansen$statistic), hansen.p.value = or_na(x$hansen$p.value),
    ar1.p.value = arp(1L), ar2.p.value = arp(2L),
    wald = or_na(x$wald$statistic), wald.p.value = or_na(x$wald$p.value),
    stringsAsFactors = FALSE
  )
}

# ---- paniv ------------------------------------------------------------------

#' @rdname tidy-methods
#' @exportS3Method generics::tidy
tidy.paniv <- function(x, conf.int = FALSE, conf.level = 0.95, ...) {
  b <- stats::coef(x)
  tidy_wald(names(b), b, tidy_vcov_se(stats::vcov(x), names(b)), conf.int, conf.level)
}

#' @rdname tidy-methods
#' @exportS3Method generics::glance
glance.paniv <- function(x, ...) {
  data.frame(
    model = x$model,
    effect = x$effect,
    vcov.type = x$vcov_type,
    sigma = if (is.null(x$sigma2)) NA_real_ else sqrt(x$sigma2),
    df.residual = or_na(x$df.residual),
    nobs = x$nobs,
    n.groups = x$n_groups,
    overid = or_na(x$overid$statistic), overid.df = or_na(x$overid$df),
    overid.p.value = or_na(x$overid$p.value),
    stringsAsFactors = FALSE
  )
}

# ---- panqr ------------------------------------------------------------------

#' @rdname tidy-methods
#' @exportS3Method generics::tidy
tidy.panqr <- function(x, conf.int = FALSE, conf.level = 0.95, ...) {
  B <- as.matrix(x$coefficients)
  rows <- lapply(seq_along(x$tau), function(j) {
    V <- if (is.list(x$vcov)) x$vcov[[j]] else NULL
    tidy_wald(rownames(B), B[, j], tidy_vcov_se(V, rownames(B)), conf.int, conf.level,
              extra = data.frame(tau = rep(x$tau[j], nrow(B))))
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' @rdname tidy-methods
#' @exportS3Method generics::glance
glance.panqr <- function(x, ...) {
  data.frame(model = x$model, vcov.type = x$vcov_type, n.tau = length(x$tau),
             nobs = x$nobs, n.groups = x$n_groups, stringsAsFactors = FALSE)
}

# ---- panselect --------------------------------------------------------------

#' @rdname tidy-methods
#' @exportS3Method generics::tidy
tidy.panselect <- function(x, conf.int = FALSE, conf.level = 0.95, ...) {
  b <- stats::coef(x)
  tidy_wald(names(b), b, tidy_vcov_se(stats::vcov(x), names(b)), conf.int, conf.level)
}

#' @rdname tidy-methods
#' @exportS3Method generics::glance
glance.panselect <- function(x, ...) {
  st <- if (x$vcov_type == "boot") x$selection_test$boot else x$selection_test$cluster
  data.frame(
    model = x$model,
    mundlak = isTRUE(x$mundlak),
    lambda = x$lambda,
    vcov.type = x$vcov_type,
    df.residual = or_na(x$df.residual),
    nobs = x$nobs,
    nobs.selection = x$nobs_selection,
    n.groups = x$n_groups,
    n.groups.selection = x$n_groups_selection,
    selection.statistic = or_na(st$statistic), selection.df = or_na(st$df),
    selection.p.value = or_na(st$p.value),
    stringsAsFactors = FALSE
  )
}

# ---- panfmm -----------------------------------------------------------------

#' @rdname tidy-methods
#' @exportS3Method generics::tidy
tidy.panfmm <- function(x, conf.int = FALSE, conf.level = 0.95, ...) {
  b <- stats::coef(x)
  nm <- names(b)
  prior <- grepl("^prior[0-9]+:", nm)
  cls <- as.integer(sub("^(class|prior)([0-9]+):.*$", "\\2", nm))
  tidy_wald(nm, b, tidy_vcov_se(stats::vcov(x), nm), conf.int, conf.level,
            extra = data.frame(component = ifelse(prior, "prior", "class"), class = cls,
                               stringsAsFactors = FALSE))
}

#' @rdname tidy-methods
#' @exportS3Method generics::glance
glance.panfmm <- function(x, ...) {
  cbind(data.frame(family = x$family$family, link = x$family$link, k = x$k,
                   vcov.type = x$vcov_type, stringsAsFactors = FALSE),
        glance_ll(x),
        data.frame(nobs = x$nobs, n.groups = x$n_groups,
                   converged = isTRUE(x$convergence == 0) || isTRUE(x$convergence)))
}

# ---- attgt ------------------------------------------------------------------

#' @rdname tidy-methods
#' @exportS3Method generics::tidy
tidy.panglm_attgt <- function(x, conf.int = FALSE, conf.level = 0.95,
                              type = c("attgt", "simple", "dynamic", "group", "calendar"),
                              ...) {
  type <- match.arg(type)
  if (type == "attgt") {
    est <- x$estimates
    return(tidy_wald(names(stats::coef(x)), est$estimate, est$se, conf.int, conf.level,
                     extra = data.frame(group = est$group, time = est$time,
                                        event = est$time - est$group)))
  }
  s <- summary(x, type = type, ...)
  if (type == "simple") {
    return(tidy_wald("ATT", s$overall$estimate, s$overall$se, conf.int, conf.level))
  }
  e <- s$estimates
  key <- c(dynamic = "event", group = "group", calendar = "time")[[type]]
  extra <- stats::setNames(data.frame(e$level), key)
  tidy_wald(sprintf("ATT(%s=%s)", substr(key, 1L, 1L), e$level), e$estimate, e$se,
            conf.int, conf.level, extra = extra)
}

#' @rdname tidy-methods
#' @exportS3Method generics::glance
glance.panglm_attgt <- function(x, ...) {
  pt <- x$pretrend
  data.frame(
    control = x$control,
    method = x$method,
    vcov.type = if (!is.null(x$boot)) "boot" else if (isTRUE(x$clustered)) "cluster" else "analytic",
    nobs = x$n,
    n.cohorts = length(unique(x$estimates$group)),
    n.periods = length(x$periods),
    pretrend.statistic = or_na(pt$statistic), pretrend.df = or_na(pt$df),
    pretrend.p.value = or_na(pt$p.value),
    stringsAsFactors = FALSE
  )
}

# ---- etwfe ------------------------------------------------------------------

#' @rdname tidy-methods
#' @exportS3Method generics::tidy
tidy.panglm_etwfe <- function(x, conf.int = FALSE, conf.level = 0.95,
                              type = c("simple", "group", "calendar", "event", "coefficients"),
                              ...) {
  type <- match.arg(type)
  if (type == "coefficients") {
    b <- stats::coef(x)
    return(tidy_wald(names(b), b, tidy_vcov_se(stats::vcov(x), names(b)), conf.int,
                     conf.level))
  }
  s <- summary(x, type = type, ...)
  if (is.null(s$estimates)) {
    return(tidy_wald("ATT", s$overall$estimate, s$overall$se, conf.int, conf.level))
  }
  e <- s$estimates
  key <- c(simple = "level", event = "event", group = "group", calendar = "time")[[type]]
  extra <- if (type == "simple") NULL else stats::setNames(data.frame(e$level), key)
  term <- if (type == "simple") rep("ATT", nrow(e)) else
    sprintf("ATT(%s=%s)", substr(key, 1L, 1L), e$level)
  if (!is.null(s$xvar)) {
    xv <- stats::setNames(data.frame(e[[s$xvar]]), s$xvar)
    extra <- if (is.null(extra)) xv else cbind(extra, xv)
    term <- sprintf("%s[%s=%s]", term, s$xvar, e[[s$xvar]])
  }
  tidy_wald(term, e$estimate, e$se, conf.int, conf.level, extra = extra)
}

#' @rdname tidy-methods
#' @exportS3Method generics::glance
glance.panglm_etwfe <- function(x, ...) {
  data.frame(family = x$family$family, control = x$control, vcov.type = x$vcov_type,
             ivar = isTRUE(x$ivar), nobs = x$nobs, n.groups = x$n_groups,
             n.cells = nrow(x$cells), stringsAsFactors = FALSE)
}

# ---- eventstudy -------------------------------------------------------------

#' @rdname tidy-methods
#' @exportS3Method generics::tidy
tidy.panglm_eventstudy <- function(x, conf.int = FALSE, conf.level = 0.95, ...) {
  e <- x$estimates
  tidy_wald(as.character(e$event), e$estimate, e$se, conf.int, conf.level,
            extra = data.frame(event = e$event))
}

#' @rdname tidy-methods
#' @exportS3Method generics::glance
glance.panglm_eventstudy <- function(x, ...) {
  pt <- x$pretrend
  data.frame(
    ref = x$ref, leads = x$leads, lags = x$lags, bin = isTRUE(x$bin),
    vcov.type = x$fit$vcov_type,
    nobs = x$fit$nobs, n.groups = x$fit$n_groups,
    pretrend.statistic = or_na(pt$statistic), pretrend.df = or_na(pt$df),
    pretrend.p.value = or_na(pt$p.value),
    stringsAsFactors = FALSE
  )
}

# ---- tests ------------------------------------------------------------------

tidy_test <- function(statistic, parameter, p.value, distribution, method) {
  out <- data.frame(statistic = as.numeric(statistic), p.value = as.numeric(p.value),
                    parameter = as.numeric(parameter[1L]), stringsAsFactors = FALSE)
  if (length(parameter) > 1L) out$den.df <- as.numeric(parameter[2L])
  out$distribution <- distribution
  out$method <- method
  out
}

#' @rdname tidy-methods
#' @exportS3Method generics::tidy
tidy.panglm_hausman <- function(x, ...) {
  tidy_test(x$statistic, x$parameter, x$p.value, x$distribution,
            if (x$method == "mundlak") "Hausman test (Mundlak)" else "Hausman test (classical)")
}

#' @rdname tidy-methods
#' @exportS3Method generics::glance
glance.panglm_hausman <- function(x, ...) tidy.panglm_hausman(x, ...)

#' @rdname tidy-methods
#' @exportS3Method generics::tidy
tidy.panglm_jointtest <- function(x, ...) {
  tidy_test(x$statistic, x$df, x$p.value, "chisq", "Joint Wald test")
}

#' @rdname tidy-methods
#' @exportS3Method generics::glance
glance.panglm_jointtest <- function(x, ...) tidy.panglm_jointtest(x, ...)

#' @rdname tidy-methods
#' @exportS3Method generics::tidy
tidy.panglm_dispersiontest <- function(x, ...) {
  out <- tidy_test(x$statistic, x$parameter, x$p.value, "chisq", "Pearson dispersion test")
  out$ratio <- x$ratio
  out
}

#' @rdname tidy-methods
#' @exportS3Method generics::glance
glance.panglm_dispersiontest <- function(x, ...) tidy.panglm_dispersiontest(x, ...)
