#' First treated period: checks and per-row relative time
#'
#' @keywords internal
#' @noRd
did_timing <- function(data, index, treat) {
  if (length(index) != 2L) stop("difference-in-differences needs index = c(individual, time)",
                                call. = FALSE)
  for (v in c(index, treat)) {
    if (!v %in% names(data)) stop("column '", v, "' not found in data", call. = FALSE)
  }
  time <- data[[index[2]]]
  if (!is.numeric(time)) stop("the time index must be numeric to define event time", call. = FALSE)
  first <- data[[treat]]
  if (!is.numeric(first)) {
    stop("'", treat, "' must hold the first treated period (numeric; NA, 0 or Inf for ",
         "never-treated units)", call. = FALSE)
  }
  first[!is.finite(first) | first == 0] <- NA
  # unused levels of a factor index (e.g. after subset()) are empty groups
  id <- data[[index[1]]]
  if (is.factor(id)) id <- droplevels(id)
  same <- tapply(first, id, function(v) length(unique(v)) == 1L)
  if (!all(same)) stop("'", treat, "' must be constant within individuals", call. = FALSE)
  list(time = time, first = first, rel = time - first, ever = !is.na(first))
}

#' Two-way fixed-effects difference-in-differences
#'
#' Adds the indicator `treated = 1(t >= first treated period)` to the model
#' and fits it with individual and time fixed effects (two-way fixed
#' effects, TWFE). With a single treatment date (the canonical two-group,
#' two-period design and its multi-period extension) the coefficient on
#' `treated` is the average treatment effect on the treated under parallel
#' trends. With staggered adoption and effects that vary over time or across
#' cohorts the TWFE coefficient is a weighted average that can put negative
#' weight on some effects (Goodman-Bacon 2021; de Chaisemartin and
#' D'Haultfoeuille 2020): use [eventstudy()] to look at the dynamics and a
#' heterogeneity-robust estimator ([attgt()] for Callaway and Sant'Anna
#' 2021, or [etwfe()] for the extended TWFE of Wooldridge 2021, 2023) when
#' adoption is staggered.
#'
#' The family can be Gaussian (effect on the mean), Poisson or negative
#' binomial (effect on the log rate: `exp(coef)` is a ratio, and Poisson FE
#' DiD needs only the conditional mean to be correct, Wooldridge 2023), or
#' binomial for a binary outcome. The binomial model is the two-way
#' conditional (Chamberlain) logit of [panglm()]: the individual effects
#' are conditioned out exactly and the time effects enter as dummies inside
#' the conditional likelihood, as in `survival::clogit(y ~ treated +
#' factor(time) + strata(id))`. Its coefficient `treated` is then a log
#' odds ratio, conditional on the individual effect (`exp(coef)` is the
#' odds ratio), parallel trends is assumed on the log-odds scale, and
#' individuals whose outcome never changes carry no information and are
#' dropped. It is not an effect on the probability scale: for average
#' effects on the treated in probability points use [etwfe()] with
#' `family = "binomial"`. Controls that do not vary within individuals are
#' absorbed by the individual effects and dropped.
#'
#' @param formula outcome and control variables, e.g. `y ~ age` or `y ~ 1`;
#'   the treatment indicator is added
#' @param data a data.frame in long format
#' @param index `c(individual, time)`; the time index must be numeric
#' @param treat name of the column holding each individual's first treated
#'   period (in the units of the time index); `NA`, `0` or `Inf` for units
#'   that are never treated
#' @param family `"gaussian"`, `"poisson"`, `"negbin"` or `"binomial"`
#'   (conditional logit, logit link only; see Details)
#' @param vcov covariance type, cluster-robust by individual by default
#' @param ... further arguments passed to [panglm()]
#' @return a [panglm()] fit whose coefficient `treated` is the DiD estimate (a
#'   log odds ratio for `family = "binomial"`), with a component `did` (number
#'   of treated and never-treated individuals, treatment cohorts). Its
#'   `call` is the `did()` call. All [panglm()] methods apply.
#' @references Goodman-Bacon, A. (2021). Difference-in-differences with
#'   variation in treatment timing. *Journal of Econometrics*, 225, 254-277.
#'   \doi{10.1016/j.jeconom.2021.03.014}
#'
#'   de Chaisemartin, C. and D'Haultfoeuille, X. (2020). Two-way fixed
#'   effects estimators with heterogeneous treatment effects. *American
#'   Economic Review*, 110, 2964-2996. \doi{10.1257/aer.20181169}
#'
#'   Wooldridge, J. M. (2023). Simple approaches to nonlinear
#'   difference-in-differences with panel data. *The Econometrics Journal*,
#'   26, C31-C66. \doi{10.1093/ectj/utad016}
#' @family policy evaluation
#' @seealso [summary.panglm()] and the other [panglm()] methods for the fitted
#'   model.
#' @section Validation:
#' Numerical accuracy is checked in the package tests and in the validation
#' suite; see the validation article on the package website.
#' @examples
#' set.seed(1)
#' d <- expand.grid(id = 1:200, year = 2015:2022)
#' d$first <- ifelse(d$id <= 100, 2019, NA)        # policy in 2019 for half
#' a <- rnorm(200)
#' d$visits <- 3 + a[d$id] + 0.1 * (d$year - 2015) +
#'   1.5 * (!is.na(d$first) & d$year >= d$first) + rnorm(nrow(d))
#' fit <- did(visits ~ 1, data = d, index = c("id", "year"), treat = "first")
#' summary(fit)
#' # binary outcome: conditional logit, the coefficient is a log odds ratio
#' d$admitted <- rbinom(nrow(d), 1, plogis(-1 + a[d$id] + 0.05 * (d$year - 2015) +
#'   0.8 * (!is.na(d$first) & d$year >= d$first)))
#' fitb <- did(admitted ~ 1, data = d, index = c("id", "year"), treat = "first",
#'             family = "binomial")
#' exp(coef(fitb))                                  # odds ratio
#' @export
did <- function(formula, data, index, treat, family = "gaussian",
                vcov = "cluster", ...) {
  timing <- did_timing(data, index, treat)
  if ("treated" %in% names(data)) stop("data already has a column 'treated'", call. = FALSE)
  data$treated <- as.numeric(timing$ever & timing$time >= timing$first)
  data$treated[is.na(data$treated)] <- 0
  f <- stats::update(formula, . ~ treated + .)
  fit <- panglm(f, data = data, index = index, model = "within", family = family,
                effect = "twoways", vcov = vcov, ...)
  # report the user's did() call, not the internal panglm() call
  fit$call <- match.call()
  ids <- data[[index[1]]]
  fit$did <- list(treat = treat,
                  treated_units = length(unique(ids[timing$ever])),
                  control_units = length(unique(ids[!timing$ever])),
                  cohorts = sort(unique(timing$first[timing$ever])))
  fit
}

#' Event-study regression with a pre-trend test
#'
#' Replaces the single treatment indicator of [did()] by indicators of event
#' time `e = t - first treated period`, from `-leads` to `lags`, with
#' `e = ref` (by default the period before treatment) as the omitted
#' reference, in a two-way fixed-effects model. Event times beyond the
#' window are binned into the end points (`bin = TRUE`) or dropped from the
#' sample (`bin = FALSE`). The lead coefficients estimate differences in
#' trends before treatment; their joint Wald test is the usual check of
#' parallel pre-trends (Freyaldenhoven et al. 2019). The lag coefficients
#' trace the effect over time since treatment.
#'
#' Identification of all event-time coefficients needs never-treated units
#' or a window narrower than the observed event times; without them a
#' second coefficient is dropped (with a warning). With staggered adoption
#' and heterogeneous effects the TWFE event-study coefficients are
#' contaminated by effects from other periods (Sun and Abraham 2021); use a
#' heterogeneity-robust estimator such as [attgt()] or [etwfe()] there.
#'
#' With `family = "binomial"` the model is the two-way conditional logit
#' (see [did()]) and the event-time coefficients are log odds ratios
#' relative to the reference event time.
#'
#' @inheritParams did
#' @param leads,lags number of pre- and post-treatment event times in the
#'   window (event times `-leads, ..., lags`)
#' @param ref omitted reference event time, default `-1`
#' @param bin logical; bin event times outside the window into the end
#'   points (`TRUE`) or drop those observations (`FALSE`)
#' @param level confidence level of the reported intervals
#' @return an object of class `"panglm_eventstudy"`: a list with `fit` (the
#'   [panglm()] fit), `estimates` (a data frame with `event`, `estimate`,
#'   `se`, `lower`, `upper`; the reference row is zero) and `pretrend` (the
#'   joint test of the leads, from [jointtest()]) and `call` (the
#'   `eventstudy()` call, also stored as the `call` of `fit`). Methods: `print`, `summary`,
#'   `coef`, `vcov`, `confint`, `logLik`, `nobs`, `plot`, `tidy` and `glance`.
#' @references Freyaldenhoven, S., Hansen, C. and Shapiro, J. M. (2019).
#'   Pre-event trends in the panel event-study design. *American Economic
#'   Review*, 109, 3307-3338. \doi{10.1257/aer.20180609}
#'
#'   Sun, L. and Abraham, S. (2021). Estimating dynamic treatment effects in
#'   event studies with heterogeneous treatment effects. *Journal of
#'   Econometrics*, 225, 175-199. \doi{10.1016/j.jeconom.2020.09.006}
#' @family policy evaluation
#' @seealso [eventstudy-methods] and [plot.panglm_eventstudy()] for the fitted
#'   model, [tidy-methods] for `tidy()` and `glance()`.
#' @section Validation:
#' Numerical accuracy is checked in the package tests and in the validation
#' suite; see the validation article on the package website.
#' @examples
#' set.seed(1)
#' d <- expand.grid(id = 1:200, year = 2015:2022)
#' d$first <- ifelse(d$id <= 100, 2019, NA)
#' a <- rnorm(200)
#' e <- d$year - d$first
#' d$visits <- 3 + a[d$id] + 0.1 * (d$year - 2015) +
#'   ifelse(!is.na(e) & e >= 0, 0.5 + 0.4 * e, 0) + rnorm(nrow(d))
#' es <- eventstudy(visits ~ 1, data = d, index = c("id", "year"),
#'                  treat = "first", leads = 4, lags = 3)
#' es
#' plot(es)
#' @export
eventstudy <- function(formula, data, index, treat, leads = 4, lags = 4, ref = -1,
                       bin = TRUE, family = "gaussian", vcov = "cluster", level = 0.95, ...) {
  leads <- as.integer(leads); lags <- as.integer(lags)
  if (leads < 1L || lags < 0L) stop("need leads >= 1 and lags >= 0", call. = FALSE)
  if (!ref %in% seq(-leads, lags)) stop("'ref' must lie in the event window", call. = FALSE)
  timing <- did_timing(data, index, treat)
  rel <- timing$rel
  if (bin) {
    rel <- pmin(pmax(rel, -leads), lags)
  } else {
    out <- timing$ever & (rel < -leads | rel > lags)
    data <- data[!out, , drop = FALSE]
    rel <- rel[!out]
  }
  events <- setdiff(seq(-leads, lags), ref)
  event_name <- function(e) ifelse(e < 0, paste0("event_lead", -e), paste0("event_lag", e))
  names_e <- event_name(events)
  clash <- intersect(names_e, names(data))
  if (length(clash)) stop("data already has column(s) ", paste(clash, collapse = ", "), call. = FALSE)
  for (j in seq_along(events)) {
    v <- as.numeric(rel == events[j]); v[is.na(v)] <- 0
    data[[names_e[j]]] <- v
  }
  rhs <- paste(c(names_e, "."), collapse = " + ")
  f <- stats::update(formula, stats::as.formula(paste(". ~", rhs)))
  fit <- panglm(f, data = data, index = index, model = "within", family = family,
                effect = "twoways", vcov = vcov, ...)
  fit$call <- match.call()
  b <- stats::coef(fit); V <- stats::vcov(fit)
  kept <- names_e %in% names(b)
  est <- rep(NA_real_, length(events)); se <- est
  est[kept] <- b[names_e[kept]]
  se[kept] <- sqrt(diag(V)[names_e[kept]])
  z <- stats::qnorm(1 - (1 - level) / 2)
  tab <- data.frame(event = c(events, ref), estimate = c(est, 0), se = c(se, 0))
  tab <- tab[order(tab$event), ]
  tab$lower <- tab$estimate - z * tab$se
  tab$upper <- tab$estimate + z * tab$se
  rownames(tab) <- NULL
  lead_terms <- names_e[events < 0 & kept]
  pretrend <- if (length(lead_terms)) jointtest(fit, terms = lead_terms) else NULL
  if (!is.null(pretrend)) pretrend$hypothesis <- "no pre-treatment differences in trends"
  structure(list(fit = fit, estimates = tab, pretrend = pretrend, ref = ref,
                 leads = leads, lags = lags, bin = bin, level = level, call = fit$call),
            class = "panglm_eventstudy")
}

#' @export
print.panglm_eventstudy <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nEvent study (two-way fixed effects, ", x$fit$family$family, ")\n", sep = "")
  cat("N =", x$fit$nobs, " groups =", x$fit$n_groups, " vcov:", x$fit$vcov_type, "\n")
  cat("Reference event time:", x$ref, if (x$bin) "(end points binned)" else "", "\n\n")
  print(format(x$estimates, digits = digits), row.names = FALSE)
  if (!is.null(x$pretrend)) {
    cat("\nPre-trend test: Chi-squared = ", format(x$pretrend$statistic, digits = digits),
        ", df = ", x$pretrend$df, ", p-value = ",
        format.pval(x$pretrend$p.value, digits = digits), "\n", sep = "")
  }
  invisible(x)
}

#' @export
coef.panglm_eventstudy <- function(object, ...) {
  stats::setNames(object$estimates$estimate, object$estimates$event)
}

#' Plot event-study coefficients
#'
#' @param x a `"panglm_eventstudy"` object from [eventstudy()]
#' @param xlab,ylab axis labels
#' @param ... passed to [graphics::plot()]
#' @return `x`, invisibly
#' @family eventstudy methods
#' @examples
#' set.seed(1)
#' d <- expand.grid(id = 1:100, year = 2015:2020)
#' d$first <- ifelse(d$id <= 50, 2018, NA)
#' e <- d$year - d$first
#' d$visits <- 3 + rnorm(100)[d$id] + ifelse(!is.na(e) & e >= 0, 1, 0) +
#'   rnorm(nrow(d))
#' es <- eventstudy(visits ~ 1, data = d, index = c("id", "year"),
#'                  treat = "first", leads = 3, lags = 2)
#' plot(es, main = "Effect of the policy by event time")
#' @export
plot.panglm_eventstudy <- function(x, xlab = "Event time", ylab = "Estimate", ...) {
  tab <- x$estimates
  ok <- is.finite(tab$estimate)
  graphics::plot(tab$event[ok], tab$estimate[ok], pch = 19, xlab = xlab, ylab = ylab,
                 ylim = range(tab$lower[ok], tab$upper[ok]), xaxt = "n", ...)
  graphics::axis(1, at = tab$event)
  graphics::segments(tab$event[ok], tab$lower[ok], tab$event[ok], tab$upper[ok])
  graphics::abline(h = 0, lty = 3)
  graphics::abline(v = x$ref + 0.5, lty = 2)
  invisible(x)
}
