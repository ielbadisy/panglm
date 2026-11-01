#' Plot a two-part model
#'
#' Coefficient estimates with Wald confidence intervals of the binary part
#' (probability of a positive outcome), of the positive part, or of both
#' side by side, in the dot-and-whisker layout of [plot.panglm()]. The
#' intervals use `vcov()`, so that they are the joint-likelihood intervals of
#' correlated two-part fits.
#'
#' @param x a fitted `"twopart"` object
#' @param which `"both"` (default), `"binary"` or `"positive"`
#' @param parm names or positions of the terms to plot; by default every
#'   term except `"(Intercept)"`. With `which = "both"`, names absent from
#'   one part are ignored in that part.
#' @param level confidence level of the intervals
#' @param ... additional arguments passed to [graphics::plot.default()]
#' @return `x`, invisibly
#' @family twopart methods
#' @examples
#' data(copd)
#' set.seed(1)
#' copd$cost <- ifelse(copd$exacerbations > 0,
#'                     exp(6 + 0.1 * copd$crp + rnorm(nrow(copd), 0, 0.5)), 0)
#' tp <- twopart(cost ~ treatment + age + crp, data = copd,
#'               index = c("id", "visit"), positive = "gamma")
#' plot(tp)
#' plot(tp, which = "positive", parm = "crp")
#' @export
plot.twopart <- function(x, which = c("both", "binary", "positive"), parm = NULL,
                         level = 0.95, ...) {
  which <- match.arg(which)
  est <- stats::coef(x)
  se <- .panglm_coef_se(est, stats::vcov(x))
  parts <- if (which == "both") c("binary", "positive") else which
  titles <- c(binary = paste0("Binary part (", x$link %||% "probit", ")"),
              positive = paste0("Positive part (", x$positive_type %||% "positive", ")"))
  panels <- lapply(parts, function(part) {
    prefix <- paste0(part, ":")
    sel <- startsWith(names(est), prefix)
    list(main = titles[[part]],
         parts = list(list(est = stats::setNames(est[sel], substring(names(est)[sel], nchar(prefix) + 1L)),
                           se = stats::setNames(se[sel], substring(names(se)[sel], nchar(prefix) + 1L)))))
  })
  .panglm_coef_panels(panels, parm, level, ...)
  invisible(x)
}

#' Plot a latent class (finite mixture) model
#'
#' `which = "coefficients"` compares the coefficients of the latent classes
#' in one dot-and-whisker panel (the classes dodged within each term, as in
#' [coefplot()]); the coefficients of the class-membership (concomitant)
#' model are not shown. `which = "posterior"` is a classification-quality
#' diagnostic: the histogram of each individual's largest posterior class
#' probability. Mass close to one means well-separated classes; the
#' subtitle reports its mean and the relative entropy
#' `1 - sum(-p log p) / (N log k)` (one for a perfect classification).
#'
#' @param x a fitted `"panfmm"` object
#' @param which `"coefficients"` (default) or `"posterior"`
#' @param parm names or positions of the terms to plot (without the class
#'   prefix); by default every term except `"(Intercept)"`
#' @param level confidence level of the intervals
#' @param main optional plot title
#' @param xlab optional horizontal-axis label
#' @param ... additional arguments passed to [graphics::plot.default()]
#' @return `x`, invisibly
#' @family panfmm methods
#' @examples
#' data(copd)
#' fit <- panfmm(exacerbations ~ crp + treatment, data = copd,
#'               index = c("id", "visit"), k = 2, family = "poisson",
#'               nstart = 2, seed = 1)
#' plot(fit)
#' plot(fit, which = "posterior")
#' @export
plot.panfmm <- function(x, which = c("coefficients", "posterior"), parm = NULL,
                        level = 0.95, main = NULL, xlab = NULL, ...) {
  which <- match.arg(which)
  if (which == "posterior") {
    post <- as.matrix(x$posterior)
    best <- apply(post, 1L, max)
    k <- ncol(post)
    ent <- -sum(post * log(pmax(post, 1e-300)))
    rel <- 1 - ent / (nrow(post) * log(k))
    graphics::hist(best, breaks = seq(1 / k, 1, length.out = 21L),
                   main = main %||% "Maximum posterior class probability",
                   xlab = xlab %||% "Largest posterior probability",
                   ylab = "Individuals", col = "grey80", border = "grey45", ...)
    graphics::mtext(sprintf("mean %.3f, relative entropy %.3f", mean(best), rel),
                    side = 3, line = 0.3, cex = 0.85 * graphics::par("cex"))
    return(invisible(x))
  }
  est <- stats::coef(x)
  se <- .panglm_coef_se(est, tryCatch(stats::vcov(x), error = function(e) NULL))
  cls <- regmatches(names(est), regexpr("^class[0-9]+:", names(est)))
  in_class <- grepl("^class[0-9]+:", names(est))
  classes <- unique(sub(":$", "", cls))
  parts <- lapply(classes, function(cl) {
    sel <- startsWith(names(est), paste0(cl, ":")) & in_class
    nm <- substring(names(est)[sel], nchar(cl) + 2L)
    list(est = stats::setNames(est[sel], nm), se = stats::setNames(se[sel], nm))
  })
  names(parts) <- sub("^class", "Class ", classes)
  tab <- .panglm_coef_matrices(parts, parm, level)
  .panglm_dotwhisker(tab, main = main %||% "Coefficients by latent class",
                     xlab = xlab %||% "Estimate", ...)
  invisible(x)
}

#' Plot the coefficients of a GMM, IV or selection model
#'
#' Coefficient estimates with Wald confidence intervals, in the
#' dot-and-whisker layout of [plot.panglm()]. For [panselect()] fits,
#' `which = "selection"` shows the coefficients of the per-period probits of
#' the selection equation (one marker per period, dodged within each term);
#' they are drawn without intervals because the first-step probits are
#' stored without standard errors. `which = "both"` shows the selection and
#' outcome equations side by side.
#'
#' @param x a fitted `"pangmm"`, `"paniv"` or `"panselect"` object
#' @param parm names or positions of the terms to plot; by default every
#'   term except `"(Intercept)"`
#' @param level confidence level of the intervals
#' @param main optional plot title
#' @param xlab optional horizontal-axis label
#' @param ... additional arguments passed to [graphics::plot.default()]
#' @return `x`, invisibly
#' @family pangmm methods
#' @examples
#' data(copd)
#' gm <- pangmm(fev1 ~ lag(fev1, 1) + crp | lag(fev1, 2:99), data = copd,
#'              index = c("id", "visit"))
#' plot(gm, parm = c("lag(fev1, 1)", "crp"))
#'
#' iv <- paniv(fev1 ~ crp + treatment | age + treatment, data = copd,
#'             index = c("id", "visit"), model = "pooling")
#' plot(iv)
#' @export
plot.pangmm <- function(x, parm = NULL, level = 0.95, main = NULL, xlab = NULL, ...) {
  .panglm_plot_single(x, parm, level, main, xlab, ...)
}

#' @rdname plot.pangmm
#' @family paniv methods
#' @export
plot.paniv <- function(x, parm = NULL, level = 0.95, main = NULL, xlab = NULL, ...) {
  .panglm_plot_single(x, parm, level, main, xlab, ...)
}

#' @rdname plot.pangmm
#' @param which for `panselect` fits, `"both"` (default), `"outcome"` or
#'   `"selection"`
#' @family panselect methods
#' @examples
#' set.seed(1)
#' d <- data.frame(id = rep(1:150, each = 4), wave = rep(1:4, 150))
#' d$treat <- rnorm(600)
#' d$reminder <- rnorm(600)
#' v <- rnorm(600)
#' d$respond <- as.numeric(0.3 + 0.5 * d$treat + d$reminder + v > 0)
#' d$qol <- ifelse(d$respond == 1, 1 + 0.5 * d$treat + 0.7 * v + rnorm(600), NA)
#' sel <- panselect(respond ~ treat + reminder, qol ~ treat, data = d,
#'                  index = c("id", "wave"), model = "pooling", boot = 0)
#' plot(sel)
#' @export
plot.panselect <- function(x, which = c("both", "outcome", "selection"), parm = NULL,
                           level = 0.95, xlab = NULL, ...) {
  which <- match.arg(which)
  panels <- list()
  if (which %in% c("both", "selection")) {
    probit <- Filter(Negate(is.null), x$probit)
    if (!length(probit)) stop("the fit stores no selection probits", call. = FALSE)
    sel_parts <- lapply(probit, function(b) {
      list(est = b, se = stats::setNames(rep(NA_real_, length(b)), names(b)))
    })
    names(sel_parts) <- paste("Period", names(probit))
    panels$selection <- list(main = "Selection probits", parts = sel_parts)
  }
  if (which %in% c("both", "outcome")) {
    est <- stats::coef(x)
    panels$outcome <- list(main = "Outcome equation",
                           parts = list(list(est = est, se = .panglm_coef_se(est, stats::vcov(x)))))
  }
  .panglm_coef_panels(panels, parm, level, xlab = xlab %||% "Estimate", ...)
  invisible(x)
}

.panglm_plot_single <- function(x, parm, level, main, xlab, ...) {
  est <- stats::coef(x)
  tab <- .panglm_coef_matrices(
    list(list(est = est, se = .panglm_coef_se(est, tryCatch(stats::vcov(x), error = function(e) NULL)))),
    parm, level
  )
  .panglm_dotwhisker(tab, main = main %||% "Coefficient estimates",
                     xlab = xlab %||% "Estimate", ...)
  invisible(x)
}
