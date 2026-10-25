#' Panel-aware lag(), lead() and diff() in model formulas
#'
#' Inside the formula of [panglm()], [paniv()], [panqr()] and [panfmm()],
#' `lag(x, k)`, `lead(x, k)` and `diff(x, lag, differences)` are taken
#' within each individual, on the grid of the sorted distinct values of
#' the time index (the convention of [pangmm()]), so the first `k` waves of
#' each individual, and waves whose predecessor is missing, get `NA` and
#' leave the estimation sample. `lag(x, 0:2)` gives `x` and its first two
#' lags as a matrix. Any expression of the data can be shifted, e.g.
#' `lag(log(income), 1)`.
#'
#' @keywords internal
#' @noRd
panel_call_names <- c("lag", "lead", "diff")

has_panel_calls <- function(...) {
  walk <- function(e) {
    if (!is.call(e)) return(FALSE)
    f <- e[[1L]]
    if (is.name(f) && as.character(f) %in% panel_call_names) return(TRUE)
    # index the call: as.list() of a terms object returns the object itself
    for (i in seq_along(e)[-1L]) {
      a <- e[[i]]
      if (!missing(a) && walk(a)) return(TRUE)
    }
    FALSE
  }
  any(vapply(list(...), function(f) !is.null(f) && walk(f), logical(1)))
}

# stop for estimators that do not shift within the panel
stop_panel_calls <- function(..., what) {
  if (has_panel_calls(...)) {
    stop("lag(), lead() and diff() inside the formula are not supported by ", what,
         ": they would not be taken within each individual. Create the shifted ",
         "columns in the data first", call. = FALSE)
  }
  invisible(NULL)
}

# environment whose lag/lead/diff shift vectors aligned with the rows of data
panel_shift_env <- function(data, index, parent) {
  if (length(index) < 2L) {
    stop("lag(), lead() and diff() in a formula need index = c(<individual>, <time>)",
         call. = FALSE)
  }
  shift <- panel_shift_function(data[[index[1L]]], data[[index[2L]]])
  n <- nrow(data)
  env <- new.env(parent = parent)
  check <- function(x) {
    if (NROW(x) != n) {
      stop("lag(), lead() and diff() in a formula apply to variables of the data (",
           n, " rows), not to a vector of length ", NROW(x), call. = FALSE)
    }
  }
  env$lag <- function(x, k = 1L, ...) {
    check(x)
    if (!is.numeric(k) || !length(k) || anyNA(k) || any(k != round(k))) {
      stop("the lag order 'k' must be a vector of integers", call. = FALSE)
    }
    if (length(k) == 1L) return(shift(x, k))
    if (!is.numeric(x) || is.matrix(x)) {
      stop("several lag orders need a numeric vector", call. = FALSE)
    }
    out <- vapply(k, function(j) shift(x, j), numeric(n))
    colnames(out) <- as.character(k)
    out
  }
  env$lead <- function(x, k = 1L, ...) env$lag(x, -k)
  env$diff <- function(x, lag = 1L, differences = 1L, ...) {
    check(x)
    if (!is.numeric(x)) stop("diff() in a formula needs a numeric variable", call. = FALSE)
    for (i in seq_len(differences)) x <- x - shift(x, lag)
    x
  }
  env
}

panel_shift_function <- function(id, time) {
  if (anyNA(id) || anyNA(time)) {
    stop("panel index columns must not contain missing values", call. = FALSE)
  }
  grid <- sort(unique(time))
  pos <- match(time, grid)
  m <- length(grid)
  key <- (match(id, unique(id)) - 1) * m + pos
  if (anyDuplicated(key)) {
    stop("lag(), lead() and diff() in a formula need one row per individual and ",
         "time period", call. = FALSE)
  }
  function(x, k) {
    k <- as.integer(k)
    target <- match(key - k, key)
    target[pos - k < 1L | pos - k > m] <- NA_integer_
    if (is.matrix(x)) x[target, , drop = FALSE] else x[target]
  }
}

#' Model frame with panel-aware lag(), lead() and diff()
#'
#' Without such calls this is `model.frame()`. With them the frame is
#' evaluated on all rows of `data` (so a lag can reach a row that is later
#' dropped), restricted to `rows`, then `na.action` is applied. The terms
#' keep the formula's own environment.
#'
#' @keywords internal
#' @noRd
panel_model_frame <- function(formula, data, index, na.action = stats::na.omit,
                              rows = NULL, xlev = NULL) {
  if (!has_panel_calls(formula)) {
    if (!is.null(rows)) data <- data[rows, , drop = FALSE]
    return(stats::model.frame(formula, data = data, na.action = na.action, xlev = xlev))
  }
  env0 <- environment(formula)
  environment(formula) <- panel_shift_env(data, index, env0 %||% globalenv())
  mf <- stats::model.frame(formula, data = data, na.action = stats::na.pass, xlev = xlev)
  terms_object <- attr(mf, "terms")
  attr(terms_object, ".Environment") <- env0
  if (!is.null(rows)) mf <- mf[rows, , drop = FALSE]
  if (is.null(na.action)) na.action <- stats::na.omit
  if (is.character(na.action)) na.action <- get(na.action, mode = "function")
  mf <- na.action(mf)
  na_attr <- attr(mf, "na.action")
  attr(mf, "terms") <- terms_object
  attr(mf, "na.action") <- na_attr
  mf
}

#' Design rows of new data for a fit with panel calls
#'
#' The rows of `newdata` are evaluated together with the rows of the
#' estimation data whose (individual, time) pair is not in `newdata`, so a
#' lag of the first new row reaches back into the estimation data.
#'
#' @keywords internal
#' @noRd
panel_newdata_frame <- function(terms_object, newdata, history, index, xlev) {
  if (is.null(history)) {
    stop("this fit has lag(), lead() or diff() in its formula but kept no estimation ",
         "data to rebuild them on new data", call. = FALSE)
  }
  missing_index <- setdiff(index, names(newdata))
  if (length(missing_index)) {
    stop("newdata must contain the panel index column(s) ",
         paste(missing_index, collapse = ", "),
         " to rebuild lag(), lead() and diff() terms", call. = FALSE)
  }
  vars <- intersect(names(history), names(newdata))
  nd <- as.data.frame(newdata)[, vars, drop = FALSE]
  key_new <- paste(nd[[index[1L]]], nd[[index[2L]]], sep = "\r")
  key_old <- paste(history[[index[1L]]], history[[index[2L]]], sep = "\r")
  extra <- history[!key_old %in% key_new, vars, drop = FALSE]
  combined <- rbind(nd, extra)
  rownames(combined) <- NULL
  panel_model_frame(terms_object, combined, index, na.action = stats::na.pass,
                    rows = seq_len(nrow(nd)), xlev = xlev)
}

#' Guard environment for formulas with panel calls
#'
#' A fit keeps its formula in this environment, so that code that rebuilds
#' the design with plain `model.frame()` stops instead of applying
#' `stats::lag()` (which does not shift) or `base::diff()` (which drops a
#' row) across individuals.
#'
#' @keywords internal
#' @noRd
panel_guard_env <- function(parent) {
  env <- new.env(parent = parent %||% globalenv())
  refuse <- function(...) {
    stop("this computation cannot rebuild the lag(), lead() or diff() terms of the ",
         "formula within each individual; create the shifted columns in the data ",
         "and refit", call. = FALSE)
  }
  env$lag <- env$lead <- env$diff <- refuse
  env
}

# give a formula with panel calls the guard environment
panel_guard_formula <- function(formula) {
  if (has_panel_calls(formula)) environment(formula) <- panel_guard_env(environment(formula))
  formula
}
