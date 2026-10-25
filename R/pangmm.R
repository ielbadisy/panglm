#' Difference and system GMM for dynamic linear panel models
#'
#' Arellano and Bond (1991) difference GMM and Blundell and Bond (1998)
#' system GMM for linear panel models with a lagged dependent variable and
#' an unobserved individual effect, `y_it = sum_j rho_j y_i,t-j + x_it'b +
#' mu_i + e_it`. Both one-step and two-step estimators are available, with
#' robust one-step and Windmeijer (2005) corrected two-step standard
#' errors, the Arellano-Bond tests for first- and second-order serial
#' correlation of the differenced residuals and the Sargan and Hansen tests
#' of the overidentifying restrictions.
#'
#' @section Formula:
#' The formula follows [plm::pgmm()] and has two or three parts separated
#' by `|`:
#'
#' `y ~ lag(y, 1) + x + lag(x, 1) | lag(y, 2:99) | z`
#'
#' * The first part lists the regressors. `lag(v, k)` is the lag of order
#'   `k` of `v` within the individual, `k` may be a vector (`lag(v, 0:2)`
#'   gives `v` and its first two lags) and `lag(v)` means `lag(v, 1)`. Any
#'   expression of the data can be lagged, e.g. `lag(log(bmi), 1)`. The
#'   intercept is not estimated: the individual effects are differenced
#'   out, and time effects or the system-GMM constant are added through
#'   `effect`.
#' * The second part lists the GMM-style instruments: `lag(v, a:b)` uses,
#'   for the differenced equation of period `t`, the levels `v_(t-a)`, ...,
#'   `v_(t-b)` that exist in the panel (`b = 99` means all available
#'   lags). For an endogenous variable `a = 2`, for a predetermined one
#'   `a = 1`.
#' * The optional third part lists standard (IV-style) instruments, with
#'   the same `lag()` syntax; each enters in first differences in the
#'   differenced equation (and in levels in the level equation of system
#'   GMM). Without a third part, every regressor whose variable does not
#'   appear among the GMM-style instruments is used as its own standard
#'   instrument, as in [plm::pgmm()].
#'
#' Lags are taken on the time grid formed by the sorted distinct values of
#' the time index, so missing waves are respected and waves need not be
#' equally spaced in calendar time (`plm` instead lags by the numeric value
#' of the time index, so biennial waves coded as years have no lag one
#' period back there); rows with a missing
#' value in the response or a regressor are dropped from the moment
#' conditions and missing instruments are set to zero, as in `plm` and
#' Stata's `xtabond2`.
#'
#' @section Estimators:
#' With `transformation = "d"` (difference GMM) the model is estimated in
#' first differences, with the GMM-style instruments of each period in
#' their own columns (block-diagonal instrument matrix) or, with
#' `collapse = TRUE`, one column per lag (Roodman 2009). The one-step
#' weight matrix is `(sum_i Z_i' H Z_i)^-1` with `H` the tridiagonal
#' matrix with 2 on the diagonal and -1 off it, the covariance of the
#' differenced errors under homoskedasticity. With `transformation = "ld"`
#' (system GMM) the level equations of the same periods are stacked under
#' the differenced ones and instrumented by the first difference of each
#' GMM-style variable lagged `a - 1` periods (`a` its smallest lag); the
#' one-step weight uses the block matrix `H = [G, D; D', I]` of
#' [plm::pgmm()] (`fsm = "full"`). The two-step estimator uses the weight
#' `(sum_i Z_i' e_i e_i' Z_i)^-1` built from the one-step residuals.
#'
#' With `effect = "twoways"` (the default) time dummies are included: in
#' differences for difference GMM, and for system GMM a constant and time
#' dummies in the level equation (which instrument themselves there). With
#' `effect = "individual"` there are no time effects; system GMM then keeps
#' a constant in the level equation (`plm::pgmm()` omits it).
#'
#' @section Inference:
#' `vcov = "robust"` gives the one-step sandwich covariance or, for the
#' two-step estimator, the Windmeijer (2005) finite-sample corrected
#' covariance. `vcov = "classical"` gives, for the one-step estimator,
#' `s2 (X'Z S1^-1 Z'X)^-1` with `S1 = sum_i Z_i' H Z_i` and
#' `s2 = sum(de^2) / (2 (n_d - k))` estimated from the `n_d` differenced
#' one-step residuals `de`; for the two-step estimator it gives the
#' conventional (downward biased) `(X'Z S2^-1 Z'X)^-1` with
#' `S2 = sum_i Z_i' e1_i e1_i' Z_i`. (The one-step `vcov()` of
#' [plm::pgmm()] is `(X'Z S1^-1 Z'X)^-1 / N`, without the error variance.)
#'
#' The Arellano-Bond statistic of order `m` tests for autocorrelation of
#' order `m` in the differenced residuals; it is standard normal under the
#' null of no autocorrelation. First-order autocorrelation is expected by
#' construction; second-order autocorrelation indicates that `y_(t-2)` is
#' not a valid instrument. It is computed with the robust covariance of the
#' fit, as in `summary(plm::pgmm())`. The Sargan statistic
#' `(Z'e1)' (sum Z_i' H Z_i)^-1 (Z'e1) / s2` uses the one-step residuals
#' and is valid only under homoskedasticity; the Hansen J statistic
#' `(Z'e2)' (sum Z_i' e1_i e1_i' Z_i)^-1 (Z'e2)` uses the two-step
#' residuals and is robust. Both are chi-squared with `L - k` degrees of
#' freedom (`L` instruments, `k` coefficients); `plm::sargan()` returns
#' the Hansen statistic by default and `N s2` times the Sargan statistic
#' with `weights = "onestep"`. Many instruments weaken the
#' Hansen test (Roodman 2009); collapsing or limiting the lag depth keeps
#' the instrument count well below the number of individuals.
#'
#' @section Finite-sample bias:
#' The estimators are consistent as the number of individuals grows, but
#' in moderate samples they can be badly biased, and matching another
#' package does not protect against it. Two mechanisms matter. Too many
#' instruments overfit the endogenous regressors and pull the estimate
#' toward pooled OLS, while making the Hansen test look too good (Roodman
#' 2009). And when the variance of the individual effect is large relative
#' to the idiosyncratic variance, the level moments of system GMM become
#' weak and the estimate again drifts toward OLS (Bun and Windmeijer 2010).
#' In a simulated BMI panel (true persistence 0.6, effect SD 1.5, error
#' SD 1, 8 waves, 63 instruments) uncollapsed system GMM averaged 0.78
#' one-step and 0.72 two-step with 400 individuals, 0.65 and 0.61 with 2000
#' individuals; collapsing the instruments gave 0.62 already with 400.
#' Practical rules: prefer the two-step estimator, use `collapse = TRUE`
#' and a small `maxlag`, keep the instrument count well below the number
#' of individuals (`pangmm()` warns when it is not), treat a Hansen
#' p-value close to 1 as a warning sign rather than as support, and check
#' that the coefficient of the lagged outcome lies between the
#' fixed-effects (downward biased) and pooled OLS (upward biased)
#' estimates (Bond 2002).
#'
#' @param formula a formula `y ~ regressors | gmm instruments` or
#'   `y ~ regressors | gmm instruments | standard instruments`; see the
#'   Formula section
#' @param data a data.frame in long format
#' @param index length-2 character vector giving the individual and time
#'   columns
#' @param transformation `"d"` (difference GMM, Arellano-Bond) or `"ld"`
#'   (system GMM, Blundell-Bond)
#' @param model `"onestep"` or `"twostep"`
#' @param effect `"twoways"` (time dummies, the default) or `"individual"`
#' @param collapse logical: collapse the GMM-style instruments to one
#'   column per lag (and, for system GMM, one column for the level
#'   equation). A named logical vector applies it to selected GMM
#'   variables, named as in the formula, e.g. `c("log(bmi)" = TRUE)`.
#' @param maxlag optional positive integer: the largest lag of the
#'   GMM-style instruments, overriding larger upper lags in the formula
#'   (`lag(y, 2:99)` with `maxlag = 4` uses lags 2 to 4)
#' @param vcov `"robust"` (one-step sandwich or Windmeijer-corrected
#'   two-step, both clustered by individual) or `"classical"`; `"cluster"`
#'   is accepted as a synonym of `"robust"`, the name used by the other
#'   estimators of the package. See [vcov.pangmm()] to switch after
#'   fitting
#' @return an object of class `"pangmm"`, a list with components
#'   `coefficients` (the final estimates; time dummies are named by their
#'   period), `coefficients_onestep`, `vcov`, `vcov_type`, `vcov_list` (the
#'   classical and robust covariances), `residuals` (one numeric vector per
#'   individual, differenced rows first and, for system GMM, level rows
#'   after; rows not used are zero), `artest` (data frame of the
#'   Arellano-Bond tests), `sargan`, `hansen` and `wald` (lists with
#'   `statistic`, `df` and `p.value`; `wald` tests the joint significance
#'   of the coefficients other than time dummies and the constant),
#'   `sigma2` (`s2`), `n_instruments`, `nobs` (number of equations used),
#'   `nobs_diff`, `nobs_level`, `n_groups`, `n_periods`, `timedummies`
#'   (names of the time-effect coefficients), `transformation`, `model`,
#'   `effect`, `collapse`, `call` and `formula`. It has `print`,
#'   `summary`, `coef`, `vcov`, `nobs`, `residuals` and `confint` methods.
#' @references
#' Bond, S. (2002). Dynamic panel data models: a guide to micro data
#' methods and practice. *Portuguese Economic Journal*, 1, 141-162.
#'
#' Bun, M. J. G. and Windmeijer, F. (2010). The weak instrument problem of
#' the system GMM estimator in dynamic panel data models. *The
#' Econometrics Journal*, 13, 95-126.
#'
#' Arellano, M. and Bond, S. (1991). Some tests of specification for panel
#' data: Monte Carlo evidence and an application to employment equations.
#' *Review of Economic Studies*, 58, 277-297.
#'
#' Blundell, R. and Bond, S. (1998). Initial conditions and moment
#' restrictions in dynamic panel data models. *Journal of Econometrics*,
#' 87, 115-143.
#'
#' Windmeijer, F. (2005). A finite sample correction for the variance of
#' linear efficient two-step GMM estimators. *Journal of Econometrics*,
#' 126, 25-51.
#'
#' Roodman, D. (2009). How to do xtabond2: An introduction to difference
#' and system GMM in Stata. *Stata Journal*, 9, 86-136.
#' @examples
#' # Habit persistence in smoking: cigarettes per day depend on last
#' # wave's consumption (state dependence) and on the price of tobacco,
#' # with a permanent individual propensity to smoke. The lagged outcome is
#' # correlated with that propensity, so pooled OLS overstates persistence
#' # and the within estimator understates it (Nickell bias); GMM is
#' # consistent. True values: persistence 0.5, price -0.3.
#' set.seed(1)
#' N <- 300; TT <- 7
#' d <- data.frame(id = rep(seq_len(N), each = TT), wave = rep(seq_len(TT), N))
#' propensity <- rnorm(N)
#' d$price <- rnorm(N * TT)
#' d$cigs <- NA
#' for (i in seq_len(N)) {
#'   r <- (i - 1) * TT + seq_len(TT)
#'   y <- propensity[i] / 0.5 + rnorm(1)
#'   for (t in seq_len(TT)) {
#'     y <- 0.5 * y - 0.3 * d$price[r[t]] + propensity[i] + rnorm(1)
#'     d$cigs[r[t]] <- y
#'   }
#' }
#'
#' # Arellano-Bond difference GMM, two-step with Windmeijer standard errors
#' ab <- pangmm(cigs ~ lag(cigs, 1) + price | lag(cigs, 2:99), data = d,
#'              index = c("id", "wave"), model = "twostep")
#' summary(ab)
#'
#' # Blundell-Bond system GMM with collapsed instruments and at most three
#' # lags of the outcome as instruments
#' bb <- pangmm(cigs ~ lag(cigs, 1) + price | lag(cigs, 2:99), data = d,
#'              index = c("id", "wave"), transformation = "ld",
#'              model = "twostep", collapse = TRUE, maxlag = 4)
#' summary(bb)
#'
#' # Employment equation of Arellano and Bond (1991)
#' if (requireNamespace("plm", quietly = TRUE)) {
#'   data("EmplUK", package = "plm")
#'   emp <- pangmm(log(emp) ~ lag(log(emp), 1:2) + lag(log(wage), 0:1) +
#'                   log(capital) + lag(log(output), 0:1) | lag(log(emp), 2:99),
#'                 data = EmplUK, index = c("firm", "year"), model = "twostep")
#'   summary(emp)
#' }
#' @export
pangmm <- function(formula, data, index,
                   transformation = c("d", "ld"),
                   model = c("onestep", "twostep"),
                   effect = c("twoways", "individual"),
                   collapse = FALSE, maxlag = NULL,
                   vcov = c("robust", "classical", "cluster")) {
  transformation <- match.arg(transformation)
  model <- match.arg(model)
  effect <- match.arg(effect)
  vcov <- match.arg(vcov)
  if (vcov == "cluster") vcov <- "robust"
  if (missing(index)) stop("'index' is required, e.g. index = c(\"id\", \"time\")", call. = FALSE)
  if (!is.character(index) || length(index) != 2L || anyDuplicated(index)) {
    stop("'index' must contain two distinct column names, c(id, time)", call. = FALSE)
  }
  data <- as.data.frame(data)
  missing_index <- setdiff(index, names(data))
  if (length(missing_index)) {
    stop("index column(s) not found in data: ", paste(missing_index, collapse = ", "),
         call. = FALSE)
  }
  if (!is.null(maxlag) &&
      (!is.numeric(maxlag) || length(maxlag) != 1L || is.na(maxlag) || maxlag < 1)) {
    stop("'maxlag' must be a positive integer", call. = FALSE)
  }

  parts <- gmm_split_formula(formula)
  env <- environment(formula) %||% parent.frame()
  main <- gmm_terms(parts$main, merge = FALSE)
  if (!length(main)) stop("the formula has no regressors", call. = FALSE)
  gmmv <- gmm_terms(parts$gmm, merge = TRUE)
  if (!length(gmmv)) stop("the second part of the formula must list GMM-style instruments",
                          call. = FALSE)
  if (any(vapply(gmmv, function(v) min(v$lags), numeric(1)) < 0)) {
    stop("GMM-style instrument lags must be non-negative", call. = FALSE)
  }
  if (!is.null(maxlag)) {
    for (j in seq_along(gmmv)) {
      if (min(gmmv[[j]]$lags) > maxlag) {
        stop("'maxlag' is smaller than the first lag of ", names(gmmv)[j], call. = FALSE)
      }
      gmmv[[j]]$lags <- gmmv[[j]]$lags[gmmv[[j]]$lags <= maxlag]
    }
  }
  ivv <- if (!is.null(parts$iv)) {
    gmm_terms(parts$iv, merge = TRUE)
  } else {
    keep <- !names(main) %in% names(gmmv)
    if (any(keep)) gmm_merge(main[keep]) else list()
  }
  collapse <- gmm_collapse_arg(collapse, names(gmmv))

  # balanced N x T grid of every variable used
  id <- data[[index[1]]]
  tm <- data[[index[2]]]
  if (anyNA(id) || anyNA(tm)) stop("panel index columns must not contain missing values",
                                   call. = FALSE)
  uid <- sort(unique(id))
  utm <- sort(unique(tm))
  ii <- match(id, uid)
  tt <- match(tm, utm)
  if (anyDuplicated(cbind(ii, tt))) {
    stop("each individual-time index pair must identify a unique observation", call. = FALSE)
  }
  N <- length(uid)
  TT <- length(utm)
  lhs <- parts$lhs
  exprs <- c(list(lhs), lapply(main, `[[`, "expr"), lapply(gmmv, `[[`, "expr"),
             lapply(ivv, `[[`, "expr"))
  labels <- vapply(exprs, gmm_deparse, character(1))
  exprs <- exprs[!duplicated(labels)]
  labels <- labels[!duplicated(labels)]
  grids <- lapply(exprs, function(e) {
    v <- eval(e, data, env)
    if (is.logical(v)) v <- as.numeric(v)
    if (!is.numeric(v) || NCOL(v) != 1L || length(v) != nrow(data)) {
      stop("'", gmm_deparse(e), "' must evaluate to a numeric vector with one value per row",
           call. = FALSE)
    }
    M <- matrix(NA_real_, N, TT)
    M[cbind(ii, tt)] <- as.numeric(v)
    M[!is.finite(M)] <- NA_real_
    M
  })
  names(grids) <- labels
  lagged <- function(label, k) {
    M <- grids[[label]]
    if (k == 0) return(M)
    out <- matrix(NA_real_, N, TT)
    if (k < TT) out[, (k + 1):TT] <- M[, seq_len(TT - k), drop = FALSE]
    out
  }

  # first period with a differenced equation: TL1 + 1, as in plm::pgmm()
  main_maxlag <- max(unlist(lapply(main, `[[`, "lags")))
  gmm_minlag <- min(unlist(lapply(gmmv, `[[`, "lags")))
  TL1 <- max(main_maxlag + 1, gmm_minlag)
  m <- TT - TL1
  if (m < 1L) stop("not enough time periods for the lags in the formula", call. = FALSE)
  pd <- (TL1 + 1):TT
  pl <- if (transformation == "ld") TL1:TT else integer(0)
  R <- m + length(pl)
  sys <- transformation == "ld"
  zeros_l <- matrix(0, N, length(pl))

  # one N x R block per column: differenced rows first, then level rows
  dpart <- function(M) M[, pd, drop = FALSE] - M[, pd - 1L, drop = FALSE]
  lpart <- function(M) M[, pl, drop = FALSE]
  both <- function(M) cbind(dpart(M), lpart(M))
  stack <- function(blocks) {
    if (!length(blocks)) return(matrix(0, N * R, 0))
    out <- vapply(blocks, function(B) as.vector(t(B)), numeric(N * R))
    if (!is.matrix(out)) out <- matrix(out, nrow = N * R)
    colnames(out) <- names(blocks)
    out
  }

  ylab <- gmm_deparse(lhs)
  Yb <- both(grids[[ylab]])
  Xb <- list()
  for (term in main) {
    lab <- gmm_deparse(term$expr)
    for (k in term$lags) Xb[[gmm_lagname(lab, k)]] <- both(lagged(lab, k))
  }
  Xb <- Xb[!duplicated(names(Xb))]
  if (ylab %in% names(Xb)) stop("the response cannot be a contemporaneous regressor", call. = FALSE)

  # GMM-style instruments
  Wb <- list()
  for (g in names(gmmv)) {
    lags <- gmmv[[g]]$lags
    a <- min(lags)
    b <- max(lags)
    G <- grids[[g]]
    if (collapse[[g]]) {
      for (l in a:min(b, TT - 1)) {
        B <- matrix(0, N, m)
        for (r in seq_len(m)) if (pd[r] - l >= 1) B[, r] <- G[, pd[r] - l]
        Wb[[sprintf("%s.L%d", g, l)]] <- cbind(B, zeros_l)
      }
    } else {
      for (r in seq_len(m)) {
        p <- pd[r]
        if (min(b, p - 1) < a) next
        for (l in rev(seq(a, min(b, p - 1)))) {
          B <- matrix(0, N, m)
          B[, r] <- G[, p - l]
          Wb[[sprintf("%s.L%d.t%s", g, l, utm[p])]] <- cbind(B, zeros_l)
        }
      }
    }
    if (sys) {
      s <- max(a - 1, 0)
      Dg <- matrix(NA_real_, N, TT)
      Dg[, -1] <- G[, -1, drop = FALSE] - G[, -TT, drop = FALSE]
      avail <- pl[pl - s >= 2]
      if (collapse[[g]]) {
        B <- zeros_l
        for (p in avail) B[, match(p, pl)] <- Dg[, p - s]
        Wb[[sprintf("%s.D%d", g, s)]] <- cbind(matrix(0, N, m), B)
      } else {
        for (p in avail) {
          B <- zeros_l
          B[, match(p, pl)] <- Dg[, p - s]
          Wb[[sprintf("%s.D%d.t%s", g, s, utm[p])]] <- cbind(matrix(0, N, m), B)
        }
      }
    }
  }
  # standard instruments: differences, and levels in a separate block
  for (v in names(ivv)) {
    for (k in ivv[[v]]$lags) {
      M <- lagged(v, k)
      nm <- gmm_lagname(v, k)
      Wb[[paste0("diff.", nm)]] <- cbind(dpart(M), zeros_l)
      if (sys) Wb[[paste0("level.", nm)]] <- cbind(matrix(0, N, m), lpart(M))
    }
  }

  # time effects and the system-GMM constant
  timedummies <- character(0)
  Tb <- list()
  if (effect == "twoways") {
    if (sys) Tb[["(Intercept)"]] <- cbind(matrix(0, N, m), matrix(1, N, length(pl)))
    for (q in pd) {
      D <- matrix(0, N, TT)
      D[, q] <- 1
      Tb[[as.character(utm[q])]] <- both(D)
    }
  } else if (sys) {
    Tb[["(Intercept)"]] <- cbind(matrix(0, N, m), matrix(1, N, length(pl)))
  }
  timedummies <- names(Tb)

  y <- as.vector(t(Yb))
  X <- stack(c(Xb, Tb))
  W <- stack(Wb)
  Wt <- stack(Tb)
  if (sys && ncol(Wt)) {
    # time effects instrument themselves in the level equation only
    Wt[rep(seq_len(R), N) <= m, ] <- 0
  }
  W <- cbind(W, Wt)

  narow <- is.na(y) | rowSums(is.na(X)) > 0
  y[narow] <- 0
  X[narow, ] <- 0
  W[is.na(W)] <- 0
  W[narow, ] <- 0
  used <- !narow
  rowtype <- rep(c(rep("d", m), rep("l", length(pl))), N)
  group <- rep(seq_len(N), each = R)

  # drop instruments that are zero for every individual or redundant
  nz <- colSums(W != 0) > 0
  W <- W[, nz, drop = FALSE]
  W <- gmm_independent(W)
  if (any(colSums(X != 0) == 0)) {
    stop("regressor(s) with no usable observation: ",
         paste(colnames(X)[colSums(X != 0) == 0], collapse = ", "), call. = FALSE)
  }
  k <- ncol(X)
  L <- ncol(W)
  if (L < k) stop("insufficient number of instruments: ", L, " for ", k, " coefficients",
                  call. = FALSE)
  if (L >= N) {
    warning(L, " instruments for ", N, " individuals: too many instruments bias GMM ",
            "toward OLS and weaken the Hansen test; use collapse = TRUE or a smaller ",
            "maxlag", call. = FALSE)
  }

  H <- gmm_fsm(m, sys)
  fit <- gmm_estimate(y, X, W, H, group, R, model)
  if (!all(is.finite(fit$coefficients))) stop("GMM estimation failed", call. = FALSE)

  dres <- fit$e1[used & rowtype == "d"]
  n_d <- sum(used & rowtype == "d")
  sigma2 <- sum(dres^2) / (2 * (n_d - k))
  fit$vcov_list$classical <- if (model == "onestep") fit$M1 * sigma2 else fit$V2

  out <- list(
    coefficients = fit$coefficients,
    coefficients_onestep = fit$b1,
    vcov_list = fit$vcov_list,
    sigma2 = sigma2,
    n_instruments = L,
    nobs = sum(used),
    nobs_diff = n_d,
    nobs_level = sum(used & rowtype == "l"),
    n_groups = N,
    n_periods = TT,
    timedummies = timedummies,
    transformation = transformation,
    model = model,
    effect = effect,
    collapse = unlist(collapse),
    maxlag = maxlag,
    index = index,
    call = match.call(),
    formula = formula
  )
  e <- if (model == "onestep") fit$e1 else fit$e2
  out$residuals <- split(e, group)
  names(out$residuals) <- as.character(uid)

  # Sargan (one-step, homoskedastic) and Hansen J (two-step) statistics
  z1 <- crossprod(W, fit$e1)
  s_stat <- as.numeric(crossprod(z1, fit$S1i %*% z1)) / sigma2
  z2 <- crossprod(W, fit$e2)
  h_stat <- as.numeric(crossprod(z2, fit$S2i %*% z2))
  df <- L - k
  out$sargan <- list(statistic = s_stat, df = df,
                     p.value = stats::pchisq(s_stat, df, lower.tail = FALSE))
  out$hansen <- list(statistic = h_stat, df = df,
                     p.value = stats::pchisq(h_stat, df, lower.tail = FALSE))

  # Arellano-Bond serial-correlation tests on the differenced residuals
  vr <- fit$vcov_list$robust
  A <- if (model == "onestep") fit$S1i else fit$S2i
  B <- if (model == "onestep") fit$M1 else fit$V2
  ar <- lapply(1:2, function(ord) {
    if (ord >= m) return(NULL)
    gmm_artest(e, X, W, group, rowtype, m, R, N, ord, A, B, vr)
  })
  ar <- do.call(rbind, ar)
  out$artest <- ar
  class(out) <- "pangmm"
  out$vcov_type <- vcov
  out$vcov <- out$vcov_list[[vcov]]
  out$wald <- gmm_wald(out, out$vcov)
  out
}

#' Split `y ~ a | b` or `y ~ a | b | c` by hand (no Formula dependency)
#'
#' @keywords internal
#' @noRd
gmm_split_formula <- function(formula) {
  if (!inherits(formula, "formula") || length(formula) != 3L) {
    stop("'formula' must be a two-sided formula 'y ~ regressors | gmm instruments'",
         call. = FALSE)
  }
  bar <- function(e) is.call(e) && identical(e[[1L]], as.name("|"))
  rhs <- formula[[3L]]
  if (!bar(rhs)) {
    stop("'formula' must have at least two parts, 'y ~ regressors | gmm instruments'",
         call. = FALSE)
  }
  if (bar(rhs[[2L]])) {
    if (bar(rhs[[2L]][[2L]])) stop("'formula' has at most three parts", call. = FALSE)
    out <- list(main = rhs[[2L]][[2L]], gmm = rhs[[2L]][[3L]], iv = rhs[[3L]])
  } else {
    out <- list(main = rhs[[2L]], gmm = rhs[[3L]], iv = NULL)
  }
  if (bar(out$gmm) || (!is.null(out$iv) && bar(out$iv))) {
    stop("could not parse the parts of 'formula'", call. = FALSE)
  }
  out$lhs <- formula[[2L]]
  out
}

#' Parse the terms of one formula part into variable expressions and lags
#'
#' @keywords internal
#' @noRd
gmm_terms <- function(part, merge) {
  labs <- attr(stats::terms(stats::as.formula(call("~", part))), "term.labels")
  out <- lapply(labs, function(lab) {
    e <- str2lang(lab)
    if (is.call(e) && identical(e[[1L]], as.name("lag"))) {
      if (length(e) < 2L || length(e) > 3L) stop("lag() takes a variable and lags: ", lab,
                                                   call. = FALSE)
      lags <- if (length(e) == 3L) eval(e[[3L]], baseenv()) else 1
      if (!is.numeric(lags) || anyNA(lags) || any(lags != round(lags)) || any(lags < 0)) {
        stop("invalid lags in ", lab, call. = FALSE)
      }
      list(expr = e[[2L]], lags = as.integer(lags))
    } else {
      list(expr = e, lags = 0L)
    }
  })
  names(out) <- vapply(out, function(x) gmm_deparse(x$expr), character(1))
  if (merge) gmm_merge(out) else out
}

#' Merge the lags of repeated variables
#'
#' @keywords internal
#' @noRd
gmm_merge <- function(terms) {
  nm <- unique(names(terms))
  out <- lapply(nm, function(v) {
    same <- terms[names(terms) == v]
    list(expr = same[[1L]]$expr, lags = sort(unique(unlist(lapply(same, `[[`, "lags")))))
  })
  names(out) <- nm
  out
}

gmm_deparse <- function(e) paste(deparse(e, width.cutoff = 500L), collapse = "")

gmm_lagname <- function(lab, k) if (k == 0) lab else sprintf("lag(%s, %d)", lab, k)

#' Validate and expand the collapse argument to one flag per GMM variable
#'
#' @keywords internal
#' @noRd
gmm_collapse_arg <- function(collapse, vars) {
  if (!is.logical(collapse) || anyNA(collapse)) stop("'collapse' must be logical", call. = FALSE)
  if (is.null(names(collapse))) {
    if (length(collapse) == 1L) collapse <- rep(collapse, length(vars))
    if (length(collapse) != length(vars)) {
      stop("'collapse' must have length 1 or one value per GMM variable", call. = FALSE)
    }
    names(collapse) <- vars
  } else {
    bad <- setdiff(names(collapse), vars)
    if (length(bad)) stop("unknown names in 'collapse': ", paste(bad, collapse = ", "),
                          call. = FALSE)
    full <- stats::setNames(rep(FALSE, length(vars)), vars)
    full[names(collapse)] <- collapse
    collapse <- full
  }
  as.list(collapse)
}

#' Keep the linearly independent instrument columns (in their order)
#'
#' @keywords internal
#' @noRd
gmm_independent <- function(W, tol = 1e-9) {
  if (!ncol(W)) return(W)
  sc <- sqrt(colSums(W^2))
  dec <- qr(sweep(W, 2L, sc, "/"), tol = tol)
  if (dec$rank == ncol(W)) return(W)
  W[, sort(dec$pivot[seq_len(dec$rank)]), drop = FALSE]
}

#' One-step weighting structure H (plm fsm "G" and "full")
#'
#' @keywords internal
#' @noRd
gmm_fsm <- function(m, sys) {
  G <- diag(2, m)
  if (m > 1L) {
    G[cbind(1:(m - 1), 2:m)] <- -1
    G[cbind(2:m, 1:(m - 1))] <- -1
  }
  if (!sys) return(G)
  D <- matrix(0, m, m + 1L)
  D[cbind(1:m, 1:m)] <- -1
  D[cbind(1:m, 2:(m + 1L))] <- 1
  rbind(cbind(G, D), cbind(t(D), diag(1, m + 1L)))
}

#' Inverse of a GMM weighting matrix, generalized when it is singular
#'
#' Same rule as plm::pgmm(): a Moore-Penrose inverse (as MASS::ginv) when
#' the smallest eigenvalue is below 1e-9, which happens when there are more
#' instruments than individuals.
#'
#' @keywords internal
#' @noRd
gmm_inverse <- function(S, step) {
  ev <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) >= 1e-9) return(solve(S))
  warning("the ", step, " weighting matrix is singular; a generalized inverse is used",
          call. = FALSE)
  sv <- svd(S)
  pos <- sv$d > max(sqrt(.Machine$double.eps) * sv$d[1L], 0)
  sv$v[, pos, drop = FALSE] %*% (t(sv$u[, pos, drop = FALSE]) / sv$d[pos])
}

#' One-step and two-step GMM, robust and Windmeijer covariances
#'
#' Data are stacked individual by individual with R rows each.
#'
#' @keywords internal
#' @noRd
gmm_estimate <- function(y, X, W, H, group, R, model) {
  N <- length(y) / R
  L <- ncol(W)
  HW <- matrix(H %*% matrix(W, nrow = R), nrow = N * R)
  S1 <- crossprod(W, HW)
  S1 <- (S1 + t(S1)) / 2
  ZX <- crossprod(W, X)
  Zy <- crossprod(W, y)
  S1i <- gmm_inverse(S1, "first-step")
  S1iZX <- S1i %*% ZX
  M1 <- solve(crossprod(ZX, S1iZX))
  b1 <- as.numeric(M1 %*% crossprod(S1iZX, Zy))
  names(b1) <- colnames(X)
  e1 <- as.numeric(y - X %*% b1)
  Q1 <- rowsum(W * e1, group, reorder = FALSE)
  S2 <- crossprod(Q1)
  S2i <- gmm_inverse(S2, "second-step")
  S2iZX <- S2i %*% ZX
  V2 <- solve(crossprod(ZX, S2iZX))
  b2 <- as.numeric(V2 %*% crossprod(S2iZX, Zy))
  names(b2) <- colnames(X)
  e2 <- as.numeric(y - X %*% b2)
  V1r <- M1 %*% crossprod(S1iZX, S2 %*% S1iZX) %*% M1
  dn <- list(colnames(X), colnames(X))
  dimnames(M1) <- dimnames(V1r) <- dimnames(V2) <- dn
  out <- list(b1 = b1, b2 = b2, e1 = e1, e2 = e2, M1 = M1, V2 = V2, S1 = S1, S2 = S2,
              S1i = S1i, S2i = S2i, ZX = ZX)
  if (model == "onestep") {
    out$coefficients <- b1
    out$vcov_list <- list(robust = V1r)
  } else {
    # Windmeijer (2005) correction
    We <- crossprod(W, e2)
    S2iWe <- S2i %*% We
    D <- vapply(seq_len(ncol(X)), function(j) {
      Pk <- rowsum(W * X[, j], group, reorder = FALSE)
      dS <- -(crossprod(Pk, Q1) + crossprod(Q1, Pk))
      as.numeric(-V2 %*% crossprod(S2iZX, dS %*% S2iWe))
    }, numeric(ncol(X)))
    D <- matrix(D, ncol(X))
    Vw <- V2 + D %*% V2 + V2 %*% t(D) + D %*% V1r %*% t(D)
    Vw <- (Vw + t(Vw)) / 2
    dimnames(Vw) <- dn
    out$coefficients <- b2
    out$vcov_list <- list(robust = Vw)
  }
  out$vcov_list$onestep_robust <- V1r
  out
}

#' Arellano-Bond test of order `ord` (formula of plm::mtest())
#'
#' @keywords internal
#' @noRd
gmm_artest <- function(e, X, W, group, rowtype, m, R, N, ord, A, B, vv) {
  isd <- rowtype == "d"
  ed <- ifelse(isd, e, 0)
  # lag the differenced residuals within individual
  E <- matrix(ed, nrow = R)
  El <- matrix(0, R, N)
  El[(ord + 1):m, ] <- E[seq_len(m - ord), , drop = FALSE]
  el <- as.vector(El)
  s <- colSums(E * El)
  num <- sum(s)
  EVE <- sum(s^2)
  EX <- crossprod(el, X)
  ZVE <- crossprod(rowsum(W * ed, group, reorder = FALSE), s)
  ZX <- crossprod(W, X)
  denom <- EVE - 2 * EX %*% B %*% t(ZX) %*% A %*% ZVE + EX %*% vv %*% t(EX)
  # the variance estimate can be negative in small samples; the test is then
  # undefined (NA) rather than NaN with a warning
  stat <- if (is.finite(denom) && denom > 0) as.numeric(num / sqrt(denom)) else NA_real_
  data.frame(order = ord, statistic = stat,
             p.value = 2 * stats::pnorm(-abs(stat)))
}

#' Wald test of the coefficients other than time effects and the constant
#'
#' @keywords internal
#' @noRd
gmm_wald <- function(object, V) {
  keep <- setdiff(names(object$coefficients), object$timedummies)
  b <- object$coefficients[keep]
  stat <- as.numeric(crossprod(b, solve(V[keep, keep, drop = FALSE], b)))
  list(statistic = stat, df = length(keep),
       p.value = stats::pchisq(stat, length(keep), lower.tail = FALSE))
}

#' Covariance of a pangmm fit
#'
#' @param object a `"pangmm"` fit
#' @param type `NULL` (the covariance stored on the fit), `"robust"` (the
#'   one-step sandwich, or the Windmeijer-corrected two-step covariance) or
#'   `"classical"` (the homoskedastic one-step covariance, or the
#'   uncorrected two-step covariance); `"cluster"` is a synonym of
#'   `"robust"`
#' @param ... unused
#' @return the covariance matrix of the coefficients
#' @examples
#' set.seed(3)
#' N <- 200; TT <- 6
#' d <- data.frame(id = rep(seq_len(N), each = TT), t = rep(seq_len(TT), N))
#' mu <- rep(rnorm(N), each = TT)
#' d$y <- ave(mu + rnorm(N * TT), d$id, FUN = function(e) stats::filter(e, 0.4,
#'   method = "recursive"))
#' fit <- pangmm(y ~ lag(y, 1) | lag(y, 2:99), data = d, index = c("id", "t"),
#'               model = "twostep")
#' sqrt(diag(vcov(fit, type = "classical")))
#' sqrt(diag(vcov(fit, type = "robust")))
#' @export
vcov.pangmm <- function(object, type = NULL, ...) {
  if (is.null(type)) return(object$vcov)
  type <- match.arg(type, c("robust", "classical", "cluster"))
  if (type == "cluster") type <- "robust"
  object$vcov_list[[type]]
}

#' @export
coef.pangmm <- function(object, ...) object$coefficients

#' @export
nobs.pangmm <- function(object, ...) object$nobs

#' @export
residuals.pangmm <- function(object, ...) object$residuals

#' @export
confint.pangmm <- function(object, parm, level = 0.95, ...) {
  confint.panglm(object, parm, level = level, ...)
}

gmm_label <- function(x) {
  paste0(if (x$effect == "twoways") "two-way " else "individual ",
         if (x$model == "onestep") "one-step " else "two-step ",
         if (x$transformation == "d") "difference GMM (Arellano-Bond)" else
           "system GMM (Blundell-Bond)")
}

#' @export
print.pangmm <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("pangmm(", gmm_label(x), ")\n\n", sep = "")
  print(x$coefficients[setdiff(names(x$coefficients), x$timedummies)], digits = digits)
  cat("\nN =", x$n_groups, " T =", x$n_periods, " equations used =", x$nobs,
      " instruments =", x$n_instruments, "\n")
  invisible(x)
}

#' Summarize a pangmm fit
#'
#' @param object a `"pangmm"` fit
#' @param time.dummies logical: show the time effects (and the system-GMM
#'   constant) in the coefficient table
#' @param ... unused
#' @return an object of class `"summary.pangmm"`: a list with the
#'   coefficient table (`coefficients`: estimate, standard error, z value
#'   and p-value), the Arellano-Bond tests (`artest`), `sargan`, `hansen`,
#'   `wald` and the fit information (`call`, `label`, `vcov_type`, `nobs`,
#'   `n_groups`, `n_periods`, `n_instruments`).
#' @examples
#' set.seed(3)
#' N <- 200; TT <- 6
#' d <- data.frame(id = rep(seq_len(N), each = TT), t = rep(seq_len(TT), N))
#' mu <- rep(rnorm(N), each = TT)
#' d$y <- ave(mu + rnorm(N * TT), d$id, FUN = function(e) stats::filter(e, 0.4,
#'   method = "recursive"))
#' summary(pangmm(y ~ lag(y, 1) | lag(y, 2:99), data = d, index = c("id", "t")))
#' @export
summary.pangmm <- function(object, time.dummies = FALSE, ...) {
  coefs <- object$coefficients
  se <- sqrt(diag(object$vcov))
  z <- coefs / se
  table <- cbind(Estimate = coefs, `Std. Error` = se, `z value` = z,
                 `Pr(>|z|)` = 2 * stats::pnorm(-abs(z)))
  if (!time.dummies) table <- table[!rownames(table) %in% object$timedummies, , drop = FALSE]
  out <- list(call = object$call, label = gmm_label(object), coefficients = table,
              artest = object$artest, sargan = object$sargan, hansen = object$hansen,
              wald = object$wald, vcov_type = object$vcov_type, model = object$model,
              nobs = object$nobs, n_groups = object$n_groups, n_periods = object$n_periods,
              n_instruments = object$n_instruments)
  class(out) <- "summary.pangmm"
  out
}

#' @export
print.summary.pangmm <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nCall:\n")
  print(x$call)
  cat("\nEstimator:", x$label, "\n")
  se_lab <- if (x$vcov_type == "classical") "classical" else if (x$model == "twostep")
    "robust (Windmeijer)" else "robust"
  cat("N =", x$n_groups, " T =", x$n_periods, " equations used =", x$nobs,
      " instruments =", x$n_instruments, " vcov:", se_lab, "\n")
  cat("\nCoefficients:\n")
  stats::printCoefmat(x$coefficients, digits = digits, signif.stars = TRUE, has.Pvalue = TRUE)
  test_line <- function(name, t) {
    cat(name, ": chi2(", t$df, ") = ", format(t$statistic, digits = digits), ", p = ",
        format.pval(t$p.value, digits = digits), "\n", sep = "")
  }
  cat("\n")
  if (!is.null(x$artest)) {
    for (i in seq_len(nrow(x$artest))) {
      cat("Arellano-Bond AR(", x$artest$order[i], ") test: z = ",
          format(x$artest$statistic[i], digits = digits), ", p = ",
          format.pval(x$artest$p.value[i], digits = digits), "\n", sep = "")
    }
  }
  test_line("Sargan test (not robust)", x$sargan)
  test_line("Hansen J test (robust)", x$hansen)
  test_line("Wald test for coefficients", x$wald)
  if (x$n_instruments >= x$n_groups) {
    cat("\nNote: ", x$n_instruments, " instruments for ", x$n_groups, " individuals; ",
        "estimates are biased toward OLS. Use collapse = TRUE or a smaller maxlag.\n", sep = "")
  } else if (isTRUE(x$hansen$p.value > 0.99)) {
    cat("\nNote: a Hansen p-value this close to 1 often signals too many instruments ",
        "(Roodman 2009); compare with collapse = TRUE.\n", sep = "")
  }
  cat("\n")
  invisible(x)
}
