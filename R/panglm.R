#' Fit a generalized linear model for panel data
#'
#' Fast pooled, fixed-effects ("within") and random-effects panel GLMs,
#' backed by an Rcpp/RcppArmadillo/RcppParallel numerical core. See the
#' package vignette for the estimators currently implemented per family.
#'
#' @param formula a model formula, e.g. `y ~ x1 + x2`.
#'   `lag(x, k)`, `lead(x, k)` and `diff(x)` are taken within each
#'   individual on the grid of sorted distinct time values (needs
#'   `index = c(id, time)`); the first waves of each individual get `NA`
#'   and leave the sample, and `lag(x, 0:2)` gives several lags
#' @param data a data.frame containing the variables in `formula` and `index`
#' @param index length-1 or length-2 character vector identifying the
#'   individual (and optionally time) columns, e.g. `c("id", "time")`
#' @param model one of `"pooling"`, `"within"`, `"random"`, `"between"`.
#'   The between estimator (`family = "gaussian"` only) is least squares on
#'   the individual means of the outcome and the regressors, one row per
#'   individual, as in `plm(model = "between")`; it uses only the variation
#'   between individuals and is the cross-sectional counterpart of the
#'   within estimator. Its robust covariances (`"HC1"`, `"cluster"`) are the
#'   same, since each individual is one observation of the means regression.
#' @param family one of `"gaussian"`, `"poisson"`, `"binomial"`, `"probit"`,
#'   `"negbin"`, `"gamma"`, `"ologit"`, `"oprobit"`, `"beta"`, `"tobit"`,
#'   or a family spec from [gaussianfamily()] / [poissonfamily()] /
#'   [binomialfamily()] / [negbinfamily()] / [gammafamily()] /
#'   [orderedfamily()] / [betafamily()] / [tobitfamily()]
#' @param weights optional nonnegative observation (prior) weights: the
#'   name of a column of `data` or a numeric vector with one value per row
#'   of `data`, e.g. inverse probability weights for attrition.
#'   Nontrivial weights are supported for `model = "pooling"` (all
#'   families); the fit then maximizes the weighted (pseudo) log-likelihood
#'   and robust covariances use the weighted scores. For `model = "random"`
#'   the weights must be constant within individuals (one weight per
#'   individual, e.g. a baseline sampling weight) and the fit maximizes
#'   the weighted sum of the individual integrated log-likelihoods; the
#'   Gaussian model is then estimated by maximum likelihood instead of
#'   Swamy-Arora FGLS. Use `vcov = "cluster"` for design-based standard
#'   errors. Weighted fixed-effects fits are rejected.
#' @param offset optional model offset (the name of a column of `data` or a
#'   numeric vector with one value per row of `data`). Nonzero offsets, including offsets
#'   specified in `formula`, are rejected explicitly rather than ignored.
#' @param na.action missing-data action passed to [stats::model.frame()]
#' @param effect one of `"individual"` or `"twoways"` (individual + time
#'   fixed effects). `"twoways"` is implemented for `model = "within"`
#'   only: `family = "gaussian"` (exact alternating-projections
#'   demeaning), `"poisson"` and `"negbin"` (both via outer-IRLS /
#'   inner-weighted-FWL, the algorithm behind `fixest::feglm()`), and
#'   `"binomial"`, the two-way conditional (Chamberlain) logit: the
#'   individual effects are conditioned out exactly and the time effects
#'   are estimated as time dummies inside the conditional likelihood, as in
#'   `survival::clogit(y ~ x + factor(time) + strata(id), method = "exact")`.
#'   For the two-way logit `coef()` holds the slopes only; the time effects
#'   (first time level as reference, `NA` when not identified) and their
#'   standard errors are in `fit$time_effects` and `fit$time_effects_se`,
#'   and the classical and robust covariances of the slopes account for
#'   their estimation.
#' @param vcov one of `"classical"`, `"HC1"` (heteroskedasticity-robust) or
#'   `"cluster"` (cluster-robust, clustered by the panel individual by
#'   default). For `model = "random"` only `"classical"` and `"cluster"`
#'   are available; the cluster covariance uses the individual scores of
#'   the integrated likelihood, so clusters must contain whole individuals.
#'   See [vcov.panglm()] to recompute a different type after fitting.
#' @param cre logical. If `TRUE`, fit a correlated random effects
#'   (Mundlak) model: the individual means of every time-varying regressor
#'   are added to the design matrix, named `mean_<column>`. Available for
#'   `model = "pooling"` and `model = "random"` and every family. The slopes
#'   of the time-varying regressors then have a within interpretation, and
#'   the joint test of the mean terms is a robust Hausman-type test
#'   (reported by [summary.panglm()]).
#'   Defaults to `dynamic`.
#' @param dynamic logical. If `TRUE`, fit a dynamic model with the
#'   Wooldridge (2005) solution to the initial conditions problem: the
#'   outcome of the previous wave (`lag_<y>`) and the outcome of the first
#'   observed wave (`initial_<y>`) are added to the regressors, the first
#'   wave is dropped, and the individual means of the time-varying
#'   regressors are added (`cre = TRUE` by default; the lag is not
#'   averaged). The second element of `index` must be a numeric time
#'   variable; a lag is only defined when the previous wave `t - 1` is
#'   observed, other rows are dropped. Factor and logical outcomes (ordered
#'   and binary models) give factor or logical lags, so their partial
#'   effects are discrete changes. With `model = "random"` the coefficient
#'   of the lag measures true state dependence net of unobserved
#'   heterogeneity; [partialeffects()] gives it on the probability or mean
#'   scale. Available for `model = "pooling"` and `model = "random"`.
#' @param random one-sided formula for the random effects of
#'   `model = "random"`. The default `~ 1` is the random intercept. A
#'   formula with variables, e.g. `~ time`, adds correlated random slopes:
#'   `g(E[y_it | b_i]) = x_it' beta + z_it' b_i` with `z_it = (1, time)`
#'   and `b_i ~ N(0, Sigma)`, `Sigma` unstructured (Cholesky
#'   parameterized), the model of `lme4::lmer(y ~ x + (1 + time | id))`.
#'   The slope variables usually also appear in `formula` (their average
#'   slope). Available for `family = "gaussian"` (exact maximum likelihood,
#'   any number of slopes), and for `"poisson"`, `"negbin"` (NB2, the
#'   model of `glmmTMB::nbinom2()`, with `theta` stored as `negbin_theta`),
#'   `"gamma"` (log link, the model of `glmmTMB` with `Gamma(link = "log")`,
#'   shape stored as `shape`) and `"binomial"` (logit or probit), with any
#'   number of slopes, by adaptive Gauss-Hermite product quadrature with
#'   `R` nodes per dimension, `R^q` nodes for `q` random effects; `R = 1`
#'   is the Laplace approximation of `lme4::glmer()` and `glmmTMB`.
#'   With random slopes the random effects are always normal on the link
#'   scale, so for `"poisson"` and `"negbin"` the multiplicative effect
#'   `exp(z_it' b_i)` is lognormal. These are different models from the
#'   random-intercept fits of `random = ~ 1`: the random-intercept
#'   `"poisson"` model uses gamma mixing (the conjugate model of Hausman,
#'   Hall and Griliches, whose likelihood is in closed form, as Stata's
#'   `xtpoisson, re`) and the random-intercept `"negbin"` model is the
#'   beta-mixed negative binomial of Hausman, Hall and Griliches (Stata's
#'   `xtnbreg, re`), whereas with random slopes `"negbin"` is NB2 with a
#'   lognormal effect, as `glmmTMB::nbinom2()`. Adding a slope whose
#'   variance is estimated at zero therefore does not reproduce the
#'   `random = ~ 1` fit for these two families: the slope coefficients
#'   keep their conditional (rate ratio) meaning, but the intercept
#'   (`E[exp(b_i)] = exp(sigma^2 / 2)` under lognormal mixing against 1
#'   under gamma mixing), the variance parameters and the log-likelihood
#'   differ. For `"binomial"` and `"gamma"` the random-intercept models are
#'   already normal on the link scale and are the zero-slope special case.
#'   The fit stores `varcomp` (the
#'   covariance matrix `Sigma`), `ranef` (posterior modes, the BLUPs in
#'   the Gaussian model, one row per individual) and, for the Gaussian
#'   model, the residual variance `sigma_v2`. Weights, `dynamic = TRUE`
#'   and the other families are not supported with random slopes;
#'   `cre = TRUE` is. Response predictions and [partialeffects()] are
#'   population averages over the random effects.
#' @param R number of Gauss-Hermite quadrature nodes (random-effects
#'   binomial and gamma models only). With random slopes, the number of
#'   nodes per dimension of the adaptive product rule. When `R` is not
#'   given it is 7 for one slope, 5 for two, 3 for three to five and 1
#'   (Laplace) beyond, so that `R^q` stays moderate; a warning is given
#'   when `R^q` exceeds 2000. `R = 1` gives the Laplace approximation.
#' @param maxit maximum IRLS/Newton iterations
#' @param tol convergence tolerance
#' @return an object of class `"panglm"` (with the additional class
#'   `"panordered"` for the ordered families): a list with `coefficients`,
#'   `vcov`, `vcov_type`, `fitted.values`, `linear.predictors`, `loglik`, the
#'   variance components of random-effects fits (`sigma_mu2`, `sigma_v2`,
#'   `theta`), the panel structure, `nobs`, `n_groups`, `call` and `formula`.
#'   Methods: `print`, `summary`, `coef`, `vcov`, `confint`, `predict`,
#'   `fitted`, `residuals`, `logLik`, `nobs`, `model.matrix`, `plot`,
#'   [partialeffects()], [effectcurve()], `tidy` and `glance`.
#' @family panel models
#' @seealso [panglm-families] for the outcome families, [summary.panglm()],
#'   [vcov.panglm()], [predict.panglm()] and [partialeffects()] for the fitted
#'   model, [tidy-methods] for `tidy()` and `glance()`.
#' @section Validation:
#' Each estimator is checked against an established implementation of the same
#' model, on the same data, in the package's validation suite on GitHub and in
#' the package tests:
#' * pooled models: `stats::lm()`, `stats::glm()`, `MASS::glm.nb()`,
#'   `plm::plm()` and fixest, with `sandwich` robust and cluster covariances;
#' * fixed effects: `plm::plm()`, `fixest::feols()`, `fixest::fepois()`,
#'   `fixest::fenegbin()`, `pglm::pglm()`, `glm()` with individual dummies, and
#'   `survival::clogit()` for the conditional and two-way logit;
#' * between estimator: `plm::plm(model = "between")`;
#' * random effects: `plm::plm()` and `pglm::pglm()` (random intercept),
#'   `lme4::lmer()`, `lme4::glmer()` and glmmTMB (random intercept and slopes),
#'   and GLMMadaptive (adaptive quadrature with random slopes);
#' * ordered outcomes: `MASS::polr()`, `ordinal::clm()`, `ordinal::clmm()`, and
#'   a blow-up and cluster fit with `survival::clogit()` for the fixed-effects
#'   ordered logit;
#' * Tobit: `AER::tobit()` and censReg; beta: `betareg::betareg()` and glmmTMB;
#'   gamma: `stats::glm()` and glmmTMB.
#' @examples
#' data(copd)
#'
#' # pooled OLS on lung function (ignores the panel structure)
#' fit_pool <- panglm(fev1 ~ treatment + age + smoker + crp, data = copd,
#'                     index = c("id", "visit"), model = "pooling", family = "gaussian")
#' summary(fit_pool)
#'
#' # fixed-effects ("within") estimator: only the time-varying regressor
#' # (crp) survives demeaning; time-invariant covariates like treatment
#' # or age are not identified by a fixed-effects estimator
#' fit_fe <- panglm(fev1 ~ crp, data = copd,
#'                   index = c("id", "visit"), model = "within", family = "gaussian")
#' coef(fit_fe)
#'
#' # random-effects (Swamy-Arora) estimator
#' fit_re <- panglm(fev1 ~ treatment + age + smoker + crp, data = copd,
#'                   index = c("id", "visit"), model = "random", family = "gaussian")
#' coef(fit_re)
#'
#' # fixed-effects Poisson on exacerbation counts
#' fit_pois <- panglm(exacerbations ~ crp, data = copd,
#'                     index = c("id", "visit"), model = "within", family = "poisson")
#' coef(fit_pois)
#'
#' # correlated random effects (Mundlak) probit for hospitalization
#' fit_cre <- panglm(hospitalized ~ treatment + age + crp, data = copd,
#'                   index = c("id", "visit"), model = "random",
#'                   family = "probit", cre = TRUE)
#' summary(fit_cre)
#'
#' # patient trajectories: random intercept and random slope of FEV1
#' # decline per visit, correlated
#' fit_rs <- panglm(fev1 ~ visit + treatment + age + smoker, data = copd,
#'                  index = c("id", "visit"), model = "random",
#'                  random = ~ visit)
#' summary(fit_rs)
#' fit_rs$varcomp
#'
#' # two-way conditional logit: smoking across survey waves, with year
#' # shocks (e.g. a tobacco tax) absorbed by the time effects
#' set.seed(1)
#' smk <- expand.grid(wave = 1:5, id = 1:200)
#' smk$price <- rnorm(nrow(smk)) + 0.2 * smk$wave
#' smk$smoking <- rbinom(nrow(smk), 1, plogis(rnorm(200)[smk$id] - 0.5 * smk$price +
#'                                                c(0, 0.2, -0.3, -0.6, -0.8)[smk$wave]))
#' fit_2w <- panglm(smoking ~ price, data = smk, index = c("id", "wave"),
#'                  model = "within", family = "binomial", effect = "twoways",
#'                  vcov = "cluster")
#' summary(fit_2w)
#' fit_2w$time_effects
#'
#' # nonlinear effect through a natural spline: linearity test and curve
#' fit_ns <- panglm(fev1 ~ splines::ns(crp, 3), data = copd,
#'                  index = c("id", "visit"), model = "within", vcov = "cluster")
#' lineartest(fit_ns, "crp")
#' plot(effectcurve(fit_ns, "crp"))
#' @export
panglm <- function(formula, data, index,
                    model = c("pooling", "within", "random", "between"),
                    family = "gaussian",
                    weights = NULL, offset = NULL,
                    na.action = stats::na.omit,
                    effect = c("individual", "twoways"),
                    vcov = c("classical", "HC1", "cluster"),
                    cre = dynamic, dynamic = FALSE, random = ~ 1,
                    R = 21, maxit = 100, tol = 1e-10) {
  model <- match.arg(model)
  effect <- match.arg(effect)
  vcov <- match.arg(vcov)
  family <- resolve_family(family)
  weights <- data_column(weights, data, "weights")
  offset <- data_column(offset, data, "offset")
  slopes <- random_slopes_requested(random)
  if (slopes) {
    if (model != "random") {
      stop("'random' with slope variables applies to model = 'random'", call. = FALSE)
    }
    if (!family$family %in% c("gaussian", "poisson", "negbin", "gamma", "binomial")) {
      stop("random slopes are available for family = 'gaussian', 'poisson', 'negbin', ",
           "'gamma' and 'binomial' only", call. = FALSE)
    }
    if (!is.null(weights)) {
      stop("observation weights are not supported with random slopes", call. = FALSE)
    }
    if (isTRUE(dynamic)) {
      stop("dynamic = TRUE is not supported with random slopes", call. = FALSE)
    }
    if (missing(R)) {
      R <- NULL # chosen from the number of random effects
    } else {
      if (!is.numeric(R) || length(R) != 1L || R < 1) {
        stop("'R' must be a positive integer", call. = FALSE)
      }
      R <- as.integer(R)
    }
  }
  is_ordered <- identical(family$family, "ordered")
  if (is_ordered && model == "within" && family$link != "logit") {
    stop("the fixed-effects ordered model is available for the logit link only ",
         "(blow-up and cluster); use model = 'random' with cre = TRUE for the probit",
         call. = FALSE)
  }
  if (family$family == "binomial" && model == "within" && family$link != "logit") {
    stop("the fixed-effects binomial model is the conditional logit; there is no ",
         "consistent fixed-effects probit, use model = 'random' with cre = TRUE",
         call. = FALSE)
  }
  if (family$family == "beta" && model == "within") {
    stop("there is no fixed-effects beta regression; use model = 'random' with ",
         "cre = TRUE", call. = FALSE)
  }
  if (family$family == "tobit" && !model %in% c("pooling", "random")) {
    stop("the Tobit model is available for model = 'pooling' and model = 'random'; ",
         "there is no consistent fixed-effects Tobit, use model = 'random' with cre = TRUE",
         call. = FALSE)
  }
  if (family$family == "tobit" && effect != "individual") {
    stop("the Tobit model uses effect = 'individual'", call. = FALSE)
  }
  if (model == "between") {
    if (family$family != "gaussian") {
      stop("model = 'between' is available for family = 'gaussian' only", call. = FALSE)
    }
    if (isTRUE(dynamic) || isTRUE(cre)) {
      stop("model = 'between' cannot be combined with dynamic or cre", call. = FALSE)
    }
    if (effect != "individual") {
      stop("model = 'between' uses individual means; effect must be 'individual'",
           call. = FALSE)
    }
    if (!is.null(weights)) {
      stop("observation weights are not supported for model = 'between'", call. = FALSE)
    }
  }
  if (!is.logical(dynamic) || length(dynamic) != 1L || is.na(dynamic)) {
    stop("'dynamic' must be TRUE or FALSE", call. = FALSE)
  }
  if (dynamic && model == "within") {
    stop("dynamic = TRUE applies to model = 'pooling' or 'random': the within ",
         "estimator of a dynamic model is inconsistent for short panels (Nickell ",
         "bias); use a GMM estimator such as plm::pgmm() instead", call. = FALSE)
  }
  if (!is.logical(cre) || length(cre) != 1L || is.na(cre)) {
    stop("'cre' must be TRUE or FALSE", call. = FALSE)
  }
  if (cre && model == "within") {
    stop("cre = TRUE applies to model = 'pooling' or 'random'; the within ",
         "estimator already removes the individual effect", call. = FALSE)
  }

  if (missing(index)) stop("'index' is required, e.g. index = c(\"id\", \"time\")", call. = FALSE)
  if (!is.character(index) || !length(index) %in% 1:2 || anyDuplicated(index)) {
    stop("'index' must contain one or two distinct column names", call. = FALSE)
  }
  missing_index <- setdiff(index, names(data))
  if (length(missing_index)) {
    stop("index column(s) not found in data: ", paste(missing_index, collapse = ", "),
         call. = FALSE)
  }
  if (effect == "twoways" &&
      !(model == "within" && family$family %in% c("gaussian", "poisson", "negbin", "binomial"))) {
    stop("effect = 'twoways' is currently only implemented for ",
         "model = 'within', family = 'gaussian', 'poisson', 'negbin', or 'binomial'",
         call. = FALSE)
  }
  if (effect == "twoways" && family$family == "binomial" && family$link != "logit") {
    stop("effect = 'twoways' with family = 'binomial' is the conditional logit and ",
         "is only implemented for the logit link", call. = FALSE)
  }

  formula <- panel_guard_formula(formula)
  dynamic_info <- NULL
  frame_data <- data
  frame_rows <- NULL
  if (dynamic) {
    dyn <- dynamic_panel_data(formula, data, index, weights, offset)
    formula <- dyn$formula; data <- dyn$data; weights <- dyn$weights; offset <- dyn$offset
    dynamic_info <- dyn$info
    frame_data <- dyn$full
    frame_rows <- dyn$keep
  }

  panel_calls <- has_panel_calls(formula)
  # lag(), lead() and diff() are taken within individuals on all rows
  mf <- panel_model_frame(formula, frame_data, index, na.action = na.action,
                          rows = frame_rows)
  model_rows <- match(rownames(mf), rownames(data))
  if (anyNA(model_rows)) stop("could not align model-frame rows with data", call. = FALSE)

  used_weights <- NULL
  if (!is.null(weights)) {
    if (!is.numeric(weights) || length(weights) != nrow(data) ||
        any(!is.finite(weights)) || any(weights < 0)) {
      stop("'weights' must be a column name or a finite nonnegative numeric vector with nrow(data) values",
           call. = FALSE)
    }
    used_weights <- weights[model_rows]
    if (all(used_weights == 1)) {
      used_weights <- NULL
    } else if (model == "within" || (model == "random" && !weighted_random_ok(family))) {
      stop("nontrivial observation weights are only supported for model = 'pooling'",
           " and model = 'random'", call. = FALSE)
    } else if (!any(used_weights > 0)) {
      stop("at least one observation weight must be positive", call. = FALSE)
    }
  }
  formula_offset <- stats::model.offset(mf)
  if (!is.null(offset)) {
    if (!is.numeric(offset) || length(offset) != nrow(data) || any(!is.finite(offset))) {
      stop("'offset' must be a column name or a finite numeric vector with nrow(data) values",
           call. = FALSE)
    }
    formula_offset <- (formula_offset %||% 0) + offset[model_rows]
  }
  if (!is.null(formula_offset) && any(formula_offset != 0)) {
    stop("nonzero offsets are not supported", call. = FALSE)
  }
  y <- stats::model.response(mf)
  ordered_levels <- NULL
  if (is_ordered) {
    coded <- ordered_response(y)
    y <- coded$code
    ordered_levels <- coded$levels
  }
  if (family$family == "binomial") y <- binomial_response_to_numeric(y)
  if (!is.numeric(y) || any(!is.finite(y))) {
    stop("the response must be finite and numeric", call. = FALSE)
  }
  if (family$family %in% c("poisson", "negbin") && any(y < 0)) {
    stop("count-family responses must be nonnegative", call. = FALSE)
  }
  if (family$family == "gamma" && any(y[if (is.null(used_weights)) TRUE else used_weights > 0] <= 0)) {
    stop("family = 'gamma' requires a strictly positive response", call. = FALSE)
  }
  if (family$family == "beta" && any(y <= 0 | y >= 1)) {
    stop("family = 'beta' requires a response strictly between 0 and 1; ",
         "transform boundary values first, e.g. (y * (n - 1) + 0.5) / n", call. = FALSE)
  }
  X_full <- stats::model.matrix(formula, data = mf)
  terms_object <- stats::terms(mf)
  dynamic_columns <- if (dynamic) {
    involved <- attr(terms_object, "factors")[dynamic_info$lag, ] > 0
    colnames(X_full)[attr(X_full, "assign") %in% which(involved)]
  } else character(0)
  xlevels <- stats::.getXlevels(terms_object, mf)
  contrasts <- attr(X_full, "contrasts")

  panel_cols <- data[model_rows, index, drop = FALSE]
  if (anyNA(panel_cols)) stop("panel index columns must not contain missing values", call. = FALSE)
  panel <- build_panel_index(cbind(panel_cols, .panglm_row = seq_len(nrow(mf))), index)
  ord <- panel$data$.panglm_row

  y <- y[ord]
  X_full <- X_full[ord, , drop = FALSE]
  if (slopes) stop_panel_calls(random, what = "the 'random' formula")
  Z_random <- if (slopes) slopes_z(random, data[model_rows[ord], , drop = FALSE]) else NULL
  prior_weights <- if (is.null(used_weights)) NULL else used_weights[ord]
  group_start <- panel$group_start
  group_size <- panel$group_size
  time_id <- if (length(index) > 1) panel$data[[index[2]]] else NULL

  has_intercept <- "(Intercept)" %in% colnames(X_full)
  X_noint <- if (has_intercept) X_full[, setdiff(colnames(X_full), "(Intercept)"), drop = FALSE] else X_full

  cre_columns <- character(0)
  if (cre) {
    # the lagged outcome is not averaged (Wooldridge 2005)
    cre_means <- cre_mean_columns(X_noint[, setdiff(colnames(X_noint), dynamic_columns),
                                          drop = FALSE], group_start, group_size)
    cre_columns <- colnames(cre_means)
    X_full <- cbind(X_full, cre_means)
  }
  if (is_ordered && model != "within") {
    # the cutpoints play the role of the intercept
    X_full <- X_full[, setdiff(colnames(X_full), "(Intercept)"), drop = FALSE]
  }

  candidate_X <- if (model %in% c("pooling", "random", "between")) X_full else X_noint
  estimable <- if (is_ordered && model != "within") {
    setdiff(panglm_estimable_columns(cbind(`(Intercept)` = 1, candidate_X), model, effect,
                                     group_start, group_size, time_id), 1L) - 1L
  } else {
    panglm_estimable_columns(candidate_X, model, effect, group_start, group_size, time_id)
  }
  if (!length(estimable)) stop("the model contains no estimable covariates", call. = FALSE)
  aliased <- setdiff(colnames(candidate_X), colnames(candidate_X)[estimable])
  if (length(aliased)) {
    warning("rank-deficient model matrix: dropping unidentified column(s): ",
            paste(aliased, collapse = ", "), call. = FALSE)
  }
  if (model %in% c("pooling", "random", "between")) {
    X_full <- X_full[, estimable, drop = FALSE]
  } else {
    X_noint <- X_noint[, estimable, drop = FALSE]
  }

  fit <- if (is_ordered) {
    fit_ordered(if (model == "within") X_noint else X_full, y, ordered_levels, family$link,
                model, group_start, group_size, prior_weights, R, maxit, tol)
  } else if (family$family == "beta") {
    fit_beta(X_full, y, family, model, group_start, group_size, prior_weights, R, maxit, tol)
  } else if (family$family == "tobit") {
    fit_tobit(X_full, y, family, model, group_start, group_size, prior_weights, R, maxit, tol)
  } else switch(model,
    pooling = fit_pooled(X_full, y, family, maxit, tol, prior_weights),
    between = fit_between(X_full, y, group_start, group_size),
    within  = if (effect == "twoways" && family$family == "gaussian") {
      fit_within_twoways_gaussian(X_noint, y, group_start, group_size, time_id, maxit, tol)
    } else if (effect == "twoways" && family$family == "poisson") {
      fit_within_twoways_poisson(X_noint, y, group_start, group_size, time_id, maxit, tol)
    } else if (effect == "twoways" && family$family == "negbin") {
      fit_within_twoways_negbin(X_noint, y, group_start, group_size, time_id, maxit, tol)
    } else if (effect == "twoways" && family$family == "binomial") {
      fit_within_twoways_binomial(X_noint, y, group_start, group_size, time_id, maxit, tol)
    } else {
      fit_within(X_noint, y, family, group_start, group_size, maxit, tol)
    },
    random  = if (slopes) {
      fit_random_slopes(X_full, y, Z_random, family, group_size, R, maxit, tol)
    } else if (is.null(prior_weights)) {
      fit_random(X_full, y, family, group_start, group_size, R, maxit, tol)
    } else {
      fit_random_weighted(X_full, y, family, group_start, group_size, R, maxit, tol,
                          prior_weights)
    }
  )

  fit$X <- if (model %in% c("pooling", "random", "between")) X_full else X_noint
  fit$y <- y
  fit$weights <- prior_weights
  fit$cre <- cre
  fit$dynamic <- dynamic_info
  fit$cre_columns <- intersect(cre_columns, colnames(fit$X))
  data_vars <- unique(c(intersect(c(all.vars(terms_object), all.vars(random)), names(data)),
                        index))
  if (slopes) {
    fit$random <- random
    fit$Z <- Z_random
    rownames(fit$ranef) <- as.character(panel$group_id)
  }
  fit$data <- data[model_rows[ord], data_vars, drop = FALSE]
  if (panel_calls) {
    fit$panel_history <- frame_data[, intersect(data_vars, names(frame_data)), drop = FALSE]
  }
  fit$group_start <- group_start
  fit$group_size <- group_size
  fit$time_id <- time_id
  fit$cluster_id <- rep(panel$group_id, group_size)
  fit$group_levels <- panel$group_id
  fit$time_levels <- if (is.null(time_id)) NULL else sort(unique(time_id))
  fit$panel_order <- ord
  fit$inverse_order <- order(ord)
  fit$terms <- terms_object
  fit$xlevels <- xlevels
  fit$contrasts <- contrasts
  fit$aliased <- aliased
  fit$na.action <- attr(mf, "na.action")
  if (!is.null(fit$individual_effects)) {
    effect_index <- if (!is.null(names(fit$individual_effects))) {
      suppressWarnings(as.integer(names(fit$individual_effects)))
    } else {
      seq_along(fit$individual_effects)
    }
    if (all(is.finite(effect_index)) && all(effect_index %in% seq_along(panel$group_id))) {
      names(fit$individual_effects) <- as.character(panel$group_id[effect_index])
    }
  }
  if (!is.null(fit$time_effects) && length(fit$time_effects) == length(fit$time_levels)) {
    names(fit$time_effects) <- as.character(fit$time_levels)
  }
  if (is.null(fit$linear.predictors)) {
    fit$linear.predictors <- if (model == "random") {
      as.numeric(X_full %*% fit$coefficients)
    } else if (family$family == "gaussian") {
      fit$fitted.values
    } else if (!is.null(fit$fitted.values)) {
      linkfun_r(fit$fitted.values, family$link_id)
    } else {
      NULL
    }
  }
  fit$call <- match.call()
  fit$formula <- formula
  fit$model <- model
  fit$effect <- effect
  fit$family <- family
  fit$index <- index
  fit$nobs <- length(y)
  fit$n_groups <- panel$n_groups
  fit$npar <- fit$npar %||% panglm_parameter_count(
    fit, model, family$family, effect, panel$n_groups,
    if (is.null(time_id)) 0L else length(unique(time_id))
  )
  fit$nobs_likelihood <- if (!is.null(fit$keep_rows)) sum(fit$keep_rows) else length(y)
  fit$vcov_type <- "classical"
  class(fit) <- c(if (is_ordered) "panordered", "panglm")

  # The model-based covariance is kept so that vcov(fit, type = "classical")
  # stays available after a robust covariance replaces fit$vcov.
  fit$vcov_classical <- fit$vcov
  if (is_ordered && model == "within" && vcov == "classical") vcov <- "cluster"
  if (vcov != "classical") {
    fit$vcov <- stats::vcov(fit, type = vcov)
    fit$vcov_type <- vcov
  }
  fit
}

panglm_estimable_columns <- function(X, model, effect, group_start, group_size, time) {
  if (!ncol(X)) stop("the model contains no estimable covariates", call. = FALSE)
  transformed <- if (model == "between") {
    group_means_cpp(X, rep(0, nrow(X)), group_start, group_size)$Xbar
  } else if (model != "within") {
    X
  } else if (effect == "twoways") {
    demean_twoway(X, rep(0, nrow(X)), group_start, group_size, time)$X
  } else {
    within_demean_cpp(X, rep(0, nrow(X)), group_start, group_size)$X
  }
  # A column removed by the transformation up to rounding (e.g. a time
  # trend under two-way demeaning) keeps a norm of order 1e-13 and would
  # pass qr()'s test, which is relative to the transformed column itself,
  # so the norm is first compared with that of the original column.
  raw_norm <- sqrt(colSums(X^2))
  kept_norm <- sqrt(colSums(transformed^2))
  candidate <- which(kept_norm > 1e-7 * raw_norm)
  if (!length(candidate)) {
    stop("the model contains no covariates identified by the requested estimator",
         call. = FALSE)
  }
  decomposition <- qr(transformed[, candidate, drop = FALSE], tol = 1e-7)
  if (decomposition$rank == 0L) {
    stop("the model contains no covariates identified by the requested estimator",
         call. = FALSE)
  }
  sort(candidate[decomposition$pivot[seq_len(decomposition$rank)]])
}

panglm_parameter_count <- function(fit, model, family, effect, n_groups, n_time) {
  k <- length(fit$coefficients)
  if (model == "pooling") {
    return(k + as.integer(family %in% c("gaussian", "negbin", "gamma")))
  }
  if (model == "between") return(k + 1L)
  if (model == "random") {
    extra <- switch(family, gaussian = 2L, poisson = 1L, negbin = 2L,
                    binomial = 1L, gamma = 2L, 0L)
    return(k + extra)
  }
  if (effect == "twoways") {
    fixed_effect_rank <- n_groups + n_time - 1L
    return(k + fixed_effect_rank + as.integer(family %in% c("gaussian", "negbin")))
  }
  if (family %in% c("poisson", "binomial")) return(k)
  used_groups <- fit$n_used_groups %||% n_groups
  k + used_groups + as.integer(family %in% c("gaussian", "negbin", "gamma"))
}

#' Between estimator: least squares on the individual means
#'
#' One row per individual (unweighted, as in plm). The log-likelihood is the
#' Gaussian likelihood of the means regression; residuals and fitted values
#' are returned at the observation level for predict() and residuals().
#'
#' @keywords internal
#' @noRd
fit_between <- function(X, y, group_start, group_size) {
  gm <- group_means_cpp(X, y, group_start, group_size)
  Xb <- gm$Xbar
  yb <- as.numeric(gm$ybar)
  colnames(Xb) <- colnames(X)
  G <- nrow(Xb); k <- ncol(Xb)
  if (G <= k) stop("model = 'between' needs more individuals than coefficients", call. = FALSE)
  decomposition <- qr(Xb)
  coefs <- as.numeric(qr.coef(decomposition, yb))
  names(coefs) <- colnames(X)
  resid_means <- yb - as.numeric(Xb %*% coefs)
  rss <- sum(resid_means^2)
  sigma2 <- rss / (G - k)
  bread <- chol2inv(qr.R(decomposition))
  bread <- bread[order(decomposition$pivot), order(decomposition$pivot), drop = FALSE]
  dimnames(bread) <- list(colnames(X), colnames(X))
  fitted <- as.numeric(X %*% coefs)
  list(coefficients = coefs, vcov = bread * sigma2, bread = bread,
       fitted.values = fitted, linear.predictors = fitted,
       loglik = -0.5 * G * (log(2 * pi) + log(rss / G) + 1),
       dispersion = sigma2, df.residual = G - k, iterations = 1L,
       between_X = Xb, between_resid = resid_means)
}

fit_pooled <- function(X, y, family, maxit, tol, weights = NULL) {
  n_pos <- if (is.null(weights)) nrow(X) else sum(weights > 0)
  if (family$family == "negbin") {
    res <- negbin_irls_fit_cpp(X, y, maxit, tol, weights)
    coefs <- as.numeric(res$coefficients)
    names(coefs) <- colnames(X)
    vcov <- res$vcov_unscaled
    dimnames(vcov) <- list(colnames(X), colnames(X))
    fitted <- as.numeric(exp(X %*% res$coefficients))
    return(list(coefficients = coefs, vcov = vcov, bread = vcov, fitted.values = fitted,
                loglik = res$loglik, theta = res$theta, dispersion = 1,
                df.residual = n_pos - ncol(X) - 1, iterations = res$iterations))
  }
  res <- irls_fit_cpp(X, y, family$family_id, family$link_id, maxit, tol, weights)
  vcov <- res$vcov_unscaled * res$dispersion
  bread <- res$vcov_unscaled
  coefs <- as.numeric(res$coefficients)
  names(coefs) <- colnames(X)
  dimnames(vcov) <- dimnames(bread) <- list(colnames(X), colnames(X))
  fitted <- as.numeric(res$fitted.values)
  loglik <- res$loglik
  if (family$family == "gamma") loglik <- gamma_loglik(y, fitted, weights)
  list(coefficients = coefs, vcov = vcov, bread = bread, fitted.values = fitted,
       loglik = loglik, dispersion = res$dispersion, df.residual = n_pos - ncol(X),
       iterations = res$iterations)
}

#' Gamma log-likelihood with the shape estimated from the deviance
#'
#' Uses the convention of `stats::Gamma()$aic`: dispersion = deviance / sum(w),
#' so that logLik() and AIC() agree with `stats::glm(family = Gamma)`.
#'
#' @keywords internal
#' @noRd
gamma_loglik <- function(y, mu, weights = NULL) {
  w <- if (is.null(weights)) rep(1, length(y)) else weights
  keep <- w > 0 & is.finite(mu)
  y <- y[keep]; mu <- mu[keep]; w <- w[keep]
  dev <- 2 * sum(w * (-log(y / mu) + (y - mu) / mu))
  disp <- dev / sum(w)
  sum(w * stats::dgamma(y, shape = 1 / disp, scale = mu * disp, log = TRUE))
}

#' Individual means of the time-varying columns (Mundlak device)
#'
#' @keywords internal
#' @noRd
cre_mean_columns <- function(X, group_start, group_size) {
  if (!ncol(X)) stop("cre = TRUE needs at least one regressor", call. = FALSE)
  deviation <- within_demean_cpp(X, rep(0, nrow(X)), group_start, group_size)$X
  scale <- pmax(1, apply(abs(X), 2, max))
  varying <- apply(abs(deviation), 2, max) > 1e-8 * scale
  if (!any(varying)) {
    stop("cre = TRUE needs at least one regressor that varies within individuals",
         call. = FALSE)
  }
  means <- X[, varying, drop = FALSE] - deviation[, varying, drop = FALSE]
  colnames(means) <- paste0("mean_", colnames(X)[varying])
  means
}

fit_within <- function(X, y, family, group_start, group_size, maxit, tol) {
  if (family$family == "gaussian") {
    dm <- within_demean_cpp(X, y, group_start, group_size)
    res <- irls_fit_cpp(dm$X, dm$y, family$family_id, family$link_id, maxit, tol)
    n <- nrow(X); k <- ncol(X); G <- length(group_start)
    df_resid <- n - G - k
    resid <- dm$y - dm$X %*% res$coefficients
    sigma2 <- sum(resid^2) / max(1, df_resid)
    vcov <- res$vcov_unscaled * sigma2
    bread <- res$vcov_unscaled
    coefs <- as.numeric(res$coefficients)
    names(coefs) <- colnames(X)
    dimnames(vcov) <- dimnames(bread) <- list(colnames(X), colnames(X))
    gm <- group_means_cpp(X, y, group_start, group_size)
    alpha_i <- gm$ybar - as.numeric(gm$Xbar %*% res$coefficients)
    alpha_obs <- rep(alpha_i, group_size)
    fitted <- alpha_obs + as.numeric(X %*% res$coefficients)

    ## By the Frisch-Waugh-Lovell theorem, the within-demeaned residuals
    ## (`resid` above) are algebraically identical to the residuals of the
    ## full least-squares-dummy-variable regression (y on X plus one
    ## intercept per group) - so the Gaussian log-likelihood of that full
    ## model (npar = k regressors + G group intercepts + 1 residual
    ## variance, exactly what panglm_parameter_count() already reports for
    ## this branch) is the standard OLS ML log-likelihood evaluated on
    ## `resid` at n = total observations, not the demeaned df:
    ##   logLik = -n/2 * (log(2*pi) + log(RSS/n) + 1)
    ## (same formula as stats:::logLik.lm()). Recomputed here from
    ## `fitted`/`y` directly rather than assumed, so it is exact even if
    ## `resid` and the LSDV residual ever diverge numerically.
    full_resid <- y - fitted
    rss <- sum(full_resid^2)
    loglik <- if (rss > 0) -0.5 * n * (log(2 * pi) + log(rss / n) + 1) else NA_real_

    return(list(coefficients = coefs, vcov = vcov, bread = bread, fitted.values = fitted,
                individual_effects = alpha_i,
                loglik = loglik, dispersion = sigma2, df.residual = df_resid,
                iterations = res$iterations))
  }
  if (family$family == "poisson") {
    scale <- apply(X, 2, function(col) { s <- stats::sd(col); if (s > 0) s else 1 })
    Xs <- sweep(X, 2, scale, "/")
    res <- within_poisson_fit_cpp(Xs, y, group_start, group_size, maxit, tol)
    coefs <- as.numeric(res$coefficients) / scale
    names(coefs) <- colnames(X)
    vcov <- res$vcov_unscaled / outer(scale, scale)
    dimnames(vcov) <- list(colnames(X), colnames(X))

    # Closed-form group intercept for the conditional Poisson MLE:
    # alpha_i = log(sum(y_i) / sum(exp(x_it'beta))), so mu_it = (Yi/Li) * exp(x_it'beta).
    eta <- as.numeric(X %*% coefs)
    lit <- exp(eta)
    group <- rep(seq_along(group_size), group_size)
    Li <- as.numeric(rowsum(lit, group))
    Yi <- as.numeric(rowsum(y, group))
    fitted <- (Yi / Li)[group] * lit
    alpha_i <- log(Yi / Li)

    return(list(coefficients = coefs, vcov = vcov, bread = vcov, fitted.values = fitted,
                individual_effects = alpha_i,
                loglik = res$loglik, dispersion = 1, df.residual = nrow(X) - ncol(X) - length(group_start),
                iterations = res$iterations))
  }
  if (family$family == "negbin") return(fit_within_negbin_dummy(X, y, group_start, group_size, maxit, tol))
  if (family$family == "binomial") return(fit_within_binomial(X, y, group_start, group_size, maxit, tol))
  if (family$family == "gamma") return(fit_within_gamma(X, y, group_start, group_size, maxit, tol))
  stop("model = 'within' is not yet implemented for family '", family$family, "'", call. = FALSE)
}

#' Unconditional dummy-variable fixed-effects negative binomial
#'
#' Augments the design matrix with one dummy column per panel individual
#' and jointly estimates (covariate slopes, individual intercepts, NB2
#' dispersion) by NB2 maximum likelihood. This is the unconditional
#' dummy-variable specification discussed by Allison and Waterman (2002),
#' not the Hausman-Hall-Griliches conditional negative binomial model. It
#' represents individual heterogeneity through explicit intercepts, but it
#' remains a nonlinear fixed-effects estimator. With short panels, finite-T
#' incidental-parameter bias can therefore be relevant. The augmented
#' design also grows with the number of individuals, so the estimator is
#' intended for small-to-moderate panels.
#'
#' @keywords internal
#' @noRd
fit_within_negbin_dummy <- function(X, y, group_start, group_size, maxit, tol) {
  n0 <- nrow(X); k <- ncol(X); G0 <- length(group_size)
  group <- rep(seq_along(group_size), group_size)
  colnames_X <- colnames(X)

  # Groups with an all-zero outcome have an unbounded dummy-variable MLE
  # (the intercept diverges to -Inf under the log link), the same boundary
  # case fit_within_binomial() screens out for the conditional logit. Left
  # in, the divergent direction can contaminate the shared beta via the
  # joint Newton/IRLS solve, not just that group's own intercept.
  zero_group <- vapply(seq_len(G0), function(g) all(y[group == g] == 0), logical(1))
  n_dropped <- sum(zero_group)
  kept_groups <- which(!zero_group)
  keep_rows <- !zero_group[group]
  if (n_dropped > 0) {
    message(n_dropped, " group(s) with all-zero outcomes dropped -- their fixed effect ",
            "is unbounded under the dummy-variable NB2 MLE and carries no information ",
            "about the shared slope.")
    X <- X[keep_rows, , drop = FALSE]
    y <- y[keep_rows]
    group <- match(group[keep_rows], kept_groups)
  }
  n <- length(y); G <- length(kept_groups)

  dummies <- matrix(0, n, G)
  dummies[cbind(seq_len(n), group)] <- 1
  X_aug <- cbind(X, dummies)

  res <- negbin_irls_fit_cpp(X_aug, y, maxit, tol)
  coefs <- as.numeric(res$coefficients[seq_len(k)])
  alpha_i <- as.numeric(res$coefficients[k + seq_len(G)])
  names(coefs) <- colnames_X

  vcov_full <- res$vcov_unscaled
  vcov <- vcov_full[seq_len(k), seq_len(k), drop = FALSE]
  dimnames(vcov) <- list(colnames_X, colnames_X)

  # Pad back to the original (pre-screening) row count with NA for dropped
  # rows, so fitted.values stays aligned with the full-length y stored on
  # the panglm object (residuals.panglm() does object$y - fitted.values).
  fitted <- rep(NA_real_, n0)
  fitted[keep_rows] <- as.numeric(exp(X_aug %*% res$coefficients))

  list(coefficients = coefs, vcov = vcov, bread = vcov, fitted.values = fitted,
       loglik = res$loglik, theta = res$theta, dispersion = 1,
       individual_effects = stats::setNames(alpha_i, kept_groups),
       n_used_groups = G, n_dropped_groups = n_dropped,
       # Kept for robust/cluster vcov (see robust_vcov_within_negbin() in
       # vcov.R): the full (k+G)-parameter augmented design/bread, and
       # keep_rows to align a full-length cluster/offset vector to it.
       X_aug = X_aug, keep_rows = keep_rows, vcov_full_aug = vcov_full,
       df.residual = n - k - G - 1, iterations = res$iterations)
}

fit_within_binomial <- function(X, y, group_start, group_size, maxit, tol) {
  k <- ncol(X)

  # BFGS calls the objective and gradient separately at (usually) the same
  # beta; cache the last evaluation so the O(sum T_i^2) DP isn't run twice
  # per iteration.
  cache <- new.env(parent = emptyenv())
  eval_at <- function(beta) {
    if (is.null(cache$beta) || !identical(cache$beta, beta)) {
      cache$beta <- beta
      cache$res <- conditional_logit_loglik_grad_cpp(beta, X, y, group_start, group_size)
    }
    cache$res
  }
  objective <- function(beta) -eval_at(beta)$loglik
  gradient <- function(beta) -eval_at(beta)$gradient

  opt <- stats::optim(rep(0, k), objective, gradient, method = "BFGS",
                       control = list(maxit = maxit, reltol = tol), hessian = TRUE)

  info <- conditional_logit_loglik_grad_cpp(opt$par, X, y, group_start, group_size)
  n_dropped <- length(group_size) - info$n_used_groups
  if (n_dropped > 0) {
    message(n_dropped, " group(s) with no within-group outcome variation dropped ",
            "(all-0 or all-1) -- they carry no information for the conditional likelihood.")
  }

  coefs <- as.numeric(opt$par)
  names(coefs) <- colnames(X)
  vcov <- tryCatch(solve(opt$hessian), error = function(e) matrix(NA_real_, k, k))
  bread <- vcov
  dimnames(vcov) <- dimnames(bread) <- list(colnames(X), colnames(X))

  # df.residual isn't well-defined for conditional logistic regression (the
  # per-group fixed effects are never estimated, only conditioned out), so
  # it's left NA; see n_used_groups/n_dropped_groups instead.
  conditional_fitted <- y - as.numeric(info$score_eta)
  list(coefficients = coefs, vcov = vcov, bread = bread,
       fitted.values = conditional_fitted,
       loglik = -opt$value, dispersion = 1, df.residual = NA_real_,
       n_used_groups = info$n_used_groups, n_dropped_groups = n_dropped,
       iterations = opt$counts[[1]], convergence = opt$convergence)
}

fit_random <- function(X, y, family, group_start, group_size, R, maxit, tol) {
  if (family$family == "gaussian") return(fit_random_gaussian(X, y, group_start, group_size, maxit, tol))
  if (family$family == "poisson") {
    scale <- apply(X, 2, function(col) { s <- stats::sd(col); if (s > 0) s else 1 })
    Xs <- sweep(X, 2, scale, "/")

    # Better starting values than zeros/1: pooled Poisson MLE for beta (on
    # the scaled design, matching the units random_poisson_fit_cpp works
    # in), and a method-of-moments dispersion from the pooled residual
    # over/under-dispersion. Both are cheap and make the Newton path much
    # less likely to wander into a near-singular Hessian region on sparse
    # or skewed panels.
    pooled <- irls_fit_cpp(Xs, y, 1L, 1L, maxit, tol)
    beta_start <- as.numeric(pooled$coefficients)
    mu_pooled <- as.numeric(pooled$fitted.values)
    resid_var <- mean((y - mu_pooled)^2)
    mu_bar <- mean(mu_pooled)
    d_start <- if (resid_var > mu_bar && mu_bar > 0) mu_bar^2 / (resid_var - mu_bar) else 1.0
    if (!is.finite(d_start) || d_start <= 0) d_start <- 1.0

    res <- random_poisson_fit_cpp(Xs, y, group_start, group_size, maxit, tol, beta_start, d_start)
    if (!is.null(res$damping_used) && res$damping_used > 0) {
      warning("random-effects Poisson fit: the Hessian was ill-conditioned and ",
              "required Levenberg-Marquardt damping (max ridge = ",
              format(res$damping_used, digits = 3), "); standard errors may be ",
              "less reliable than usual -- check for sparse or highly skewed panels.",
              call. = FALSE)
    }
    coefs <- as.numeric(res$coefficients) / scale
    names(coefs) <- colnames(X)
    scale_full <- c(scale, 1) # dispersion parameter is unscaled
    vcov_full <- res$vcov_unscaled / outer(scale_full, scale_full)
    vcov <- vcov_full[seq_len(ncol(X)), seq_len(ncol(X)), drop = FALSE]
    dimnames(vcov) <- list(colnames(X), colnames(X))
    fitted <- as.numeric(exp(X %*% coefs))
    return(list(coefficients = coefs, vcov = vcov, fitted.values = fitted,
                loglik = res$loglik, dispersion_param = res$dispersion_param,
                damping_used = res$damping_used,
                df.residual = nrow(X) - ncol(X) - 1, iterations = res$iterations))
  }
  if (family$family == "negbin") {
    scale <- apply(X, 2, function(col) { s <- stats::sd(col); if (s > 0) s else 1 })
    Xs <- sweep(X, 2, scale, "/")
    res <- tryCatch(random_negbin_fit_cpp(Xs, y, group_start, group_size, maxit, tol),
                    error = function(e) NULL)
    if (is.null(res) || !is.finite(res$loglik) || any(!is.finite(res$vcov_unscaled))) {
      # the compiled Newton iteration can fail on strongly overdispersed
      # panels; maximize the same likelihood with the generic optimizer
      return(fit_random_weighted(X, y, family, group_start, group_size, R, maxit, tol,
                                 rep(1, nrow(X)), compiled_start = FALSE))
    }
    coefs <- as.numeric(res$coefficients) / scale
    names(coefs) <- colnames(X)
    scale_full <- c(scale, 1, 1) # a, b are unscaled
    vcov_full <- res$vcov_unscaled / outer(scale_full, scale_full)
    vcov <- vcov_full[seq_len(ncol(X)), seq_len(ncol(X)), drop = FALSE]
    dimnames(vcov) <- list(colnames(X), colnames(X))
    negbin_random_warning(res$a, vcov_full)
    mean_multiplier <- if (res$a > 1) res$b / (res$a - 1) else NA_real_
    fitted <- mean_multiplier * as.numeric(exp(X %*% coefs))
    return(list(coefficients = coefs, vcov = vcov, fitted.values = fitted,
                loglik = res$loglik, negbin_a = res$a, negbin_b = res$b,
                marginal_mean_multiplier = mean_multiplier,
                df.residual = nrow(X) - ncol(X) - 2, iterations = res$iterations))
  }
  if (family$family == "binomial") return(fit_random_binomial(X, y, family, group_start, group_size, R, maxit, tol))
  if (family$family == "gamma") return(fit_random_gamma(X, y, group_start, group_size, R, maxit, tol))
  stop("model = 'random' is not yet implemented for family '", family$family, "'", call. = FALSE)
}

#' Warn about a random-effects negative binomial fit at the boundary
#'
#' In the beta-negative-binomial model of Hausman, Hall and Griliches
#' (1984) the counts given the individual effect are Poisson when
#' `a -> Inf` with `exp(intercept) / a` fixed: the model then reduces to
#' the random-effects Poisson with a gamma effect, the slopes are those of
#' `family = "poisson"` and the intercept and `a` are not separately
#' identified (their joint Hessian is singular). Warns in that case, and
#' whenever the Hessian is not negative definite.
#'
#' @keywords internal
#' @noRd
negbin_random_warning <- function(a, vcov) {
  v <- diag(vcov)
  bad_vcov <- any(!is.finite(vcov)) || any(v <= 0)
  if (is.finite(a) && a > 1e4) {
    warning("random-effects negative binomial: the estimate of 'a' diverges (a = ",
            format(a, digits = 3), "), so the counts show no overdispersion beyond ",
            "the individual effect. The model reduces to the random-effects Poisson ",
            "(family = \"poisson\"), whose slopes it reproduces; the intercept is not ",
            "identified and its standard error is not reliable", call. = FALSE)
  } else if (bad_vcov) {
    warning("random-effects negative binomial: the Hessian is not negative definite, ",
            "so some standard errors are not available", call. = FALSE)
  }
  invisible(NULL)
}

fit_random_gaussian <- function(X, y, group_start, group_size, maxit, tol) {
  has_intercept <- "(Intercept)" %in% colnames(X)
  Xk <- if (has_intercept) X[, setdiff(colnames(X), "(Intercept)"), drop = FALSE] else X
  k <- ncol(Xk); n <- nrow(X); G <- length(group_start)
  Ti <- group_size

  # The within (FE) and between (BW) regressions below only feed the
  # Swamy-Arora variance-component estimate (sigma_v2, sigma_mu2), not the
  # final reported coefficients. They're exactly rank-deficient at common
  # panel shapes: a time-invariant covariate (e.g. baseline sex alongside a
  # pre/post indicator) demeans to a column of zeros in the within
  # regression, and a balanced/deterministic covariate (e.g. a 0/1 time
  # indicator present in every group) has an identical group mean
  # everywhere - collinear with the intercept - in the between regression.
  # stats::lm.fit() handles rank deficiency the same way lm() does (QR with
  # pivoting, aliased coefficients dropped), giving correct residuals
  # without hitting a near-singular arma::solve() in irls_fit_cpp().
  dm <- within_demean_cpp(Xk, y, group_start, group_size)
  fe <- stats::lm.fit(dm$X, dm$y)
  df_fe <- n - G - fe$rank
  sigma_v2 <- sum(fe$residuals^2) / max(1, df_fe)

  gm <- group_means_cpp(X, y, group_start, group_size)
  bw <- stats::lm.fit(gm$Xbar, gm$ybar)
  df_bw <- G - bw$rank
  sigma_1_2 <- sum(bw$residuals^2) / max(1, df_bw)

  sigma_mu2 <- max(0, sigma_1_2 - sigma_v2 / mean(Ti))
  theta <- 1 - sqrt(sigma_v2 / (Ti * sigma_mu2 + sigma_v2))

  qd <- quasi_demean_cpp(X, y, group_start, group_size, theta)
  res <- irls_fit_cpp(qd$X, qd$y, 0L, 0L, maxit, tol)

  coefs <- as.numeric(res$coefficients)
  names(coefs) <- colnames(X)
  df_resid <- n - ncol(X)
  resid_final <- qd$y - qd$X %*% res$coefficients
  sigma2 <- sum(resid_final^2) / max(1, df_resid)
  vcov <- res$vcov_unscaled * sigma2
  dimnames(vcov) <- list(colnames(X), colnames(X))

  ## Closed-form Gaussian log-likelihood of the random-intercept model
  ## y_it = x_it'beta + u_i + e_it, u_i ~ N(0, sigma_mu2), e_it ~ N(0,
  ## sigma_v2), evaluated at the GLS coefficients and the Swamy-Arora
  ## variance-component estimates already computed above (sigma_v2,
  ## sigma_mu2). The within-group covariance is Omega_i = sigma_v2 * I_Ti +
  ## sigma_mu2 * J_Ti (J = all-ones), so det(Omega_i) and the quadratic form
  ## e_i' Omega_i^{-1} e_i both have closed forms via the Sherman-Morrison
  ## identity (standard panel-data result, e.g. Baltagi, *Econometric
  ## Analysis of Panel Data*, ch. 2): with e_i the raw (non-demeaned)
  ## residuals for group i of size Ti,
  ##   log det(Omega_i)   = (Ti - 1) log(sigma_v2) + log(sigma_v2 + Ti * sigma_mu2)
  ##   e_i' Omega_i^{-1} e_i = [sum(e_i^2) - (sigma_mu2 / (sigma_v2 + Ti * sigma_mu2)) * sum(e_i)^2] / sigma_v2
  ## Summed over groups this gives the exact (not REML) log-likelihood at
  ## the fitted beta/variance components; used for AIC/BIC via logLik.panglm().
  resid_raw <- as.numeric(y - X %*% coefs)
  loglik <- if (sigma_v2 > 0) {
    idx <- 1L
    ll <- 0
    for (g in seq_along(group_size)) {
      Ti <- group_size[g]
      ei <- resid_raw[idx:(idx + Ti - 1L)]
      quad <- (sum(ei^2) - (sigma_mu2 / (sigma_v2 + Ti * sigma_mu2)) * sum(ei)^2) / sigma_v2
      logdet <- (Ti - 1) * log(sigma_v2) + log(sigma_v2 + Ti * sigma_mu2)
      ll <- ll - 0.5 * (Ti * log(2 * pi) + logdet + quad)
      idx <- idx + Ti
    }
    ll
  } else {
    NA_real_
  }

  list(coefficients = coefs, vcov = vcov, fitted.values = as.numeric(X %*% res$coefficients),
       loglik = loglik, dispersion = sigma2,
       sigma_v2 = sigma_v2, sigma_mu2 = sigma_mu2, theta = theta,
       df.residual = df_resid, iterations = res$iterations)
}

fit_random_binomial <- function(X, y, family, group_start, group_size, R, maxit, tol) {
  gh <- gauss_hermite_quadrature(R)
  k <- ncol(X)

  start <- irls_fit_cpp(X, y, 2L, family$link_id, maxit, tol)
  theta0 <- c(as.numeric(start$coefficients), log(0.5))

  # BFGS calls the objective and gradient separately at (usually) the same
  # theta; cache the last quadrature evaluation to avoid running it twice.
  cache <- new.env(parent = emptyenv())
  eval_at <- function(theta) {
    if (is.null(cache$theta) || !identical(cache$theta, theta)) {
      cache$theta <- theta
      beta <- theta[seq_len(k)]
      sigma <- exp(theta[k + 1])
      cache$res <- random_binomial_loglik_grad_cpp(beta, sigma, X, y, group_start, group_size,
                                                     gh$nodes, gh$weights, family$link_id)
      cache$sigma <- sigma
    }
    cache
  }
  objective <- function(theta) -eval_at(theta)$res$loglik
  gradient <- function(theta) {
    c <- eval_at(theta)
    g <- c$res$gradient
    g[k + 1] <- g[k + 1] * c$sigma # chain rule: d/d(log sigma) = d/dsigma * sigma
    -g
  }

  opt <- stats::optim(theta0, objective, gradient, method = "BFGS",
                       control = list(maxit = maxit, reltol = tol),
                       hessian = TRUE)

  coefs <- opt$par[seq_len(k)]
  names(coefs) <- colnames(X)
  sigma <- exp(opt$par[k + 1])
  vcov_full <- tryCatch(solve(opt$hessian), error = function(e) matrix(NA_real_, k + 1, k + 1))
  vcov <- vcov_full[seq_len(k), seq_len(k), drop = FALSE]
  dimnames(vcov) <- list(colnames(X), colnames(X))

  eta <- as.numeric(X %*% coefs)
  fitted <- random_binomial_marginal_mean(eta, sigma, family, gh)
  list(coefficients = coefs, vcov = vcov, fitted.values = fitted,
       loglik = -opt$value, sigma = sigma, dispersion = 1,
       quadrature_nodes = R,
       df.residual = nrow(X) - k - 1, iterations = opt$counts[[1]], convergence = opt$convergence)
}

random_binomial_marginal_mean <- function(eta, sigma, family, quadrature = 21L) {
  gh <- if (is.list(quadrature)) quadrature else gauss_hermite_quadrature(quadrature)
  shifts <- sqrt(2) * sigma * gh$nodes
  as.numeric(linkinv_r(outer(eta, shifts, "+"), family$link_id) %*% gh$weights)
}

weighted_random_ok <- function(family) {
  family$family %in% c("ordered", "gaussian", "poisson", "negbin", "binomial", "gamma", "beta", "tobit")
}

#' Build the lagged and initial outcomes of a dynamic panel model
#'
#' @keywords internal
#' @noRd
dynamic_panel_data <- function(formula, data, index, weights, offset) {
  if (length(index) != 2L) {
    stop("dynamic = TRUE needs index = c(<individual>, <time>)", call. = FALSE)
  }
  lhs <- if (length(formula) == 3L) formula[[2L]] else NULL
  if (!is.name(lhs)) {
    stop("dynamic = TRUE needs the response to be a single column of data", call. = FALSE)
  }
  yname <- as.character(lhs)
  if (!yname %in% names(data)) stop("response '", yname, "' not found in data", call. = FALSE)
  time <- data[[index[2]]]
  if (!is.numeric(time)) {
    stop("dynamic = TRUE needs a numeric time index (consecutive waves differ by 1)",
         call. = FALSE)
  }
  id <- data[[index[1]]]
  if (anyDuplicated(data.frame(id, time))) {
    stop("dynamic = TRUE needs one row per individual and time period", call. = FALSE)
  }
  y <- data[[yname]]
  lag_name <- paste0("lag_", yname)
  initial_name <- paste0("initial_", yname)
  o <- order(as.character(id), time)
  ys <- y[o]; ids <- as.character(id)[o]; ts <- time[o]
  n <- length(o)
  previous <- c(NA, seq_len(n - 1L))
  ok <- !is.na(previous) & c(FALSE, ids[-1] == ids[-n]) & c(FALSE, ts[-1] - ts[-n] == 1)
  first <- !duplicated(ids)
  first_row <- match(ids, ids[first])
  if (is.character(ys)) ys <- factor(ys)
  if (is.ordered(ys)) ys <- factor(ys, levels = levels(ys), ordered = FALSE)
  lag <- ys[previous]
  lag[!ok] <- NA
  initial <- ys[first][first_row]
  # undo the sort
  back <- order(o)
  data[[lag_name]] <- lag[back]
  data[[initial_name]] <- initial[back]
  keep <- ok[back] & !first[back]
  if (!any(keep)) stop("no row has an observed previous wave", call. = FALSE)
  term_labels <- attr(stats::terms(formula), "term.labels")
  new_formula <- stats::reformulate(c(term_labels, lag_name, initial_name), response = yname,
                                    intercept = attr(stats::terms(formula), "intercept") == 1L)
  environment(new_formula) <- environment(formula)
  list(formula = new_formula, data = data[keep, , drop = FALSE], full = data, keep = keep,
       weights = if (is.null(weights)) NULL else weights[keep],
       offset = if (is.null(offset)) NULL else offset[keep],
       info = list(response = yname, lag = lag_name, initial = initial_name,
                   dropped = sum(!keep)))
}
