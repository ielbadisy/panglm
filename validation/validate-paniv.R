# Panel instrumental variables: paniv() with model = "within", "pooling"
# and "random" against plm::plm() (same estimators, same inst.method),
# AER::ivreg() (pooled 2SLS) and fixest::feols() (FE-2SLS).
# Data: plm::Crime (90 North Carolina counties, 1981-1987), with the
# police per capita (lpolpc) and, in the just-identified model, the
# probability of arrest (lprbarr) treated as endogenous, instrumented by
# tax revenue per capita (ltaxpc) and the offence mix (lmix). An
# unbalanced version drops 60 random county-years.
source("validation/common.R")
suppressPackageStartupMessages(library(plm))
data("Crime", package = "plm")
set.seed(20261025)
CrimeU <- Crime[-sample(nrow(Crime), 60), ]
idx <- c("county", "year")
exog <- "lprbconv + lprbpris + lavgsen + ldensity + lwcon + lpctymle"
f_just <- stats::as.formula(paste("lcrmrte ~ lprbarr + lpolpc +", exog, "| ltaxpc + lmix +", exog))
f_over <- stats::as.formula(paste("lcrmrte ~ lpolpc + lprbarr +", exog,
                                  "| ltaxpc + lmix + lprbarr +", exog))
se <- function(v) sqrt(diag(v))
out <- list()
add <- function(...) out[[length(out) + 1L]] <<- vcompare(...)

for (dname in c("Crime", "Crime unbalanced")) {
  d <- if (dname == "Crime") Crime else CrimeU
  for (spec in c("just", "over")) {
    f <- if (spec == "just") f_just else f_over
    lab <- function(m) sprintf("paniv %s (%s-identified)", m, spec)

    # FE-2SLS (one-way and two-way) against plm
    for (eff in c("individual", "twoways")) {
      a <- paniv(f, d, index = idx, model = "within", effect = eff)
      b <- plm::plm(f, d, index = idx, model = "within", effect = eff)
      m <- lab(if (eff == "twoways") "within twoways" else "within")
      add(m, dname, "coefficients", coef(a), list(plm = coef(b)[names(coef(a))]), 1e-8)
      add(m, dname, "classical SE", se(vcov(a)), list(plm = se(vcov(b))[names(coef(a))]), 1e-8)
      vb <- plm::vcovHC(b, method = "arellano", type = "sss", cluster = "group")
      add(m, dname, "cluster SE", se(vcov(a, type = "cluster")),
          list(plm = se(vb)[names(coef(a))]), 1e-8)
    }
    # fixest FE-2SLS: same estimator; fixest's default "nonnested" rule keeps
    # one fixed-effect parameter in K, so K.fixef = "none" is set to obtain the
    # factor G/(G-1) (n-1)/(n-k) used by panglm and plm::vcovHC(type = "sss")
    a <- paniv(f, d, index = idx, model = "within", vcov = "cluster")
    endo <- a$endogenous
    rhs_exog <- gsub(" ", "", strsplit(exog, "\\+")[[1]])
    excl <- if (spec == "just") c("ltaxpc", "lmix") else c("ltaxpc", "lmix")
    exog_all <- if (spec == "just") rhs_exog else c("lprbarr", rhs_exog)
    ff <- stats::as.formula(paste("lcrmrte ~", paste(exog_all, collapse = "+"), "| county |",
                                  paste(endo, collapse = "+"), "~", paste(excl, collapse = "+")))
    fx <- fixest::feols(ff, data = d, cluster = ~county, notes = FALSE,
                        ssc = fixest::ssc(K.fixef = "none"))
    cf <- coef(fx); names(cf) <- sub("^fit_", "", names(cf))
    sf <- se(vcov(fx)); names(sf) <- names(cf)
    add(lab("within"), dname, "coefficients", coef(a), list(fixest = cf[names(coef(a))]), 1e-8)
    add(lab("within"), dname, "cluster SE", se(vcov(a)), list(fixest = sf[names(coef(a))]), 1e-8)
    if (spec == "over") {
      # Sargan from 2SLS residuals with n = number of observations: fixest's
      # fitstat 'sargan' uses the same definition
      a0 <- paniv(f, d, index = idx, model = "within")
      fx0 <- fixest::feols(ff, data = d, notes = FALSE)
      add(lab("within"), dname, "Sargan statistic", a0$overid$statistic,
          list(fixest = fixest::fitstat(fx0, "sargan")$sargan$stat), 1e-8)
    }

    # pooled 2SLS against plm and AER::ivreg
    a <- paniv(f, d, index = idx, model = "pooling")
    b <- plm::plm(f, d, index = idx, model = "pooling")
    iv <- AER::ivreg(f, data = d)
    add(lab("pooling"), dname, "coefficients", coef(a),
        list(plm = coef(b), AER = coef(iv)[names(coef(a))]), 1e-8)
    add(lab("pooling"), dname, "classical SE", se(vcov(a)),
        list(plm = se(vcov(b)), AER = se(vcov(iv))[names(coef(a))]), 1e-8)
    add(lab("pooling"), dname, "cluster SE", se(vcov(a, type = "cluster")),
        list(plm = se(plm::vcovHC(b, method = "arellano", type = "sss", cluster = "group"))), 1e-8)
    diag_iv <- summary(iv, diagnostics = TRUE)$diagnostics
    weak <- diag_iv[grepl("^Weak", rownames(diag_iv)), "statistic"]
    add(lab("pooling"), dname, "first-stage F", a$first_stage$F, list(AER = weak), 1e-8)
    if (spec == "over") {
      add(lab("pooling"), dname, "Sargan statistic", a$overid$statistic,
          list(AER = diag_iv["Sargan", "statistic"]), 1e-8)
    }

    # random effects G2SLS and EC2SLS against plm (Swamy-Arora components)
    for (im in c("bvk", "baltagi")) {
      a <- paniv(f, d, index = idx, model = "random", inst.method = im)
      b <- plm::plm(f, d, index = idx, model = "random", inst.method = im)
      m <- lab(paste("random", im))
      add(m, dname, "coefficients", coef(a), list(plm = coef(b)), 1e-8)
      add(m, dname, "classical SE", se(vcov(a)), list(plm = se(vcov(b))), 1e-8)
      add(m, dname, "variance components", a$ercomp$sigma2,
          list(plm = b$ercomp$sigma2[c("idios", "id")]), 1e-8)
      # plm stores theta once for balanced panels and per observation otherwise
      th <- as.numeric(b$ercomp$theta)
      th_a <- a$ercomp$theta[match(as.character(d$county[order(d$county, d$year)]),
                                   as.character(sort(unique(d$county))))]
      add(m, dname, "theta", if (length(th) == 1L) th_a[1] else th_a,
          list(plm = th), 1e-8)
      if (im == "bvk") {
        # plm's vcovHC uses the quasi-demeaned instruments, i.e. the G2SLS
        # first stage; for EC2SLS it reuses them (and the residuals scaled by
        # 1 / sigma_nu), so the EC2SLS cluster covariance is not compared
        add(m, dname, "cluster SE", se(vcov(a, type = "cluster")),
            list(plm = se(plm::vcovHC(b, method = "arellano", type = "sss", cluster = "group"))),
            1e-8)
      }
    }
  }
}
do.call(vsave, c(out, name = "paniv"))
