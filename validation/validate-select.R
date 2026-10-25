# Panel sample selection: panselect() (Wooldridge 1995).
# No R package implements this panel estimator, so it is validated
# component-wise:
#   (a) the per-period probits against stats::glm(binomial("probit"))
#   (b) the inverse Mills ratio against dnorm(xb) / pnorm(xb) by hand
#   (c) the second step against plm::plm(model = "within") and stats::lm()
#       on the hand-built design (lambda_it * d_t, Mundlak means), with the
#       cluster covariance against plm::vcovHC(type = "sss") and
#       sandwich::vcovCL(type = "HC1"); the selection Wald test by hand
#   (d) a single cross-section (one period) against the two-step
#       sampleSelection::heckit() on the Mroz (1987) data
#   (e) a simulation in which selection depends on the idiosyncratic error:
#       the naive within estimator is biased, panselect() recovers the truth
# Data: simulated attrition panel (500 individuals, 5 waves, everyone
# observed in wave 1, response afterwards depends on the health shock),
# and a second panel with intermittent (non-absorbing) non-response.
source("validation/common.R")
suppressPackageStartupMessages({
  library(plm); library(sandwich); library(sampleSelection)
})
se <- function(v) sqrt(diag(v))
out <- list()
add <- function(...) out[[length(out) + 1L]] <<- vcompare(...)

sim_panel <- function(N, TT, rho = 0.7, attrition = TRUE, seed = 1) {
  set.seed(seed)
  d <- data.frame(id = rep(seq_len(N), each = TT), t = rep(seq_len(TT), N))
  xi <- rnorm(N)
  d$x <- rnorm(N * TT) + rep(xi, each = TT)
  d$z <- rnorm(N * TT)
  d$w <- rep(rbinom(N, 1, 0.5), each = TT)          # time invariant
  xbar <- ave(d$x, d$id)
  alpha <- 0.8 * xbar + rep(rnorm(N, sd = 0.3), each = TT)
  v <- rnorm(N * TT)
  d$u <- rho * v + sqrt(1 - rho^2) * rnorm(N * TT)
  d$s <- as.numeric(0.3 + 0.1 * d$t + 0.5 * d$x + d$z + 0.4 * xbar - 0.3 * d$w + v > 0)
  if (attrition) {
    d$s[d$t == 1] <- 1
    # absorbing attrition: once a patient drops out, never observed again
    d$s <- ave(d$s, d$id, FUN = cumprod)
  }
  d$y <- ifelse(d$s == 1, 1 + 0.5 * d$x + 0.3 * d$w + alpha + d$u, NA)
  d
}

# hand implementation of both steps
hand <- function(d, model, mundlak, lambda) {
  d$mean_x <- ave(d$x, d$id); d$mean_z <- ave(d$z, d$id)
  rhs <- if (mundlak) "x + z + w + mean_x + mean_z" else "x + z + w"
  d$imr <- 0
  probits <- list()
  for (p in sort(unique(d$t))) {
    r <- d$t == p
    if (all(d$s[r] == 1) || all(d$s[r] == 0)) next
    g <- glm(as.formula(paste("s ~", rhs)), data = d[r, ], family = binomial(link = "probit"),
             control = glm.control(epsilon = 1e-12, maxit = 100))
    xb <- predict(g, type = "link")
    d$imr[r] <- dnorm(xb) / pnorm(xb)
    probits[[as.character(p)]] <- coef(g)
  }
  pp <- as.numeric(names(probits))
  # syntactic names lam_<p> here, renamed lambda:<p> as in panselect()
  lam <- if (lambda == "common") "lambda" else paste0("lam_", pp)
  if (lambda == "common") d$lambda <- d$imr
  for (p in pp) d[[paste0("lam_", p)]] <- d$imr * (d$t == p)
  ds <- d[d$s == 1, ]
  rhs2 <- paste(c("x", "w", if (model == "pooling" && mundlak) c("mean_x", "mean_z"), lam),
                collapse = " + ")
  f <- as.formula(paste("y ~", rhs2))
  fit <- if (model == "within") plm::plm(f, ds, index = c("id", "t"), model = "within") else lm(f, ds)
  list(d = d, probits = probits, fit = fit, lam = sub("^lam_", "lambda:", lam))
}

for (dname in c("attrition", "intermittent")) {
  d <- sim_panel(500, 5, attrition = dname == "attrition", seed = 11)
  for (model in c("within", "pooling")) for (mundlak in c(TRUE, FALSE)) for (lam in c("period", "common")) {
    if (model == "pooling" && !mundlak && lam == "common") next
    lab <- sprintf("panselect %s%s, lambda %s", model, if (mundlak) " + Mundlak" else "", lam)
    a <- panselect(s ~ x + z + w, y ~ x + w, data = d, index = c("id", "t"),
                   model = model, mundlak = mundlak, lambda = lam, boot = 0)
    h <- hand(d, model, mundlak, lam)
    # (a) probits
    pa <- unlist(a$probit[names(h$probits)]); ph <- unlist(h$probits)
    add(lab, dname, "probit coefficients", pa, list(stats = ph), 1e-7)
    # (b) inverse Mills ratio by hand
    add(lab, dname, "inverse Mills ratio", a$imr, list(hand = h$d$imr), 1e-7)
    # (c) second step
    cf <- coef(h$fit); names(cf) <- sub("^lam_", "lambda:", names(cf))
    refname <- if (model == "within") "plm" else "lm/sandwich"
    add(lab, dname, "coefficients", coef(a),
        stats::setNames(list(cf[names(coef(a))]), refname), 1e-8)
    vc_cl <- if (model == "within") {
      plm::vcovHC(h$fit, method = "arellano", type = "sss", cluster = "group")
    } else {
      sandwich::vcovCL(h$fit, cluster = ~id, type = "HC1")
    }
    v_cl <- se(vc_cl); names(v_cl) <- sub("^lam_", "lambda:", names(v_cl))
    v_cs <- se(vcov(h$fit)); names(v_cs) <- sub("^lam_", "lambda:", names(v_cs))
    add(lab, dname, "classical SE", se(vcov(a, type = "classical")),
        stats::setNames(list(v_cs[names(coef(a))]), refname), 1e-8)
    add(lab, dname, "cluster SE", se(vcov(a, type = "cluster")),
        stats::setNames(list(v_cl[names(coef(a))]), refname), 1e-8)
    # selection test by hand: Wald with the cluster covariance
    dimnames(vc_cl) <- lapply(dimnames(vc_cl), function(x) sub("^lam_", "lambda:", x))
    b <- cf[h$lam]
    W <- as.numeric(t(b) %*% solve(vc_cl[h$lam, h$lam, drop = FALSE], b))
    add(lab, dname, "selection Wald statistic", a$selection_test$cluster$statistic,
        list(hand = W), 1e-7)
  }
}

# (d) one period: two-step Heckman (heckit) on the Mroz (1987) data
data("Mroz87", package = "sampleSelection")
Mroz87$kids <- as.numeric(Mroz87$kids5 + Mroz87$kids618 > 0)
Mroz87$id <- seq_len(nrow(Mroz87)); Mroz87$t <- 1
hk <- sampleSelection::heckit(lfp ~ age + I(age^2) + faminc + kids + educ,
                              wage ~ exper + I(exper^2) + educ + city, data = Mroz87)
for (lam in c("common", "period")) {
  a <- panselect(lfp ~ age + I(age^2) + faminc + kids + educ,
                 wage ~ exper + I(exper^2) + educ + city, data = Mroz87,
                 index = c("id", "t"), model = "pooling", lambda = lam, boot = 0)
  lab <- sprintf("panselect pooling, one period, lambda %s", lam)
  hc <- coef(hk)
  add(lab, "Mroz87", "probit coefficients", a$probit[["1"]],
      list(sampleSelection = hc[1:6]), 1e-6)
  add(lab, "Mroz87", "outcome coefficients", coef(a),
      list(sampleSelection = hc[7:12]), 1e-6)
  add(lab, "Mroz87", "inverse Mills ratio", a$imr[Mroz87$lfp == 1],
      list(sampleSelection = hk$invMillsRatio[Mroz87$lfp == 1]), 1e-6)
}

# (e) simulation: selection on the idiosyncratic error biases the naive
# within estimator; panselect recovers the true slope 0.5
R <- 40
est <- t(vapply(seq_len(R), function(r) {
  d <- sim_panel(1000, 5, rho = 0.8, attrition = FALSE, seed = 100 + r)
  naive <- coef(plm::plm(y ~ x, d[d$s == 1, ], index = c("id", "t"), model = "within"))[["x"]]
  ps <- coef(panselect(s ~ x + z + w, y ~ x, data = d, index = c("id", "t"),
                       model = "within", boot = 0))[["x"]]
  c(naive = naive, panselect = ps)
}, numeric(2)))
m <- colMeans(est)
cat(sprintf("simulation (R = %d): mean naive FE = %.4f, mean panselect = %.4f, truth 0.5\n",
            R, m[["naive"]], m[["panselect"]]))
if (abs(m[["naive"]] - 0.5) < 0.03) stop("the naive within estimator is not biased in the simulation")
add("panselect within (simulation, R = 40)", "selection on u_it", "mean slope",
    c(x = m[["panselect"]]), list(truth = 0.5), 0.01)

# bootstrap: reproducible with a seed, and the global RNG state is untouched
d <- sim_panel(300, 4, seed = 5)
set.seed(9); before <- runif(1); set.seed(9)
b1 <- panselect(s ~ x + z + w, y ~ x + w, data = d, index = c("id", "t"), boot = 30, seed = 3)
after <- runif(1)
b2 <- panselect(s ~ x + z + w, y ~ x + w, data = d, index = c("id", "t"), boot = 30, seed = 3)
if (!identical(before, after)) stop("panselect() changed the global random number state")
add("panselect within, bootstrap", "attrition", "bootstrap SE (same seed twice)",
    se(vcov(b1)), list(panglm = se(vcov(b2))), 1e-12)

vsave(do.call(rbind, out), name = "select")
