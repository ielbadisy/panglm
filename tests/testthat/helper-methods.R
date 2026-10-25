# Small example fits of every model class, shared by the tests of the
# generic methods (tidy, glance, confint, predict, plot). Built once on
# first use.
.methods_cache <- new.env(parent = emptyenv())

methods_data <- function(name) {
  if (!is.null(.methods_cache[[paste0("d_", name)]])) return(.methods_cache[[paste0("d_", name)]])
  d <- switch(name,
    copd = {
      data("copd", package = "panglm", envir = environment())
      set.seed(1)
      copd$cost <- ifelse(copd$exacerbations > 0,
                          exp(6 + 0.1 * copd$crp + rnorm(nrow(copd), 0, 0.5)), 0)
      copd
    },
    gmm = {
      set.seed(3)
      N <- 100; TT <- 6
      d <- data.frame(id = rep(seq_len(N), each = TT), t = rep(seq_len(TT), N))
      mu <- rep(rnorm(N), each = TT)
      d$x <- rnorm(N * TT)
      d$y <- stats::ave(mu + 0.3 * d$x + rnorm(N * TT), d$id,
                        FUN = function(e) stats::filter(e, 0.4, method = "recursive"))
      d
    },
    iv = {
      set.seed(2)
      d <- data.frame(id = rep(1:50, each = 4), t = rep(1:4, 50))
      d$z <- rnorm(200); d$z2 <- rnorm(200); d$u <- rnorm(200)
      d$w <- rnorm(200) + rep(rnorm(50), each = 4)
      d$x <- d$z + 0.5 * d$z2 + d$u + rnorm(200)
      d$y <- 1 + d$x + 0.5 * d$w + d$u + rep(rnorm(50), each = 4) + rnorm(200)
      d
    },
    qr = {
      set.seed(2)
      d <- data.frame(id = rep(1:100, each = 5), t = rep(1:5, 100))
      d$x <- rnorm(500) + rep(rnorm(100), each = 5)
      d$y <- d$x + rep(rnorm(100), each = 5) + (1 + 0.2 * abs(d$x)) * rnorm(500)
      d
    },
    sel = {
      set.seed(3)
      d <- data.frame(id = rep(1:150, each = 3), t = rep(1:3, 150))
      d$x <- rnorm(450); d$z <- rnorm(450); v <- rnorm(450)
      d$s <- as.numeric(0.5 + d$x + d$z + v > 0)
      d$y <- ifelse(d$s == 1, 1 + d$x + 0.5 * v + rnorm(450), NA)
      d
    },
    fmm = {
      set.seed(3)
      d <- data.frame(id = rep(1:150, each = 3), t = rep(1:3, 150), x = rnorm(450))
      d$y <- rpois(450, exp(ifelse(rep(rbinom(150, 1, 0.4), each = 3) == 1, 1.5, 0) + 0.3 * d$x))
      d
    },
    did = {
      set.seed(2)
      d <- expand.grid(id = 1:150, t = 1:8)
      d$first <- rep(c(4, 6, NA), each = 50)[d$id]
      d$y <- rnorm(150)[d$id] + 0.2 * d$t +
        ifelse(!is.na(d$first) & d$t >= d$first, 1, 0) + rnorm(nrow(d))
      d
    })
  assign(paste0("d_", name), d, envir = .methods_cache)
  d
}

methods_fit <- function(name) {
  if (!is.null(.methods_cache[[name]])) return(.methods_cache[[name]])
  copd <- methods_data("copd")
  fit <- switch(name,
    twopart = twopart(cost ~ treatment + age + crp, data = copd, index = c("id", "visit"),
                      positive = "gamma"),
    pangmm = pangmm(y ~ lag(y, 1) + x | lag(y, 2:99), data = methods_data("gmm"),
                    index = c("id", "t"), model = "twostep"),
    paniv = paniv(y ~ x + w | z + z2 + w, data = methods_data("iv"), index = c("id", "t"),
                  model = "within"),
    panqr = panqr(y ~ x, data = methods_data("qr"), index = c("id", "t"),
                  tau = c(0.25, 0.5, 0.9), model = "mmqr"),
    panselect = panselect(s ~ x + z, y ~ x, data = methods_data("sel"), index = c("id", "t"),
                          boot = 0),
    panfmm = panfmm(y ~ x, data = methods_data("fmm"), index = c("id", "t"), nstart = 2,
                    seed = 1),
    attgt = attgt(y ~ 1, methods_data("did"), index = c("id", "t"), treat = "first"),
    etwfe = etwfe(y ~ 1, methods_data("did"), index = c("id", "t"), treat = "first",
                  control = "never"),
    eventstudy = eventstudy(y ~ 1, data = methods_data("did"), index = c("id", "t"),
                            treat = "first", leads = 3, lags = 2))
  assign(name, fit, envir = .methods_cache)
  fit
}
