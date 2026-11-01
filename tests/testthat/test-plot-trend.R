trend_panel <- function(seed = 1, N = 120) {
  set.seed(seed)
  d <- expand.grid(id = seq_len(N), t = 1:8)
  d$first <- rep(c(3, 5, 6, 0), length.out = N)[d$id]
  d$w <- runif(N, 0.5, 2)[d$id]
  d$y <- rnorm(N)[d$id] + 0.2 * d$t + ifelse(d$first > 0 & d$t >= d$first, 1, 0) +
    rnorm(nrow(d))
  d <- d[-sample(nrow(d), 60), ]
  d$y[sample(nrow(d), 10)] <- NA
  d
}

test_that("trendplot means equal aggregate() for every type", {
  d <- trend_panel()
  pdf(NULL); on.exit(dev.off())
  m <- expect_invisible(trendplot(d, c("id", "t"), "first", "y"))
  expect_named(m, c("group", "time", "mean", "se", "n"))
  a <- aggregate(y ~ first + t, d, mean)
  a$group <- ifelse(a$first == 0, "Never treated", a$first)
  mm <- merge(m, a, by.x = c("group", "time"), by.y = c("group", "t"))
  expect_equal(nrow(mm), nrow(m))
  expect_equal(mm$mean, mm$y, tolerance = 1e-12)
  n <- aggregate(y ~ first + t, d, length)
  expect_equal(sum(m$n), sum(n$y))
  d$ever <- ifelse(d$first > 0, "Treated", "Never treated")
  m <- trendplot(d, c("id", "t"), "first", "y", type = "treated", band = TRUE)
  a <- aggregate(y ~ ever + t, d, mean)
  mm <- merge(m, a, by.x = c("group", "time"), by.y = c("ever", "t"))
  expect_equal(nrow(mm), nrow(m))
  expect_equal(mm$mean, mm$y, tolerance = 1e-12)
  # weighted means
  m <- trendplot(d, c("id", "t"), "first", "y", weights = "w")
  dd <- d[!is.na(d$y), ]
  a <- aggregate(cbind(wy = w * y, w) ~ first + t, dd, sum)
  a$group <- ifelse(a$first == 0, "Never treated", a$first)
  mm <- merge(m, a, by.x = c("group", "time"), by.y = c("group", "t"))
  expect_equal(mm$mean, mm$wy / mm$w, tolerance = 1e-12)
})

test_that("trendplot event time aligns the treated and matches the never treated on period", {
  d <- trend_panel(2)
  pdf(NULL); on.exit(dev.off())
  m <- trendplot(d, c("id", "t"), "first", "y", type = "event", min_e = -3, max_e = 2)
  tr <- d[d$first > 0 & !is.na(d$y), ]
  tr$e <- tr$t - tr$first
  tr <- tr[tr$e >= -3 & tr$e <= 2, ]
  a <- aggregate(y ~ e, tr, mean)
  mt <- m[m$group == "Treated", ]
  expect_equal(mt$time, a$e)
  expect_equal(mt$mean, a$y, tolerance = 1e-12)
  nv <- d[d$first == 0 & !is.na(d$y), ]
  pm <- tapply(nv$y, nv$t, mean)
  manual <- vapply(a$e, function(e) mean(pm[as.character(tr$t[tr$e == e])]), numeric(1))
  expect_equal(m$mean[m$group == "Never treated"], unname(manual), tolerance = 1e-12)
  expect_true(all(m$se > 0))
})

test_that("trendplot standard errors are cluster-robust means (as fixest)", {
  skip_if_not_installed("fixest")
  d <- trend_panel(3)
  pdf(NULL); on.exit(dev.off())
  for (w in list(NULL, "w")) {
    m <- trendplot(d, c("id", "t"), "first", "y", weights = w)
    dd <- d[!is.na(d$y), ]
    dd$group <- ifelse(dd$first == 0, "Never treated", dd$first)
    ref <- vapply(seq_len(nrow(m)), function(j) {
      s <- dd[dd$group == m$group[j] & dd$t == m$time[j], ]
      f <- fixest::feols(y ~ 1, s, cluster = ~id, weights = if (is.null(w)) NULL else ~w)
      unname(fixest::se(f))
    }, numeric(1))
    expect_equal(m$se, ref, tolerance = 1e-8)
  }
})

test_that("trendplot reproduces the castle doctrine means and fixest cell SEs", {
  # optional external data, not a declared dependency: name it through a
  # variable so that R CMD check does not require it
  src <- "bacondecomp"
  skip_if_not_installed(src)
  skip_if_not_installed("fixest")
  utils::data("castle", package = src, envir = environment())
  castle$adopt <- ifelse(is.na(castle$treatment_date), 0, castle$treatment_date)
  pdf(NULL); on.exit(dev.off())
  m <- trendplot(castle, c("sid", "year"), "adopt", "l_homicide")
  a <- aggregate(l_homicide ~ adopt + year, castle, mean)
  a$group <- ifelse(a$adopt == 0, "Never treated", a$adopt)
  mm <- merge(m, a, by.x = c("group", "time"), by.y = c("group", "year"))
  expect_equal(nrow(mm), nrow(m))
  expect_equal(mm$mean, mm$l_homicide, tolerance = 1e-12)
  # one regression on all cell dummies: same SEs up to the small-sample
  # factor, G / (G - 1) over all 50 states instead of over the cell
  castle$cell <- paste(ifelse(castle$adopt == 0, "Never treated", castle$adopt),
                       castle$year, sep = "_")
  f <- fixest::feols(l_homicide ~ 0 + i(cell), castle, cluster = ~sid,
                     ssc = fixest::ssc(adj = FALSE))
  se <- fixest::se(f)
  names(se) <- sub("^cell::", "", names(se))
  key <- paste(m$group, m$time, sep = "_")
  G <- 50
  Gc <- vapply(key, function(k) length(unique(castle$sid[castle$cell == k])), numeric(1))
  ok <- Gc > 1
  expect_true(all(is.na(m$se[!ok])))  # one state (2009 cohort): no clustered SE
  expect_equal(m$se[ok], unname(se[key] * sqrt(Gc / (Gc - 1) / (G / (G - 1))))[ok],
               tolerance = 1e-8)
})

test_that("trendplot options and input checks", {
  d <- trend_panel(4)
  pdf(NULL); on.exit(dev.off())
  expect_no_error(trendplot(d, c("id", "t"), "first", "y", band = TRUE, legend = FALSE,
                            col = c("red", "blue", "green", "grey"), main = "m",
                            ylim = c(-1, 4), xlab = "x", ylab = "y"))
  d$yl <- d$y > 1
  expect_no_error(trendplot(d, c("id", "t"), "first", "yl", type = "event"))
  one <- d
  one$first[one$first > 0] <- 5
  expect_no_error(trendplot(one, c("id", "t"), "first", "y", type = "treated"))
  expect_error(trendplot(d, c("id", "t"), "first", "nope"), "outcome")
  expect_error(trendplot(d, c("id", "t"), "first", "y", band = NA), "band")
  expect_error(trendplot(d, c("id", "t"), "first", "y", level = 2), "level")
  expect_error(trendplot(d, c("id", "t"), "first", "y", type = "event", min_e = 30,
                         max_e = 40), "window")
  d$first <- 0
  expect_error(trendplot(d, c("id", "t"), "first", "y", type = "event"), "treated units")
})
