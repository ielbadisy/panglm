panel_plot_data <- function() {
  set.seed(5)
  d <- expand.grid(id = 1:12, t = 1:6)
  d$first <- rep(c(3, 5, 0), each = 4)[d$id]
  d$y <- rnorm(nrow(d))
  d$y[c(2, 9)] <- NA
  d
}

test_that("panelplot counts observed periods and missing values per individual", {
  d <- panel_plot_data()
  u <- d[-c(1, 14, 15, 30), ]
  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off())
  out <- expect_invisible(panelplot(u, index = c("id", "t")))
  expect_equal(out$id, 1:12)
  expect_equal(out$n_periods, as.vector(table(u$id)))
  expect_equal(sum(out$n_missing), 0)
  ov <- panelplot(u, index = c("id", "t"), value = "y")
  expect_equal(ov$n_missing, as.vector(tapply(is.na(u$y), u$id, sum)))
})

test_that("panelplot sorts individuals by treatment cohort", {
  d <- panel_plot_data()
  d$first[d$id %in% c(2, 7)] <- 5
  d$first[d$id == 11] <- 3
  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off())
  out <- panelplot(d, index = c("id", "t"), value = "y", treat = "first", main = "m",
                   ylab = "Unit", labels = TRUE)
  ufirst <- d$first[match(out$id, d$id)]
  ufirst[ufirst == 0] <- NA
  expect_equal(out$cohort, ufirst)
  expect_false(is.unsorted(out$cohort, na.rm = TRUE))
  expect_true(all(is.na(utils::tail(out$cohort, sum(is.na(out$cohort))))))
  expect_equal(out$id[out$cohort %in% 3], c(1, 3, 4, 11))
})

test_that("panelplot handles discrete values, many units and bad input", {
  d <- panel_plot_data()
  d$g <- ifelse(is.na(d$y), NA, ifelse(d$y > 0, "up", "down"))
  d$b <- d$y > 0
  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off())
  expect_no_error(panelplot(d, c("id", "t"), value = "g"))
  expect_no_error(panelplot(d, c("id", "t"), value = "b", col = c("grey60", "black"),
                            legend = FALSE))
  expect_no_error(panelplot(d, c("id", "t"), value = "y", zlim = c(-1, 1),
                            col = grDevices::gray.colors(10)))
  big <- expand.grid(id = 1:700, t = 1:12)
  big$v <- rnorm(nrow(big))
  expect_equal(nrow(panelplot(big, c("id", "t"), value = "v")), 700)
  expect_error(panelplot(d, c("id", "t"), value = "nope"), "nope")
  expect_error(panelplot(d, "id"), "index")
  expect_error(panelplot(rbind(d, d[1, ]), c("id", "t")), "more than one row")
  expect_error(panelplot(d, c("id", "t"), value = "g", col = "red"), "one colour per level")
  expect_error(panelplot(d, c("id", "t"), value = "y", zlim = c(1, 0)), "zlim")
  expect_error(panelplot(d, c("id", "t"), legend = NA), "legend")
  expect_error(panelplot(d, c("id", "t"), col = c("a", "b")), "one colour")
})
