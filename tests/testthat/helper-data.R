skip_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) testthat::skip(paste(pkg, "not installed"))
}

# Built at top-level (outside any test_that()/local() block) so that pglm's
# match.call()/parent.frame()-based NSE can resolve the augmented columns
# reliably when the comparison tests below call pglm() from inside a
# test_that() block.
.panglm_test_data <- local({
  out <- list()
  if (requireNamespace("plm", quietly = TRUE)) {
    data("Grunfeld", package = "plm", envir = environment())
    set.seed(1)
    Grunfeld$count <- rpois(nrow(Grunfeld), lambda = exp(0.5 + 0.01 * Grunfeld$value / 100))
    out$Grunfeld_poisson <- Grunfeld

    data("Grunfeld", package = "plm", envir = environment())
    set.seed(42)
    Grunfeld$u <- rep(rgamma(10, shape = 5, rate = 5), each = 20)
    Grunfeld$count <- rpois(nrow(Grunfeld), lambda = Grunfeld$u * exp(1 + 0.0002 * Grunfeld$value))
    out$Grunfeld_negbin <- Grunfeld
  }
  out
})
