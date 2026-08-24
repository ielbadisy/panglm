arguments <- commandArgs(trailingOnly = FALSE)
file_argument <- arguments[grepl("^--file=", arguments)]
script_path <- if (length(file_argument)) {
  normalizePath(sub("^--file=", "", file_argument[1]))
} else {
  normalizePath("paper/replication.R")
}
paper_directory <- dirname(script_path)
old_directory <- setwd(paper_directory)
on.exit(setwd(old_directory), add = TRUE)

options(width = 70, useFancyQuotes = FALSE)

library("panglm")
library("plm")

data("Grunfeld", package = "plm")

pool_fit <- panglm(
  inv ~ value + capital, data = Grunfeld,
  index = c("firm", "year"), model = "pooling", family = "gaussian"
)
within_fit <- panglm(
  inv ~ value + capital, data = Grunfeld,
  index = c("firm", "year"), model = "within", family = "gaussian"
)
random_fit <- panglm(
  inv ~ value + capital, data = Grunfeld,
  index = c("firm", "year"), model = "random", family = "gaussian"
)

plm_fit <- plm(
  inv ~ value + capital, data = Grunfeld,
  index = c("firm", "year"), model = "within"
)
stopifnot(max(abs(coef(within_fit) - coef(plm_fit))) < 1e-6)

cluster_covariance <- vcov(pool_fit, type = "cluster")
stopifnot(all(is.finite(cluster_covariance)))

set.seed(42)
Grunfeld$count <- rpois(
  nrow(Grunfeld), lambda = exp(0.5 + 0.0002 * Grunfeld$value)
)
poisson_fe <- panglm(
  count ~ value + capital, data = Grunfeld,
  index = c("firm", "year"), model = "within", family = "poisson"
)
stopifnot(all(is.finite(coef(poisson_fe))))

pdf("replication-coefficients.pdf", width = 6.5, height = 4.6)
plot(within_fit, which = "coefficients")
dev.off()

sessionInfo()
