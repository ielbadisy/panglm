#' @export
print.panglm <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("panglm(", x$model, ", family = ", x$family$family, ")\n\n", sep = "")
  print(stats::coef(x), digits = digits)
  cat("\nN =", x$nobs, " groups =", x$n_groups, "\n")
  invisible(x)
}
