make_cost_data <- function() {
  data(copd, package = "panglm", envir = environment())
  set.seed(11)
  u <- stats::rnorm(length(unique(copd$id)), sd = 0.4)
  copd$cost <- exp(1 + 0.3 * copd$treatment + 0.01 * copd$age + 0.2 * copd$crp +
                     u[match(copd$id, sort(unique(copd$id)))]) *
    stats::rgamma(nrow(copd), 2, 2)
  copd$cost0 <- ifelse(stats::runif(nrow(copd)) < stats::plogis(0.5 + 0.3 * copd$crp),
                       copd$cost, 0)
  copd
}
