make_ordered_data <- function(G = 300, Tt = 5, seed = 3) {
  set.seed(seed)
  d <- data.frame(id = rep(seq_len(G), each = Tt), t = rep(seq_len(Tt), G))
  d$x1 <- stats::rnorm(nrow(d))
  d$z <- rep(stats::rbinom(G, 1, 0.5), each = Tt)
  u <- rep(stats::rnorm(G), each = Tt)
  ys <- 0.8 * d$x1 - 0.5 * d$z + u + stats::rlogis(nrow(d))
  d$y <- cut(ys, c(-Inf, -1, 0.5, 2, Inf), labels = c("poor", "fair", "good", "vgood"),
             ordered_result = TRUE)
  d$w <- rep(stats::runif(G, 0.5, 2), each = Tt)
  d
}
