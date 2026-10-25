local_threads_env <- function(env = parent.frame()) {
  vars <- c("OMP_THREAD_LIMIT", "_R_CHECK_LIMIT_CORES_")
  old <- as.list(Sys.getenv(vars, unset = NA))
  oldopt <- options(panglm.threads = NULL)
  Sys.unsetenv(vars)
  restore <- function() {
    options(oldopt)
    for (v in names(old)) {
      if (is.na(old[[v]])) Sys.unsetenv(v)
      else do.call(Sys.setenv, old[v])
    }
  }
  do.call(on.exit, list(bquote(.(restore)()), add = TRUE), envir = env)
  invisible()
}

test_that("options(panglm.threads) and the check env vars cap the threads", {
  local_threads_env()
  nt <- function() panglm:::panglm_num_threads_cpp()
  expect_identical(nt(), 0L)
  options(panglm.threads = 3)
  expect_identical(nt(), 3L)
  Sys.setenv(OMP_THREAD_LIMIT = "2")
  expect_identical(nt(), 2L)
  Sys.unsetenv("OMP_THREAD_LIMIT")
  Sys.setenv(`_R_CHECK_LIMIT_CORES_` = "TRUE")
  expect_identical(nt(), 2L)
  options(panglm.threads = 1)
  expect_identical(nt(), 1L)
  Sys.setenv(`_R_CHECK_LIMIT_CORES_` = "false")
  options(panglm.threads = NULL)
  expect_identical(nt(), 0L)
})

test_that("within and random estimates do not depend on the number of threads", {
  local_threads_env()
  data(copd, package = "panglm")
  f <- function(k, model, family = "gaussian") {
    options(panglm.threads = k)
    coef(panglm(exacerbations ~ crp + fev1, copd, c("id", "visit"),
                model = model, family = family))
  }
  for (m in c("within", "random")) {
    expect_equal(f(1, m), f(4, m), tolerance = 1e-10)
  }
  expect_equal(f(1, "within", "poisson"), f(4, "within", "poisson"),
               tolerance = 1e-10)
})
