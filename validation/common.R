# Shared helpers for the panglm validation scripts.
#
# Each script validation/validate-<model>.R fits one model on one dataset
# with panglm and with one or more established packages, then compares the
# estimates (coefficients, standard errors, log-likelihood, variance
# components, test statistics). Presentation outputs (plots, printed
# tables) are not compared: they have no counterpart in other packages.
#
# Run a script from the package root: Rscript validation/validate-<model>.R
# Run all of them and write validation/REPORT.md: Rscript validation/run-all.R

suppressPackageStartupMessages(library(panglm))

# Compare one quantity between panglm and each reference package.
# panglm: named numeric vector; refs: named list of numeric vectors in the
# same order (list names are package names). A value passes when
# |panglm - reference| <= tol * max(1, |reference|).
vcompare <- function(model, data, quantity, panglm, refs, tol) {
  p <- as.numeric(panglm)
  terms <- names(panglm)
  if (is.null(terms)) terms <- if (length(p) == 1L) quantity else as.character(seq_along(p))
  rows <- lapply(names(refs), function(pkg) {
    ref <- as.numeric(refs[[pkg]])
    if (length(ref) != length(p)) {
      stop(sprintf("%s / %s: panglm has %d values, %s has %d",
                   model, quantity, length(p), pkg, length(ref)), call. = FALSE)
    }
    d <- abs(p - ref)
    data.frame(model = model, data = data, quantity = quantity, term = terms,
               package = pkg, panglm = p, reference = ref, absdiff = d, tol = tol,
               pass = d <= tol * pmax(1, abs(ref)), stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

# Bind comparison tables, save validation/results/<name>.csv, print a short
# summary, and stop with an error if any comparison fails.
vsave <- function(..., name) {
  tab <- do.call(rbind, list(...))
  dir.create("validation/results", showWarnings = FALSE, recursive = TRUE)
  utils::write.csv(tab, file.path("validation/results", paste0(name, ".csv")), row.names = FALSE)
  agg <- stats::aggregate(cbind(absdiff, fail = !pass) ~ model + data + quantity + package,
                          data = tab, FUN = max)
  agg$fail <- agg$fail > 0
  print(agg, row.names = FALSE)
  if (any(!tab$pass)) {
    print(tab[!tab$pass, ], row.names = FALSE)
    stop(sprintf("validation '%s' failed for %d comparison(s)", name, sum(!tab$pass)), call. = FALSE)
  }
  invisible(tab)
}
