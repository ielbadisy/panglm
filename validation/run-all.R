# Run every validation/validate-*.R script and write validation/REPORT.md.
# Run from the package root: Rscript validation/run-all.R
scripts <- sort(list.files("validation", pattern = "^validate-.*\\.R$", full.names = TRUE))
status <- vapply(scripts, function(s) {
  message("== ", s)
  res <- system2("Rscript", s, stdout = "", stderr = "")
  if (res == 0) "pass" else "FAIL"
}, character(1))
files <- file.path("validation/results", sub("^validate-(.*)\\.R$", "\\1.csv", basename(scripts)))
tabs <- lapply(files[file.exists(files)], utils::read.csv, stringsAsFactors = FALSE)
tab <- do.call(rbind, tabs)
agg <- stats::aggregate(cbind(absdiff, n = 1, fail = !pass) ~ model + data + quantity + package,
                        data = tab, FUN = function(x) if (length(x)) max(x) else NA)
n <- stats::aggregate(pass ~ model + data + quantity + package, data = tab, FUN = length)
agg$n <- n$pass
fmt <- function(x) formatC(x, format = "e", digits = 1)
lines <- c("# panglm validation report", "",
           sprintf("Generated %s with panglm %s and R %s.", format(Sys.Date()),
                   utils::packageVersion("panglm"), getRversion()), "",
           "Each row compares panglm with a reference package on the same model and data.",
           "`max abs diff` is the largest absolute difference over the compared values.", "",
           "| Model | Data | Quantity | Reference | Values | Max abs diff | Result |",
           "|---|---|---|---|---|---|---|",
           sprintf("| %s | %s | %s | %s | %d | %s | %s |", agg$model, agg$data, agg$quantity,
                   agg$package, agg$n, fmt(agg$absdiff), ifelse(agg$fail > 0, "FAIL", "pass")),
           "", "## Scripts", "",
           sprintf("- `%s`: %s", basename(scripts), status))
refs <- sort(unique(tab$package))
vers <- vapply(refs, function(p) tryCatch(as.character(utils::packageVersion(p)), error = function(e) "(not a package)"), "")
lines <- c(lines, "", "## Reference package versions", "", sprintf("- %s %s", refs, vers))
writeLines(lines, "validation/REPORT.md")
if (any(status != "pass")) quit(status = 1)
