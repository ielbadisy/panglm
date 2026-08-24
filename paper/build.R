arguments <- commandArgs(trailingOnly = FALSE)
file_argument <- arguments[grepl("^--file=", arguments)]
script_path <- if (length(file_argument)) {
  normalizePath(sub("^--file=", "", file_argument[1]))
} else {
  normalizePath("paper/build.R")
}
paper_directory <- dirname(script_path)
package_directory <- dirname(paper_directory)
library_directory <- tempfile("panglm-paper-library-")
dir.create(library_directory)

status <- system2(
  file.path(R.home("bin"), "R"),
  c("CMD", "INSTALL", paste0("--library=", library_directory),
    shQuote(package_directory))
)
if (status != 0L) stop("package installation failed", call. = FALSE)

.libPaths(c(library_directory, .libPaths()))
old_directory <- setwd(paper_directory)
on.exit(setwd(old_directory), add = TRUE)

knitr::render_sweave()
knitr::knit(
  "article.Rnw", output = "panglm-jss.tex",
  quiet = TRUE, envir = new.env(parent = globalenv())
)
knitr::purl("article.Rnw", output = "article-tangled.R", quiet = TRUE)
tools::texi2pdf("panglm-jss.tex", clean = FALSE)
