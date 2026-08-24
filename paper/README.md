# JSS manuscript

This directory contains the JSS-style manuscript and standalone replication
materials for `panglm`.

Files:

- `article.Rnw`: executable R/LaTeX manuscript.
- `refs.bib`: JSS-formatted bibliography.
- `replication.R`: standalone replication of the manuscript's numerical and
  graphical results.
- `build.R`: installs the package into a temporary library, weaves and tangles
  the manuscript, and compiles the PDF.
- `jss.cls`, `jss.bst`, and `jsslogo.jpg`: official JSS template assets.

From the repository root, build the submission PDF with:

```sh
Rscript paper/build.R
```

Run the standalone replication script after installing `panglm` and `plm`:

```sh
Rscript paper/replication.R
```

The build requires `knitr`, a working LaTeX installation, and BibTeX.
