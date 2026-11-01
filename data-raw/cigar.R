# Build the `cigar` dataset: Baltagi and Levin US state cigarette demand
# panel, 46 states (45 states plus the District of Columbia), 1963-1992.
#
# Source
#   Original file `Cigar.txt` from the data archive of Baltagi, Econometric
#   Analysis of Panel Data (Wiley companion web site):
#   https://www.wiley.com/legacy/wileychi/baltagi/supp/Cigar.txt
#   listed on https://www.wiley.com/legacy/wileychi/baltagi/datasets.html
#   ("Panel Data, 46 U.S. States over the period 1963-1992"). Variable
#   description from https://www.wiley.com/legacy/wileychi/baltagi/supp/readme.doc:
#   "CIGAR.TXT Source: Baltagi and Levin (1992) and Baltagi, Griffin and
#   Xiong (2000). Description: Panel Data, 46 U.S. States over the period
#   1963-1992." Downloaded 2026-11-01 (1380 lines, md5
#   17318a02974301c0e714983dab02c139).
#
# Redistribution (verified 2026-11-01)
#   The Wiley page states no licence. The same 1380 x 9 table has been
#   redistributed for many years in two CRAN packages licensed
#   "GPL (>= 2)": plm (data `Cigar`, CRAN page
#   https://cran.r-project.org/package=plm, "License: GPL-2 | GPL-3
#   [expanded from: GPL (>= 2)]"; help page source: "Online complements to
#   Baltagi (2001): https://www.wiley.com/legacy/wileychi/baltagi/") and
#   Ecdat (data `Cigar`, https://cran.r-project.org/package=Ecdat, same
#   licence; help page source: Baltagi and Levin (1992), Baltagi, Griffin
#   and Xiong (2000)). The values are public statistics compiled from the
#   Tobacco Institute's "The Tax Burden on Tobacco", US Census population
#   estimates and BLS/BEA price and income series. The build checks below
#   that the shipped copy is identical to both CRAN copies.
#
# Raw files are cached outside the package (never committed) in `rawdir`.

rawdir <- Sys.getenv("CIGAR_RAW", "../archive/raw/datasets-2026-11-01/cigar")
dir.create(rawdir, recursive = TRUE, showWarnings = FALSE)
rawfile <- file.path(rawdir, "Cigar.txt")
if (!file.exists(rawfile)) {
  utils::download.file(
    "https://www.wiley.com/legacy/wileychi/baltagi/supp/Cigar.txt",
    rawfile, mode = "wb", quiet = TRUE
  )
}

raw <- utils::read.table(
  rawfile,
  col.names = c("code", "year", "price", "pop", "pop16", "cpi", "ndi",
                "sales", "pimin")
)
stopifnot(nrow(raw) == 1380, ncol(raw) == 9)

## ---- state codes -----------------------------------------------------------
## Baltagi numbers the 50 states plus DC alphabetically by name (1 = Alabama,
## 9 = District of Columbia, 51 = Wyoming). Codes 2 (Alaska), 6 (Colorado),
## 12 (Hawaii), 34 (North Carolina) and 38 (Oregon) are absent. Checked
## against the textbook text (Baltagi 2005, p. 157): Utah (45) sold 55 packs
## per head in 1988 and Nevada (29) 142 packs.
abb <- c("AL", "AK", "AZ", "AR", "CA", "CO", "CT", "DE", "DC", "FL", "GA",
         "HI", "ID", "IL", "IN", "IA", "KS", "KY", "LA", "ME", "MD", "MA",
         "MI", "MN", "MS", "MO", "MT", "NE", "NV", "NH", "NJ", "NM", "NY",
         "NC", "ND", "OH", "OK", "OR", "PA", "RI", "SC", "SD", "TN", "TX",
         "UT", "VT", "VA", "WA", "WV", "WI", "WY")
stopifnot(setequal(setdiff(1:51, raw$code), c(2, 6, 12, 34, 38)))
stopifnot(raw$sales[raw$code == 45 & raw$year == 88] == 55,
          round(raw$sales[raw$code == 29 & raw$year == 88]) == 142)

## ---- assemble --------------------------------------------------------------
present <- sort(unique(raw$code))
cigar <- data.frame(
  state = factor(abb[raw$code], levels = abb[present]),
  year  = as.integer(1900L + raw$year),
  price = raw$price,
  pop   = raw$pop,
  pop16 = raw$pop16,
  cpi   = raw$cpi,
  ndi   = raw$ndi,
  sales = raw$sales,
  pimin = raw$pimin
)
## Baltagi's variables (textbook eq. 8.48): C = real per capita sales by
## persons of smoking age, P = real price, Pn = minimum real price in
## neighbouring states, Y = real per capita disposable income. The source
## column `sales` is packs per head of total population; C rescales it to
## the population aged 16 and over (`pop16`), which reproduces the
## published coefficients exactly (see validation).
cigar$lnc  <- log(cigar$sales * cigar$pop / cigar$pop16)
cigar$lnp  <- log(cigar$price / cigar$cpi * 100)
cigar$lnpn <- log(cigar$pimin / cigar$cpi * 100)
cigar$lny  <- log(cigar$ndi / cigar$cpi * 100)
cigar <- cigar[order(cigar$state, cigar$year), ]
rownames(cigar) <- NULL

## ---- validation 1: identical to the plm and Ecdat copies --------------------
srcmat <- unname(as.matrix(raw[order(raw$code, raw$year), ]))
if (requireNamespace("plm", quietly = TRUE)) {
  e <- new.env(); utils::data("Cigar", package = "plm", envir = e)
  stopifnot(isTRUE(all.equal(srcmat, unname(as.matrix(e$Cigar)))))
  cat("identical to plm::Cigar\n")
}
if (requireNamespace("Ecdat", quietly = TRUE)) {
  e <- new.env(); utils::data("Cigar", package = "Ecdat", envir = e)
  stopifnot(isTRUE(all.equal(srcmat, unname(as.matrix(e$Cigar)))))
  cat("identical to Ecdat::Cigar\n")
}

## ---- validation 2: Baltagi (2005) Table 8.1 ---------------------------------
## Baltagi, B.H. (2005) Econometric Analysis of Panel Data, 3rd ed., Wiley,
## Table 8.1 "Pooled Estimation Results. Cigarette Demand Equation 1963-92",
## p. 157 (same table in later editions, section "Empirical Example: Dynamic
## Demand for Cigarettes"). Columns ln C(t-1), ln P, ln Pn, ln Y; "All
## regressions except OLS and 2SLS include time dummies". Published values
## are rounded to 2-3 decimals; t statistics in parentheses.
cigar$lnclag <- ave(cigar$lnc, cigar$state, FUN = function(v) c(NA, v[-length(v)]))
ols <- stats::lm(lnc ~ lnclag + lnp + lnpn + lny, data = cigar)
wit <- stats::lm(lnc ~ lnclag + lnp + lnpn + lny + factor(state) + factor(year),
                 data = cigar)
tstat <- function(m) (coef(m) / sqrt(diag(stats::vcov(m))))[2:5]
gmm1 <- panglm::pangmm(lnc ~ lag(lnc, 1) + lnp + lnpn + lny | lag(lnc, 2:99),
                       data = cigar, index = c("state", "year"),
                       transformation = "d", model = "onestep")
cigar$lnclag <- NULL

tab81 <- data.frame(
  estimator = rep(c("OLS", "Within", "GMM-one-step"), each = 4),
  term = rep(c("lnC(-1)", "lnP", "lnPn", "lnY"), 3),
  published = c(0.97, -0.090, 0.024, -0.03,
                0.83, -0.299, 0.034, 0.10,
                0.84, -0.377, -0.016, 0.14),
  digits = rep(c(2, 3, 3, 2), 3),
  ours = c(coef(ols)[2:5], coef(wit)[2:5], coef(gmm1)[1:4])
)
tab81$ok <- round(tab81$ours, tab81$digits) == tab81$published
print(transform(tab81, ours = round(ours, 4)), row.names = FALSE)
stopifnot(all(tab81$ok))
## t statistics of OLS and Within, published (157.7, 6.2, 1.8, 5.1) and
## (66.3, 12.7, 1.2, 4.2) in absolute value
tpub <- c(157.7, 6.2, 1.8, 5.1, 66.3, 12.7, 1.2, 4.2)
tours <- abs(c(tstat(ols), tstat(wit)))
print(round(cbind(published = tpub, ours = tours), 2))
stopifnot(all(abs(round(tours, 1) - tpub) <= 0.1))

save(cigar, file = "data/cigar.rda", compress = "xz")
