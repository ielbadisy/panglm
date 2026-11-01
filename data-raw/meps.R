# MEPS Panel 23 (2018-2019) two-year adult person panel.
#
# Source: Agency for Healthcare Research and Quality (AHRQ), Medical
# Expenditure Panel Survey (MEPS) Household Component, public use file
# HC-217 "MEPS Panel 23 Longitudinal Data File" (release December 2021).
#   File page: https://meps.ahrq.gov/mepsweb/data_stats/download_data_files_detail.jsp?cboPufNumber=HC-217
#   Data (Stata): https://meps.ahrq.gov/mepsweb/data_files/pufs/h217/h217dta.zip
#   Documentation: https://meps.ahrq.gov/mepsweb/data_stats/download_data/pufs/h217/h217doc.shtml
#   Codebook: https://meps.ahrq.gov/mepsweb/data_stats/download_data_files_codebook.jsp?PUFId=H217
#
# Licence / terms. MEPS is produced by AHRQ, a US federal agency; the public
# use files carry no copyright statement and no redistribution clause. They
# are released under the following Data Use Agreement (section A of the
# HC-217 documentation, verified 2026-11-01 at the documentation URL above):
#
#   "Individual identifiers have been removed from the micro-data contained
#   in these files. Nevertheless, under sections 308 (d) and 903 (c) of the
#   Public Health Service Act (42 U.S.C. 242m and 42 U.S.C. 299 a-1), data
#   collected by the Agency for Healthcare Research and Quality (AHRQ)
#   and/or the National Center for Health Statistics (NCHS) may not be used
#   for any purpose other than for the purpose for which they were supplied;
#   any effort to determine the identity of any reported cases is
#   prohibited by law.
#   Therefore in accordance with the above referenced Federal Statute, it is
#   understood that:
#   1. No one is to use the data in this data set in any way except for
#      statistical reporting and analysis; and
#   2. If the identity of any person or establishment should be discovered
#      inadvertently, then (a) no use will be made of this knowledge, (b)
#      the Director Office of Management AHRQ will be advised of this
#      incident, (c) the information that would identify any individual or
#      establishment will be safeguarded or destroyed, as requested by
#      AHRQ, and (d) no one else will be informed of the discovered
#      identity; and
#   3. No one will attempt to link this data set with individually
#      identifiable records from any data sets other than the Medical
#      Expenditure Panel Survey or the National Health Interview Survey.
#   By using these data you signify your agreement to comply with the above
#   stated statutorily based requirements [...]
#   The Agency for Healthcare Research and Quality requests that users cite
#   AHRQ and the Medical Expenditure Panel Survey as the data source in any
#   publications or research based upon these data."
#
# The agreement restricts use (statistical analysis only, no
# re-identification, no linkage) but does not prohibit redistribution. The
# shipped extract therefore (i) replaces DUPERSID by a new sequential id so
# rows cannot be matched back to the public file or linked to NHIS, (ii)
# keeps only a random subsample of persons and a small set of variables,
# and (iii) reproduces the agreement in the help page (?meps).
#
# Build steps
#   1. download h217dta.zip if absent (directory: env var PANGLM_RAW,
#      default tempdir()) and read h217.dta with haven;
#   2. exact validation of the FULL file against the HC-217 codebook
#      (unweighted frequencies) and documentation (sum of LONGWT);
#   3. official validation: LONGWT-weighted mean total expenditure and
#      percent with any expense per year vs AHRQ Statistical Briefs #533
#      (2018) and #540 (2019);
#   4. restrict to persons with data for all five rounds (ALL5RDS = 1, so
#      in scope for the whole of 2018 and 2019) aged 18+ in 2018, draw a
#      stratified random sample of 3,000 persons (strata = 2018 age group,
#      proportional allocation, set.seed(20261101)), reshape to long format
#      (2 rows per person) and save data/meps.rda.

rawdir <- Sys.getenv("PANGLM_RAW", tempdir())
dir.create(rawdir, showWarnings = FALSE, recursive = TRUE)
zipfile <- file.path(rawdir, "h217dta.zip")
dtafile <- file.path(rawdir, "h217.dta")
if (!file.exists(dtafile)) {
  if (!file.exists(zipfile)) {
    utils::download.file(
      "https://meps.ahrq.gov/mepsweb/data_files/pufs/h217/h217dta.zip",
      zipfile, mode = "wb")
  }
  utils::unzip(zipfile, exdir = rawdir)
}
h <- as.data.frame(haven::zap_labels(haven::read_dta(dtafile)))

# ---------------------------------------------------------------------------
# 1. Exact checks on the full file (unweighted counts from the HC-217
#    codebook, variable pages; sum of weights from documentation 2.1.3).
# ---------------------------------------------------------------------------
cnt <- function(x, v) as.integer(table(factor(x, levels = v)))
exact <- list(
  nrow = c(14067L, nrow(h)),
  all5rds1 = c(13044L, sum(h$ALL5RDS == 1)),
  inscovy1 = rbind(ref = c(161L, 8386L, 4554L, 966L),
                   got = cnt(h$INSCOVY1, c(-1, 1, 2, 3))),
  racethx = rbind(ref = c(3299L, 7427L, 2081L, 725L, 535L),
                  got = cnt(h$RACETHX, 1:5)),
  povcaty2 = rbind(ref = c(140L, 2304L, 690L, 2044L, 4075L, 4814L),
                   got = cnt(h$POVCATY2, c(-1, 1:5))),
  totexpy1 = rbind(ref = c(161L, 1818L), got = cnt(h$TOTEXPY1, c(-1, 0))),
  obtotvy2 = rbind(ref = c(140L, 3646L), got = cnt(h$OBTOTVY2, c(-1, 0)))
)
print(exact)
stopifnot(all(vapply(exact, function(m) {
  m <- matrix(m, nrow = 2); all(m[1, ] == m[2, ])
}, logical(1))))
# Weighted totals: codebook TOTAL row 330,644,961 (all persons);
# documentation 2.1.3: LONGWT over the 13,044 ALL5RDS = 1 cases gives
# "305.7 million".
wsum <- c(all = sum(h$LONGWT), all5 = sum(h$LONGWT[h$ALL5RDS == 1]))
print(wsum)
stopifnot(abs(wsum[["all"]] - 330644961) < 1,
          abs(wsum[["all5"]] / 1e6 - 305.7) < 0.05)

# ---------------------------------------------------------------------------
# 2. Official estimates. AHRQ Statistical Brief #533 (Mitchell, Jan 2021),
#    2018: overall mean expenditure per person $6,063 (figure 2), "about 13
#    percent ... had no personal healthcare expenditures". Statistical Brief
#    #540 (Mitchell, Feb 2022), 2019: overall mean $6,252 (figure 2), "About
#    14 percent ... had no personal healthcare expenditures".
#    These come from the annual full-year files (HC-209, HC-216) with
#    annual person weights, pooling two panels and covering everyone in
#    scope at any time in the year. Here: panel 23 only, LONGWT, persons
#    in scope for the whole of 2018-2019 (ALL5RDS = 1). That excludes 2018
#    decedents and people who entered later, so the 2018 mean is expected
#    to be somewhat lower. Tolerance: 10 percent relative for the mean,
#    2.5 percentage points for the share with any expense (the briefs
#    round it to whole percents).
# ---------------------------------------------------------------------------
f5 <- h[h$ALL5RDS == 1, ]
off <- data.frame(
  year = c(2018L, 2019L),
  mean_ref = c(6063, 6252),
  mean_got = c(weighted.mean(f5$TOTEXPY1, f5$LONGWT),
               weighted.mean(f5$TOTEXPY2, f5$LONGWT)),
  pany_ref = c(87, 86),
  pany_got = 100 * c(weighted.mean(f5$TOTEXPY1 > 0, f5$LONGWT),
                     weighted.mean(f5$TOTEXPY2 > 0, f5$LONGWT))
)
off$mean_reldiff <- off$mean_got / off$mean_ref - 1
print(off, digits = 4)
stopifnot(all(abs(off$mean_reldiff) < 0.10),
          all(abs(off$pany_got - off$pany_ref) < 2.5))

# ---------------------------------------------------------------------------
# 3. Sample: adults (18+ in 2018) with all five rounds, stratified by age.
# ---------------------------------------------------------------------------
pop <- f5[f5$AGEY1X >= 18 & f5$LONGWT > 0, ]
pop <- pop[order(pop$DUPERSID), ]
pop$agegrp <- cut(pop$AGEY1X, c(17, 34, 49, 64, Inf))
n_target <- 3000L
alloc <- round(n_target * table(pop$agegrp) / nrow(pop))
set.seed(20261101)
keep <- unlist(lapply(names(alloc), function(g) {
  ix <- which(pop$agegrp == g)
  ix[sort(sample.int(length(ix), alloc[[g]]))]
}))
s <- pop[sort(keep), ]
cat("persons sampled:", nrow(s), "\n"); print(rbind(pop = table(pop$agegrp),
                                                    sample = alloc))

neg2na <- function(x) replace(x, x < 0, NA)
chronic <- function(d, y) {
  v <- c("HIBPDX", "CHDDX", "STRKDX", "EMPHDX", "CHOLDX", "DIABDX", "ARTHDX",
         "ASTHDX")
  nm <- paste0(v, "Y", y)
  nm[v == "DIABDX"] <- paste0("DIABDXY", y, "_M18")
  as.integer(rowSums(sapply(nm, function(k) d[[k]] == 1)))
}
one_year <- function(d, y) {
  yy <- c("Y1", "Y2")[y]
  g <- function(v) d[[paste0(v, yy)]]
  hl <- neg2na(d[[c("RTHLTH3", "RTHLTH5")[y]]])
  data.frame(
    id = seq_len(nrow(d)),
    year = c(2018L, 2019L)[y],
    totexp = as.numeric(g("TOTEXP")),
    oopexp = as.numeric(g("TOTSLF")),
    anyexp = as.integer(g("TOTEXP") > 0),
    officevis = as.integer(g("OBTOTV")),
    ervis = as.integer(g("ERTOT")),
    ipstays = as.integer(g("IPDIS")),
    insurance = factor(g("INSCOV"), 1:3,
                       c("private", "public", "uninsured")),
    age = as.integer(d[[paste0("AGE", yy, "X")]]),
    sex = factor(d$SEX, 1:2, c("male", "female")),
    race = factor(d$RACETHX, 1:5,
                  c("hispanic", "white", "black", "asian", "other")),
    region = factor(neg2na(g("REGION")), 1:4,
                    c("northeast", "midwest", "south", "west")),
    povcat = factor(g("POVCAT"), 1:5,
                    c("poor", "nearpoor", "low", "middle", "high")),
    povlev = round(as.numeric(g("POVLEV")), 1),
    health = factor(hl, 1:5, c("excellent", "verygood", "good", "fair",
                               "poor"), ordered = TRUE),
    nchronic = chronic(d, y),
    longwt = as.numeric(d$LONGWT),
    varstr = as.integer(d$VARSTR),
    varpsu = as.integer(d$VARPSU)
  )
}
meps <- rbind(one_year(s, 1), one_year(s, 2))
meps <- meps[order(meps$id, meps$year), ]
rownames(meps) <- NULL
str(meps)
stopifnot(nrow(meps) == 2 * n_target, !anyDuplicated(meps[c("id", "year")]),
          !anyNA(meps[setdiff(names(meps), c("health", "region"))]),
          all(meps$totexp >= meps$oopexp))

# Loose check: sample vs full-file adults (same restriction), weighted.
# The sample is self-weighting within age strata only approximately, so
# a 15 percent tolerance is used for the mean and 3 points for the share.
loose <- data.frame(
  year = c(2018L, 2019L),
  mean_pop = c(weighted.mean(pop$TOTEXPY1, pop$LONGWT),
               weighted.mean(pop$TOTEXPY2, pop$LONGWT)),
  mean_smp = tapply(meps$totexp * meps$longwt, meps$year, sum) /
    tapply(meps$longwt, meps$year, sum),
  pany_pop = 100 * c(weighted.mean(pop$TOTEXPY1 > 0, pop$LONGWT),
                     weighted.mean(pop$TOTEXPY2 > 0, pop$LONGWT)),
  pany_smp = 100 * tapply(meps$anyexp * meps$longwt, meps$year, sum) /
    tapply(meps$longwt, meps$year, sum)
)
print(loose, digits = 4)
stopifnot(all(abs(loose$mean_smp / loose$mean_pop - 1) < 0.15),
          all(abs(loose$pany_smp - loose$pany_pop) < 3))

save(meps, file = "data/meps.rda", compress = "xz")
cat("size (KB):", file.size("data/meps.rda") / 1024, "\n")
