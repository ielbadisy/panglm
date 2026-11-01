# Build the `docvisits` dataset: German health-care demand panel
# (Riphahn, Wambach and Million 2003, GSOEP 1984-1994).
#
# Source (original file)
#   Journal of Applied Econometrics Data Archive, legacy site:
#   http://qed.econ.queensu.ca/jae/2003-v18.4/riphahn-wambach-million/
#   files readme.rwm.txt and rwm-data.zip (rwm.data, 27,326 rows x 25 columns).
#   New archive (ZBW): https://journaldata.zbw.eu/dataset/c9622542-10d8-4eda-b296-3ec48edce6c4
#   (DOI 10.15456/jae.2022314.1312212303; page protected by a bot wall, licence
#   field could not be read on 2026-11-01).
#
# Licence / redistribution basis (verified 2026-11-01)
#   * The identical 27,326 x 25 data are already distributed on CRAN under
#     GPL (>= 2) as `Health` in package Rchoice (0.3-6), whose help page cites
#     the source: "Riphahn, R. T., Wambach, A., & Million, A. (2003). Incentive
#     effects in the demand for health care ..." (man/Health.Rd,
#     https://cran.r-project.org/package=Rchoice, DESCRIPTION "License: GPL (>= 2)").
#     This script checks that our build from the JAE file equals that copy.
#   * Subsets are also on CRAN under GPL-2 in package COUNT (rwm, rwm5yr,
#     rwm1984; https://cran.r-project.org/package=COUNT).
#   * The JAE readme (readme.rwm.txt) states: "The distributors of the data
#     require that the data are only used for the purpose of scientific
#     studies." The data are shipped here for teaching and research use only,
#     and the help page repeats this condition.
#
# Cleaning (documented in ?docvisits)
#   * hsat: 40 imputed non-integer values (6.58 to 6.91) set to 7, as in
#     Greene (2008, "Modeling ordered choices", section 3.5) and Rchoice hsat2.
#   * handdum -> handicap: the readme (2013 addendum, item 5) says handdum is
#     swapped (0 <-> 1) in 1987; we recode 1987 as 1 - handdum. The 18 imputed
#     (gender mean) values are rounded to the nearest integer.
#   * handper (degree of handicap) dropped: erroneously imputed in 1987
#     (readme item 5).
#   * The five schooling-degree dummies are dropped (they overlap: 1,947 rows
#     have two degrees flagged, 556 none); `educ` (years) is kept.
#   * hhninc: the file holds monthly net household income in DM (the readme's
#     "/ 1000" does not match the file); we store it in 1000 DM.
#     Greene's examples use DM / 10000, so his income coefficients are 10
#     times ours.
#   * doctor = 1(docvis > 0), hospital = 1(hospvis > 0) added as in Greene.
#   No rows are dropped and no sampling is done.

raw_dir <- file.path("..", "..", "archive", "raw", "datasets-2026-11-01", "docvisits")
# absolute fallback used when building from the maintainer's machine
if (!dir.exists(raw_dir)) {
  raw_dir <- "/home/imad-el-badisy/Desktop/BIOSTAT/PROJECTS/SOFTWARE/project-panglm/archive/raw/datasets-2026-11-01/docvisits"
}
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

jae <- "http://qed.econ.queensu.ca/jae/2003-v18.4/riphahn-wambach-million/"
zipf <- file.path(raw_dir, "rwm-data.zip")
if (!file.exists(zipf)) {
  # the legacy JAE server has an incomplete certificate chain; plain http works
  utils::download.file(paste0(jae, "rwm-data.zip"), zipf, mode = "wb")
  utils::download.file(paste0(jae, "readme.rwm.txt"),
                       file.path(raw_dir, "readme.rwm.txt"), mode = "wb")
}
txt <- readLines(unz(zipf, "rwm.data"))
raw <- utils::read.table(text = txt)
names(raw) <- c("id", "female", "year", "age", "hsat", "handdum", "handper",
                "hhninc", "hhkids", "educ", "married", "haupts", "reals",
                "fachhs", "abitur", "univ", "working", "bluec", "whitec",
                "self", "beamt", "docvis", "hospvis", "public", "addon")
stopifnot(nrow(raw) == 27326L, ncol(raw) == 25L,
          length(unique(raw$id)) == 7293L,
          !anyDuplicated(raw[, c("id", "year")]),
          !anyNA(raw))

# ---------------------------------------------------------------------------
# Cross-check 1: identical to Rchoice::Health (CRAN, GPL >= 2)
# ---------------------------------------------------------------------------
rch_tar <- file.path(raw_dir, "Rchoice_0.3-6.tar.gz")
if (!file.exists(rch_tar)) {
  for (u in c("https://cran.r-project.org/src/contrib/Rchoice_0.3-6.tar.gz",
              "https://cran.r-project.org/src/contrib/Archive/Rchoice/Rchoice_0.3-6.tar.gz")) {
    ok <- tryCatch(utils::download.file(u, rch_tar, mode = "wb", quiet = TRUE) == 0,
                   error = function(e) FALSE, warning = function(w) FALSE)
    if (isTRUE(ok)) break
  }
}
if (file.exists(rch_tar)) {
  ex <- tempfile()
  utils::untar(rch_tar, files = "Rchoice/data/Health.rda", exdir = ex)
  e <- new.env()
  load(file.path(ex, "Rchoice", "data", "Health.rda"), envir = e)
  H <- e$Health
  names(H)[names(H) == "hhinc"] <- "hhninc"
  for (v in names(raw)) {
    # Rchoice stores some columns in single precision (Stata float)
    stopifnot(isTRUE(all.equal(as.numeric(raw[[v]]), as.numeric(H[[v]]),
                               tolerance = 1e-6)))
  }
  cat("Cross-check: JAE file == Rchoice::Health on all 25 columns\n")
}

# ---------------------------------------------------------------------------
# Clean
# ---------------------------------------------------------------------------
hsat <- raw$hsat
stopifnot(sum(hsat != round(hsat)) == 40L, all(hsat[hsat != round(hsat)] > 6.5))
hsat[hsat != round(hsat)] <- 7

handicap <- ifelse(raw$year == 1987, 1 - raw$handdum, raw$handdum)
handicap <- round(handicap)

docvisits <- data.frame(
  id          = as.integer(raw$id),
  year        = as.integer(raw$year),
  docvis      = as.integer(raw$docvis),
  hospvis     = as.integer(raw$hospvis),
  doctor      = as.integer(raw$docvis > 0),
  hospital    = as.integer(raw$hospvis > 0),
  hsat        = as.integer(hsat),
  female      = as.integer(raw$female),
  age         = as.integer(raw$age),
  educ        = round(raw$educ, 4),
  married     = as.integer(raw$married),
  hhkids      = as.integer(raw$hhkids),
  hhninc      = round(raw$hhninc / 1000, 5),
  working     = as.integer(raw$working),
  bluecollar  = as.integer(raw$bluec),
  whitecollar = as.integer(raw$whitec),
  self        = as.integer(raw$self),
  civil       = as.integer(raw$beamt),
  handicap    = as.integer(handicap),
  public      = as.integer(raw$public),
  addon       = as.integer(raw$addon)
)
docvisits <- docvisits[order(docvisits$id, docvisits$year), ]
rownames(docvisits) <- NULL

# ---------------------------------------------------------------------------
# Validation against published numbers
# ---------------------------------------------------------------------------
chk <- function(what, ours, ref, tol, where) {
  data.frame(check = what, ours = signif(ours, 7), reference = ref,
             tol = tol, ok = abs(ours - ref) <= tol, source = where)
}
res <- list()

# (a) RWM (2003) Table 1, corrected values for male hospital visits given in
#     the JAE readme (2013 addendum, item 2): percent of male person-years
#     with 3 visits = 0.15 and with 4-9 visits = 0.203.
m <- docvisits[docvisits$female == 0, ]
res[[1]] <- chk("RWM T1 men: % hospvis == 3", 100 * mean(m$hospvis == 3),
                0.15, 0.005, "JAE readme.rwm.txt, item (2)")
res[[2]] <- chk("RWM T1 men: % hospvis 4-9",
                100 * mean(m$hospvis >= 4 & m$hospvis <= 9),
                0.203, 0.0015, "JAE readme.rwm.txt, item (2)")

# (b) Greene, Discrete Choice Modeling 2014, Part 3 (Panel data binary
#     choice), slides 19 and 28: pooled probit of DOCTOR on AGE, EDUC, HHNINC
#     (HHNINC in DM/10000; mean .35208); N = 27,326, log L = -17701.08500.
#     https://pages.stern.nyu.edu/wgreene/DiscreteChoice/2014/DC2014-Part3-PanelDataBinaryChoice.pdf
g <- docvisits
g$inc10 <- g$hhninc / 10
pp <- stats::glm(doctor ~ age + educ + inc10, family = binomial("probit"), data = g)
src_b <- "Greene DC2014 Part 3, slide 19/28"
res[[3]] <- chk("mean HHNINC (DM/10000)", mean(g$inc10), 0.35208, 5e-6, src_b)
res[[4]] <- chk("pooled probit: Constant", coef(pp)[[1]], 0.02159, 5e-5, src_b)
res[[5]] <- chk("pooled probit: AGE", coef(pp)[["age"]], 0.01532, 5e-6, src_b)
res[[6]] <- chk("pooled probit: EDUC", coef(pp)[["educ"]], -0.02793, 5e-6, src_b)
res[[7]] <- chk("pooled probit: HHNINC", coef(pp)[["inc10"]], -0.10204, 5e-5, src_b)
res[[8]] <- chk("pooled probit: logL", as.numeric(logLik(pp)), -17701.08500, 0.01, src_b)

# (c) Same source, slide 23: pooled probit adding FEMALE, log L = -17457.21899.
pf <- stats::glm(doctor ~ age + educ + inc10 + female,
                 family = binomial("probit"), data = g)
src_c <- "Greene DC2014 Part 3, slide 23"
res[[9]] <- chk("pooled probit + FEMALE: FEMALE", coef(pf)[["female"]], 0.35209, 5e-5, src_c)
res[[10]] <- chk("pooled probit + FEMALE: logL", as.numeric(logLik(pf)), -17457.21899, 0.01, src_c)

# (d) Greene (2008), "Modeling Ordered Choices", Table 3.1: ordered logit of
#     health satisfaction collapsed to 5 levels ((0-2)=0, (3-5)=1, (6-8)=2,
#     9=3, 10=4), year 1988 only (n = 4,483), x = (Age, Income, Education,
#     Married, Kids); log L = -5749.157, AGE -.0321, EDUC .0645, INCOME .4263.
#     https://pages.stern.nyu.edu/wgreene/DiscreteChoice/ModelsForOrderedChoices-Handbook.pdf
d88 <- g[g$year == 1988, ]
d88$health <- cut(d88$hsat, c(-Inf, 2, 5, 8, 9, 10), labels = FALSE) - 1L
ol <- MASS::polr(factor(health) ~ age + inc10 + educ + married + hhkids,
                 data = d88, method = "logistic")
src_d <- "Greene (2008) Modeling Ordered Choices, Table 3.1"
res[[11]] <- chk("1988 n", nrow(d88), 4483, 0, src_d)
res[[12]] <- chk("ologit 1988: AGE", coef(ol)[["age"]], -0.0321, 5e-5, src_d)
res[[13]] <- chk("ologit 1988: EDUC", coef(ol)[["educ"]], 0.0645, 5e-5, src_d)
res[[14]] <- chk("ologit 1988: INCOME", coef(ol)[["inc10"]], 0.4263, 5e-4, src_d)
res[[15]] <- chk("ologit 1988: logL", as.numeric(logLik(ol)), -5749.157, 0.01, src_d)

# (e) Cross-check 2: COUNT::rwm (CRAN, GPL-2) has docvis, age, educ and
#     hhninc (in 1000 DM) for the same 27,326 rows, in the same order.
cnt_tar <- file.path(raw_dir, "COUNT_1.3.5.tar.gz")
if (!file.exists(cnt_tar)) {
  for (u in c("https://cran.r-project.org/src/contrib/COUNT_1.3.5.tar.gz",
              "https://cran.r-project.org/src/contrib/Archive/COUNT/COUNT_1.3.5.tar.gz")) {
    ok <- tryCatch(utils::download.file(u, cnt_tar, mode = "wb", quiet = TRUE) == 0,
                   error = function(e) FALSE, warning = function(w) FALSE)
    if (isTRUE(ok)) break
  }
}
if (file.exists(cnt_tar)) {
  ex <- tempfile()
  utils::untar(cnt_tar, files = "COUNT/data/rwm.rda", exdir = ex)
  e <- new.env()
  load(file.path(ex, "COUNT", "data", "rwm.rda"), envir = e)
  rwm <- e$rwm
  stopifnot(nrow(rwm) == nrow(raw),
            all(as.integer(rwm$docvis) == raw$docvis),
            all(as.integer(rwm$age) == raw$age),
            max(abs(as.numeric(rwm$hhninc) - raw$hhninc / 1000)) < 1e-4)
  cat("Cross-check: docvis, age, hhninc agree with COUNT::rwm\n")
}

# (f) Optional (about 1.5 minutes): Greene DC2014 Part 3, slide 28,
#     random-effects probit by quadrature: Constant -.11819, AGE .02232,
#     EDUC -.03307, HHNINC .00660, Rho .44990, log L = -16290.72192.
#     Run with PANGLM_FULL_VALIDATION=1 and panglm installed.
if (identical(Sys.getenv("PANGLM_FULL_VALIDATION"), "1") &&
    requireNamespace("panglm", quietly = TRUE)) {
  rp <- panglm::panglm(doctor ~ age + educ + inc10, data = g,
                       index = c("id", "year"), model = "random",
                       family = "probit")
  src_f <- "Greene DC2014 Part 3, slide 28"
  rho <- rp$sigma^2 / (1 + rp$sigma^2)
  res[[length(res) + 1]] <- chk("RE probit: AGE", coef(rp)[["age"]], 0.02232, 5e-5, src_f)
  res[[length(res) + 1]] <- chk("RE probit: EDUC", coef(rp)[["educ"]], -0.03307, 5e-5, src_f)
  res[[length(res) + 1]] <- chk("RE probit: rho", rho, 0.44990, 5e-4, src_f)
  res[[length(res) + 1]] <- chk("RE probit: logL", as.numeric(logLik(rp)), -16290.72192, 0.05, src_f)
}

val <- do.call(rbind, res)
print(val, row.names = FALSE)
stopifnot(all(val$ok))

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
save(docvisits, file = file.path("data", "docvisits.rda"), compress = "xz")
cat("Saved data/docvisits.rda:", file.size(file.path("data", "docvisits.rda")), "bytes\n")
