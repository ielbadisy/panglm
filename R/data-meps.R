#' Health care expenditure and use of US adults, MEPS panel 23 (2018-2019)
#'
#' A real public two-year person panel from the US Medical Expenditure Panel
#' Survey (MEPS) Household Component, longitudinal file HC-217 (panel 23).
#' The unit is a person, the time index is the calendar year (2018, 2019):
#' each of 3,000 adults contributes exactly two rows. Expenditure outcomes
#' have the classic health-economics shape (a mass of zeros plus a long
#' right tail), utilisation outcomes are overdispersed counts with many
#' zeros, and insurance, income and chronic conditions vary over time
#' within person, which makes the data a natural testbed for two-part,
#' hurdle, tobit, gamma, quantile, selection and correlated random effects
#' (Mundlak) panel models.
#'
#' The extract keeps persons with data for all five interview rounds
#' (`ALL5RDS = 1`, in scope for the whole of 2018 and 2019) who were aged
#' 18 or older in 2018, and draws a stratified random sample of 3,000 of
#' the 10,084 eligible persons (strata: 2018 age 18-34, 35-49, 50-64, 65+;
#' proportional allocation; `set.seed(20261101)`). The original person
#' identifier `DUPERSID` is replaced by a new sequential `id` so that rows
#' cannot be matched back to the public file or linked to other sources.
#' The build script, with its validation against the codebook and against
#' official AHRQ estimates, is `data-raw/meps.R` in the package's GitHub
#' repository.
#'
#' Because the sample is drawn with equal probability within age strata,
#' `longwt` remains approximately valid for weighted descriptive
#' statistics of the adult population (the 2018 and 2019 weighted means of
#' `totexp` in the sample are within 2 percent of those of all eligible
#' adults in the full file), but the design variables `varstr` and
#' `varpsu` only describe the original MEPS design: design-based standard
#' errors computed from this subsample are approximate. The panglm
#' examples below are unweighted model-based fits, as usual for teaching.
#'
#' @format A data frame with 6,000 rows (3,000 persons x 2 years) and 20
#'   columns. The panel index pair is `c("id", "year")`.
#' \describe{
#'   \item{id}{integer person identifier (1-3000), the panel individual
#'     index; a new sequential id, not the MEPS `DUPERSID`}
#'   \item{year}{integer calendar year (2018 or 2019), the panel time
#'     index (MEPS suffix Y1, Y2)}
#'   \item{totexp}{numeric total health care expenditure in the year, US
#'     dollars, all payers (`TOTEXPY1`, `TOTEXPY2`); 14 percent zeros}
#'   \item{oopexp}{numeric out-of-pocket (self or family) expenditure, US
#'     dollars (`TOTSLFY1`, `TOTSLFY2`); never larger than `totexp`}
#'   \item{anyexp}{integer indicator (0/1) of any health care expenditure in
#'     the year, `as.integer(totexp > 0)`}
#'   \item{officevis}{integer number of office-based provider visits
#'     (`OBTOTVY1`, `OBTOTVY2`)}
#'   \item{ervis}{integer number of emergency room visits (`ERTOTY1`,
#'     `ERTOTY2`)}
#'   \item{ipstays}{integer number of hospital inpatient discharges
#'     (`IPDISY1`, `IPDISY2`)}
#'   \item{insurance}{factor health insurance coverage in the year
#'     (`INSCOVY1`, `INSCOVY2`): `private` (any private coverage),
#'     `public` (public coverage only, e.g. Medicare, Medicaid),
#'     `uninsured` (uninsured all year); reference level `private`}
#'   \item{age}{integer age in years at the end of the year (or last age in
#'     scope), top-coded at 85 (`AGEY1X`, `AGEY2X`). Age rises by one year
#'     for (almost) everyone, so within person it is collinear with the year:
#'     use it in pooled and random effects models but not as a time-varying
#'     regressor in within or Mundlak (`cre = TRUE`) fits}
#'   \item{sex}{factor `male`, `female` (`SEX`), time-invariant}
#'   \item{race}{factor race and ethnicity (`RACETHX`): `hispanic`,
#'     `white` (non-Hispanic white only), `black` (non-Hispanic black
#'     only), `asian` (non-Hispanic Asian only), `other` (non-Hispanic
#'     other or multiple race); time-invariant}
#'   \item{region}{factor Census region at the end of the year
#'     (`REGIONY1`, `REGIONY2`): `northeast`, `midwest`, `south`, `west`}
#'   \item{povcat}{factor family income as a percent of the federal poverty
#'     line, in categories (`POVCATY1`, `POVCATY2`): `poor` (below 100
#'     percent, including negative income), `nearpoor` (100 to below 125),
#'     `low` (125 to below 200), `middle` (200 to below 400), `high` (400
#'     or more)}
#'   \item{povlev}{numeric family income as a percent of the poverty line,
#'     one decimal (`POVLEVY1`, `POVLEVY2`); can be negative}
#'   \item{health}{ordered factor self-rated (perceived) health status,
#'     `excellent` < `verygood` < `good` < `fair` < `poor`, from the
#'     interview round closing the year (2018: round 3 `RTHLTH3`; 2019:
#'     round 5 `RTHLTH5`), as in the MEPS annual files; 4 values `NA`
#'     (refused, don't know or not ascertained)}
#'   \item{nchronic}{integer number (0-8) of eight MEPS priority conditions
#'     ever diagnosed as reported in that year: high blood pressure,
#'     coronary heart disease, stroke, emphysema, high cholesterol,
#'     diabetes, arthritis, asthma (`HIBPDXYt`, `CHDDXYt`, `STRKDXYt`,
#'     `EMPHDXYt`, `CHOLDXYt`, `DIABDXYt_M18`, `ARTHDXYt`, `ASTHDXYt`);
#'     refused, don't know and inapplicable answers count as not diagnosed.
#'     The longitudinal file carries the same "ever diagnosed" answers in
#'     both years, so `nchronic` is constant within person: a baseline
#'     morbidity measure, dropped by within (fixed effects) estimators}
#'   \item{longwt}{numeric MEPS panel 23 longitudinal person weight
#'     (`LONGWT`), constant within person}
#'   \item{varstr}{integer variance estimation stratum (`VARSTR`)}
#'   \item{varpsu}{integer variance estimation primary sampling unit within
#'     stratum (`VARPSU`)}
#' }
#'
#' @source Agency for Healthcare Research and Quality (AHRQ), Medical
#'   Expenditure Panel Survey, Household Component, public use file
#'   HC-217: MEPS Panel 23 Longitudinal Data File (released December
#'   2021). Rockville, MD: AHRQ.
#'   <https://meps.ahrq.gov/mepsweb/data_stats/download_data_files_detail.jsp?cboPufNumber=HC-217>.
#'   Data: <https://meps.ahrq.gov/mepsweb/data_files/pufs/h217/h217dta.zip>.
#'
#'   MEPS is produced by AHRQ, an agency of the US federal government; the
#'   public use files carry no copyright notice. AHRQ "requests that users
#'   cite AHRQ and the Medical Expenditure Panel Survey as the data source
#'   in any publications or research based upon these data." The files are
#'   released under the following Data Use Agreement (section A of the
#'   HC-217 documentation,
#'   <https://meps.ahrq.gov/mepsweb/data_stats/download_data/pufs/h217/h217doc.shtml>),
#'   which binds every user of this extract as well:
#'
#'   "Individual identifiers have been removed from the micro-data
#'   contained in these files. Nevertheless, under sections 308 (d) and 903
#'   (c) of the Public Health Service Act (42 U.S.C. 242m and 42 U.S.C. 299
#'   a-1), data collected by the Agency for Healthcare Research and Quality
#'   (AHRQ) and/or the National Center for Health Statistics (NCHS) may not
#'   be used for any purpose other than for the purpose for which they were
#'   supplied; any effort to determine the identity of any reported cases
#'   is prohibited by law. Therefore in accordance with the above
#'   referenced Federal Statute, it is understood that: (1) No one is to
#'   use the data in this data set in any way except for statistical
#'   reporting and analysis; and (2) If the identity of any person or
#'   establishment should be discovered inadvertently, then (a) no use will
#'   be made of this knowledge, (b) the Director Office of Management AHRQ
#'   will be advised of this incident, (c) the information that would
#'   identify any individual or establishment will be safeguarded or
#'   destroyed, as requested by AHRQ, and (d) no one else will be informed
#'   of the discovered identity; and (3) No one will attempt to link this
#'   data set with individually identifiable records from any data sets
#'   other than the Medical Expenditure Panel Survey or the National Health
#'   Interview Survey. By using these data you signify your agreement to
#'   comply with the above stated statutorily based requirements with the
#'   knowledge that deliberately making a false statement in any matter
#'   within the jurisdiction of any department or agency of the Federal
#'   Government violates Title 18 part 1 Chapter 47 Section 1001 and is
#'   punishable by a fine of up to $10,000 or up to 5 years in prison."
#'
#'   Validation (in `data-raw/meps.R`): on the full HC-217 file the
#'   unweighted codebook frequencies of `INSCOVY1`, `RACETHX`, `POVCATY2`,
#'   `TOTEXPY1` and `OBTOTVY2` are reproduced exactly, and `LONGWT` sums to
#'   305.7 million over the 13,044 `ALL5RDS = 1` persons as in the
#'   documentation. The `LONGWT`-weighted mean total expenditure per
#'   person is 5,676 dollars in 2018 and 6,249 dollars in 2019, against
#'   6,063 and 6,252 dollars from the annual files in AHRQ Statistical
#'   Briefs #533 and #540; the share with any expense is 87.1 and 84.4
#'   percent against about 87 and 86 percent. The 2018 gap (-6 percent) is
#'   expected: the longitudinal weight covers only people in scope for the
#'   whole two years, which excludes 2018 decedents, whose end-of-life
#'   spending is high.
#'
#' @references
#' Mitchell, E. M. (2021). Concentration of healthcare expenditures and
#' selected characteristics of persons with high expenses, U.S. civilian
#' noninstitutionalized population, 2018. Statistical Brief #533. Agency
#' for Healthcare Research and Quality, Rockville, MD.
#' <https://meps.ahrq.gov/data_files/publications/st533/stat533.shtml>
#'
#' Mitchell, E. M. (2022). Concentration of healthcare expenditures and
#' selected characteristics of persons with high expenses, U.S. civilian
#' noninstitutionalized population, 2019. Statistical Brief #540. Agency
#' for Healthcare Research and Quality, Rockville, MD.
#' <https://meps.ahrq.gov/data_files/publications/st540/stat540.shtml>
#'
#' Cohen, J. W., Cohen, S. B. and Banthin, J. S. (2009). The Medical
#' Expenditure Panel Survey: a national information resource to support
#' healthcare cost research and inform policy and practice. *Medical
#' Care*, 47(7 Suppl 1), S44-S50. \doi{10.1097/MLR.0b013e3181a23e3a}
#'
#' Deb, P., Norton, E. C. and Manning, W. G. (2017). *Health Econometrics
#' Using Stata*. College Station, TX: Stata Press. (Two-part, generalized
#' linear and count models for MEPS expenditure and use data.)
#'
#' Duan, N., Manning, W. G., Morris, C. N. and Newhouse, J. P. (1983). A
#' comparison of alternative models for the demand for medical care.
#' *Journal of Business and Economic Statistics*, 1, 115-126.
#' \doi{10.1080/07350015.1983.10509330}
#'
#' @section Teaching use:
#' All outcomes are observed for every person in both years, so the data
#' are a balanced short panel (T = 2): within and Mundlak (CRE) estimators
#' are identified only from people whose insurance or income changed
#' between 2018 and 2019 (223 insurance changes, 1,114 income-category
#' changes), and `age` and `nchronic` should enter only pooled or random
#' effects fits.
#'
#' * **Does insurance status change health spending, and through which
#'   margin?** `totexp` is zero for 14 percent of person-years and highly
#'   skewed: [twopart()] (probit + gamma, `model = "pooling"` or
#'   `"random"`) with [partialeffects()] for the dollar effect of being
#'   uninsured; compare with `panglm(family = "tobit")` on `log1p(totexp)`.
#' * **Is the cross-sectional insurance gap causal or selection?**
#'   [panglm()] `family = "gamma"` on positive spending, `model =
#'   "pooling"` or `"random"` versus `model = "within"` (identified from
#'   people who gained or lost coverage).
#' * **Who uses any care?** `anyexp` with `family = "probit"`, `model =
#'   "random"`, `cre = TRUE` (Mundlak test of the random effects
#'   assumption) and [partialeffects()].
#' * **Utilisation counts.** `officevis` (overdispersed, 24 percent
#'   zeros) or `ervis`, `ipstays` (mostly zero) with `family = "negbin"` or
#'   `"poisson"`, `model = "random"` or `"within"`; [panhurdle()] for the
#'   zero versus positive split with fixed effects.
#' * **Distributional effects.** [panqr()] (`model = "canay"` or
#'   `"pooling"`) on `log1p(totexp)` or `log(oopexp + 1)` at several
#'   quantiles; `model = "mmqr"` needs T >= 3 and does not apply here.
#' * **Selection on spending.** [panselect()] with `anyexp` as the
#'   selection equation and log spending as the outcome (the alternative
#'   to the two-part model); self-rated `health` as an ordered outcome
#'   with `family = "oprobit"`, `model = "random"`.
#'
#' @keywords datasets
#' @examples
#' data(meps)
#' str(meps)
#' # survey-weighted mean spending per adult and share with any expense
#' sapply(split(meps, meps$year), function(d)
#'   c(mean = weighted.mean(d$totexp, d$longwt),
#'     anyexp = weighted.mean(d$anyexp, d$longwt)))
#'
#' # two-part model: P(any spending) and E(spending | spending > 0)
#' tp <- twopart(totexp ~ insurance + I(age / 10) + sex + nchronic + povcat,
#'               data = meps, index = c("id", "year"), positive = "gamma")
#' summary(tp)
#' partialeffects(tp, variables = "insurance")   # dollars per person-year
#'
#' # positive spending: pooled versus within (coverage switchers)
#' pos <- subset(meps, totexp > 0)
#' fe <- panglm(totexp ~ insurance + povcat, data = pos,
#'              index = c("id", "year"), family = "gamma", model = "within")
#' summary(fe)
#'
#' # office-based visits: random effects negative binomial
#' nb <- panglm(officevis ~ insurance + I(age / 10) + sex + nchronic,
#'              data = meps, index = c("id", "year"), family = "negbin",
#'              model = "random")
#' summary(nb)
#'
#' \donttest{
#' # any use: correlated random effects probit with Mundlak test
#' cre <- panglm(anyexp ~ insurance + povcat + sex + nchronic, data = meps,
#'               index = c("id", "year"), family = "probit",
#'               model = "random", cre = TRUE)
#' summary(cre)
#'
#' # censored log spending: tobit
#' summary(panglm(log1p(totexp) ~ insurance + I(age / 10) + sex + nchronic,
#'                data = meps, index = c("id", "year"), family = "tobit",
#'                model = "pooling"))
#'
#' # quantiles of log spending with individual effects (the many tied
#' # values of log1p(totexp) make quantreg warn that the solution may be
#' # nonunique)
#' qr <- suppressWarnings(panqr(log1p(totexp) ~ insurance + povlev,
#'   data = meps, index = c("id", "year"), tau = c(0.25, 0.5, 0.75),
#'   model = "canay", boot = 50, seed = 1))
#' summary(qr)
#'
#' # selection correction for log spending observed only if positive
#' meps$logexp <- ifelse(meps$totexp > 0, log(meps$totexp), NA)
#' summary(panselect(anyexp ~ insurance + povcat, logexp ~ insurance,
#'                   data = meps, index = c("id", "year"),
#'                   model = "within", boot = 49, seed = 1))
#'
#' # self-rated health: random effects ordered probit
#' summary(panglm(health ~ nchronic + I(age / 10) + povcat, data = meps,
#'                index = c("id", "year"), family = "oprobit",
#'                model = "random"))
#' }
"meps"
