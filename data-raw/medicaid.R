# State panel of health insurance coverage and the ACA Medicaid expansion,
# 50 states + DC, 2009-2019, built from the American Community Survey (ACS)
# 1-year Summary File of the US Census Bureau.
#
# ---------------------------------------------------------------------------
# Sources and licence
# ---------------------------------------------------------------------------
# * ACS 1-year Summary File (detailed tables B27001 and B27016), state
#   summary level (logical record 1 of each state file):
#     https://www2.census.gov/programs-surveys/acs/summary_file/<year>/
#   Tables are located through the Census sequence/table lookup files
#   (documentation/.../user_tools/). The Census Data API
#   (api.census.gov) now refuses requests without a key (redirect to
#   missing_key.html, checked 2026-11-01), so the script reads the
#   Summary File, which holds the same published estimates and needs no
#   key.
# * Licence: works of the US federal government are not subject to
#   copyright (17 U.S.C. 105). The Census Bureau states, at
#   https://www.census.gov/help/topics/faq.how-do-i-cite-or-source-acs-data.html
#   (checked 2026-11-01): "Our data tables, visualizations, reports and
#   other products are freely available for public use." and "You do not
#   need permission to use American Community Survey data." Attribution:
#   "Source: U.S. Census Bureau, American Community Survey 1-year
#   estimates, 2009-2019, tables B27001 and B27016."
# * Medicaid expansion implementation dates: KFF, "Status of State Action
#   on the Medicaid Expansion Decision" (State Health Facts),
#   https://www.kff.org/affordable-care-act/state-indicator/state-activity-around-expanding-medicaid-under-the-affordable-care-act/
#   (status as of October 2026, read 2026-11-01). Dates are facts and are
#   re-typed below, not copied as a table.
# * Expansion-year convention: Miller, Johnson and Wherry (2021, QJE
#   136(3):1783-1829, doi:10.1093/qje/qjab004; NBER WP 26081, footnote 14):
#   a state is coded as expanding in the calendar year of implementation if
#   the expansion took effect before July 1 (Michigan, April 2014 -> 2014;
#   Indiana, February 2015 -> 2015), otherwise in the following year (New
#   Hampshire, August 2014 -> 2015; Alaska, September 2015 -> 2016;
#   Louisiana, July 2016 -> 2017). The same rule is applied to the later
#   adopters.
# * Pre-ACA coverage of low-income adults: Miller, Johnson and Wherry
#   (2021, footnote 6): "DE, MA, NY, and VT all expanded coverage to
#   individuals with incomes reaching the poverty line or greater prior to
#   the ACA; DC received approval to implement its ACA Medicaid expansion
#   early with enrollment starting in 2011" (coded prior = "full"). Early
#   limited expansions under the ACA option or waivers in CA, CT, MN, NJ,
#   WA (and DC): Sommers, Arntson, Kenney and Epstein (2013), Medicare and
#   Medicaid Research Review 3(4), doi:10.5600/mmrr.003.04.a02 (coded
#   prior = "early").
#
# ---------------------------------------------------------------------------
# Design choices
# ---------------------------------------------------------------------------
# * Years 2009-2019. 2008 (first ACS year with the health insurance
#   question) is dropped: its Summary File has only the 10-cell B27001 (no
#   sex detail, no table by ratio of income to poverty, no Medicaid table),
#   so the adult and low-income outcomes cannot be built. 2020 is dropped
#   because the Census Bureau did not release standard 2020 ACS 1-year
#   estimates (only experimental ones) due to pandemic data-collection
#   problems; the panel therefore ends in 2019.
# * Adults: the ACS tables used "18 to 64 years" through 2016 and
#   "19 to 64 years" from 2017 on (2017 table redesign). The adult measures
#   follow the published age band of each year; the break is common to all
#   states and is absorbed by year effects.
# * Low income: below 100 percent of poverty is available every year
#   (uninspov). The ACA cut (below 138 percent of poverty, bins under 0.50,
#   0.50 to 0.99 and 1.00 to 1.37) exists only from 2011 on (2009-2010
#   tables cut at 1.00 and 1.50), so uninslow is NA in 2009-2010.
# * No sampling: all 51 jurisdictions, all years (balanced, 561 rows).

raw_dir <- "/home/imad-el-badisy/Desktop/BIOSTAT/PROJECTS/ACTIVE/software-panglm/archive/raw/datasets-2026-11-01/medicaid"
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
base_url <- "https://www2.census.gov/programs-surveys/acs/summary_file"
years <- 2009:2019

fetch <- function(url, dest) {
  if (!file.exists(dest)) {
    utils::download.file(url, dest, mode = "wb", quiet = TRUE)
  }
  dest
}

# ---------------------------------------------------------------------------
# 1. Sequence/table lookup files
# ---------------------------------------------------------------------------
lookup_url <- function(y) {
  if (y == 2009) {
    file.path(base_url, y, "documentation/1_year/user_tools/merge_5_6.txt")
  } else if (y <= 2012) {
    file.path(base_url, y,
              "documentation/1_year/user_tools/Sequence_Number_and_Table_Number_Lookup.txt")
  } else {
    file.path(base_url, y,
              "documentation/user_tools/ACS_1yr_Seq_Table_Number_Lookup.txt")
  }
}

read_lookup <- function(y) {
  f <- fetch(lookup_url(y), file.path(raw_dir, paste0("lookup", y, ".txt")))
  lk <- utils::read.csv(f, colClasses = "character", check.names = FALSE,
                        strip.white = TRUE, fileEncoding = "latin1")
  names(lk) <- c("file", "table", "seq", "line", "start", "cells",
                 "seqcells", "title", "subject")[seq_len(ncol(lk))]
  lk
}

# position (column of the e-file) and label of every line of a table
table_lines <- function(lk, tab) {
  x <- lk[lk$table == tab, ]
  head <- x[grepl("CELL", x$cells), ]
  stopifnot(nrow(head) == 1)
  start <- as.integer(head$start)
  seq <- as.integer(head$seq)
  ln <- suppressWarnings(as.numeric(x$line))
  x <- x[!is.na(ln) & ln == round(ln), ]
  data.frame(seq = seq, line = as.integer(x$line),
             col = start + as.integer(x$line) - 1L,
             label = x$title, stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
# 2. State list (abbreviation, FIPS, Census region)
# ---------------------------------------------------------------------------
states <- data.frame(
  name = c(datasets::state.name, "District of Columbia"),
  abb = c(datasets::state.abb, "DC"),
  region = c(as.character(datasets::state.region), "South"),
  stringsAsFactors = FALSE
)
fips <- c(AL = 1, AK = 2, AZ = 4, AR = 5, CA = 6, CO = 8, CT = 9, DE = 10,
          DC = 11, FL = 12, GA = 13, HI = 15, ID = 16, IL = 17, IN = 18,
          IA = 19, KS = 20, KY = 21, LA = 22, ME = 23, MD = 24, MA = 25,
          MI = 26, MN = 27, MS = 28, MO = 29, MT = 30, NE = 31, NV = 32,
          NH = 33, NJ = 34, NM = 35, NY = 36, NC = 37, ND = 38, OH = 39,
          OK = 40, OR = 41, PA = 42, RI = 44, SC = 45, SD = 46, TN = 47,
          TX = 48, UT = 49, VT = 50, VA = 51, WA = 53, WV = 54, WI = 55,
          WY = 56)
states$fips <- as.integer(fips[states$abb])
# Census regions (base R uses "North Central" for the Midwest)
states$region[states$region == "North Central"] <- "Midwest"
states <- states[order(states$name), ]

state_dir <- function(y, name, abb) {
  if (y == 2009) abb else gsub(" ", "", name)
}

# value cells of logical record 1 (the state itself) of one sequence file
read_seq <- function(y, name, abb, seq) {
  st <- tolower(abb)
  stem <- sprintf("%d1%s%04d000", y, st, seq)
  zip <- fetch(file.path(base_url, y, "data/1_year_seq_by_state",
                         state_dir(y, name, abb), paste0(stem, ".zip")),
               file.path(raw_dir, paste0(stem, ".zip")))
  con <- unz(zip, paste0("e", stem, ".txt"))
  on.exit(close(con))
  first <- readLines(con, n = 1)
  v <- strsplit(first, ",", fixed = TRUE)[[1]]
  stopifnot(v[6] == "0000001", toupper(v[3]) == toupper(abb))
  v
}

# ---------------------------------------------------------------------------
# 3. Build one row per state and year
# ---------------------------------------------------------------------------
build_year <- function(y) {
  lk <- read_lookup(y)
  b01 <- table_lines(lk, "B27001")
  b16 <- table_lines(lk, "B27016")
  # B27001: total and every "No health insurance coverage" line
  b01_tot <- b01$col[b01$line == 1]
  b01_no <- b01$col[grepl("^No health insurance", b01$label)]
  # B27016: walk the lines, tracking the poverty bin and the age group
  bin <- age <- NA_character_
  b16$bin <- b16$age <- NA_character_
  for (i in seq_len(nrow(b16))) {
    lab <- b16$label[i]
    if (grepl("poverty threshold", lab)) {
      bin <- lab; age <- NA_character_
    } else if (grepl("years", lab)) {
      age <- lab
    }
    b16$bin[i] <- bin
    b16$age[i] <- age
  }
  adult_lab <- if (y <= 2016) "18 to 64 years:" else "19 to 64 years:"
  stopifnot(sum(b16$label == adult_lab) %in% c(7, 9))
  is_adult_total <- b16$label == adult_lab
  is_adult <- !is.na(b16$age) & b16$age == adult_lab
  is_no <- is_adult & grepl("^No health insurance", b16$label)
  is_mcd <- is_adult & grepl("Medicaid", b16$label)
  below100 <- grepl("^(Under 0.50|0.50 to .99)", b16$bin)
  below138 <- grepl("^(Under 0.50|0.50 to .99|1.00 to 1.37)", b16$bin)
  has138 <- any(grepl("^1.00 to 1.37", b16$bin))

  rows <- lapply(seq_len(nrow(states)), function(s) {
    nm <- states$name[s]; ab <- states$abb[s]
    e01 <- read_seq(y, nm, ab, b01$seq[1])
    e16 <- if (b16$seq[1] == b01$seq[1]) e01 else read_seq(y, nm, ab, b16$seq[1])
    g01 <- function(cols) sum(as.numeric(e01[cols]))
    g16 <- function(sel) sum(as.numeric(e16[b16$col[sel]]))
    data.frame(
      abb = ab, year = y,
      pop = g01(b01_tot), unins = g01(b01_no),
      adultpop = g16(is_adult_total), adultunins = g16(is_no),
      adultmcd = g16(is_mcd),
      povpop = g16(is_adult_total & below100),
      povunins = g16(is_no & below100), povmcd = g16(is_mcd & below100),
      lowpop = if (has138) g16(is_adult_total & below138) else NA_real_,
      lowunins = if (has138) g16(is_no & below138) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

counts <- do.call(rbind, lapply(years, build_year))

# ---------------------------------------------------------------------------
# 4. Medicaid expansion dates (KFF, status as of October 2026)
# ---------------------------------------------------------------------------
expdates <- c(
  AK = "2015-09-01", AZ = "2014-01-01", AR = "2014-01-01", CA = "2014-01-01",
  CO = "2014-01-01", CT = "2014-01-01", DE = "2014-01-01", DC = "2014-01-01",
  HI = "2014-01-01", ID = "2020-01-01", IL = "2014-01-01", IN = "2015-02-01",
  IA = "2014-01-01", KY = "2014-01-01", LA = "2016-07-01", ME = "2019-01-10",
  MD = "2014-01-01", MA = "2014-01-01", MI = "2014-04-01", MN = "2014-01-01",
  MO = "2021-10-01", MT = "2016-01-01", NE = "2020-10-01", NV = "2014-01-01",
  NH = "2014-08-15", NJ = "2014-01-01", NM = "2014-01-01", NY = "2014-01-01",
  NC = "2023-12-01", ND = "2014-01-01", OH = "2014-01-01", OK = "2021-07-01",
  OR = "2014-01-01", PA = "2015-01-01", RI = "2014-01-01", SD = "2023-07-01",
  UT = "2020-01-01", VT = "2014-01-01", VA = "2019-01-01", WA = "2014-01-01",
  WV = "2014-01-01"
)
nonadopters <- c("AL", "FL", "GA", "KS", "MS", "SC", "TN", "TX", "WI", "WY")
stopifnot(length(expdates) == 41, length(nonadopters) == 10,
          setequal(c(names(expdates), nonadopters), states$abb))
expdate <- as.Date(expdates[states$abb])
mo <- as.integer(format(expdate, "%m"))
yr <- as.integer(format(expdate, "%Y"))
states$expdate <- expdate
states$adopt <- ifelse(mo < 7, yr, yr + 1L)
states$prior <- ifelse(states$abb %in% c("DE", "DC", "MA", "NY", "VT"), "full",
                ifelse(states$abb %in% c("CA", "CT", "MN", "NJ", "WA"),
                       "early", "none"))

# ---------------------------------------------------------------------------
# 5. Assemble
# ---------------------------------------------------------------------------
d <- merge(counts, states, by = "abb")
pct <- function(num, den) round(100 * num / den, 2)
medicaid <- data.frame(
  state = factor(d$name, levels = states$name),
  year = as.integer(d$year),
  abb = factor(d$abb, levels = sort(states$abb)),
  fips = d$fips,
  region = factor(d$region, levels = c("Northeast", "Midwest", "South", "West")),
  pop = as.integer(d$pop),
  uninsured = pct(d$unins, d$pop),
  adultpop = as.integer(d$adultpop),
  uninsadult = pct(d$adultunins, d$adultpop),
  medicaidadult = pct(d$adultmcd, d$adultpop),
  povpop = as.integer(d$povpop),
  povbase = NA_integer_,
  uninspov = pct(d$povunins, d$povpop),
  medicaidpov = pct(d$povmcd, d$povpop),
  lowpop = as.integer(d$lowpop),
  uninslow = pct(d$lowunins, d$lowpop),
  adopt = as.integer(d$adopt),
  expdate = d$expdate,
  expanded = as.integer(!is.na(d$adopt) & d$year >= d$adopt),
  prior = factor(d$prior, levels = c("none", "early", "full")),
  stringsAsFactors = FALSE
)
medicaid <- medicaid[order(medicaid$state, medicaid$year), ]
rownames(medicaid) <- NULL
# time-invariant weight for attgt(), which needs weights constant within
# states: the 2013 (last pre-expansion year) count of adults below poverty
base13 <- medicaid$povpop[medicaid$year == 2013]
names(base13) <- as.character(medicaid$state[medicaid$year == 2013])
medicaid$povbase <- as.integer(base13[as.character(medicaid$state)])

# ---------------------------------------------------------------------------
# 6. Validation
# ---------------------------------------------------------------------------
# (a) Percent uninsured, all ages, against the PUBLISHED ACS ranking table
#     R2701 "Percent of People Who Have No Health Insurance Coverage"
#     (1-year, states), 2018 and 2019, all 51 jurisdictions:
#     https://www2.census.gov/programs-surveys/acs/summary_file/<year>/data/1_year_ranking/R2701.xlsx
#     Published to one decimal, so the tolerance is 0.05 + rounding slack.
if (requireNamespace("readxl", quietly = TRUE)) {
  rk <- do.call(rbind, lapply(2018:2019, function(y) {
    f <- fetch(file.path(base_url, y, "data/1_year_ranking/R2701.xlsx"),
               file.path(raw_dir, paste0("R2701_", y, ".xlsx")))
    x <- as.data.frame(readxl::read_excel(f, col_names = FALSE, .name_repair = "minimal"))
    x <- x[as.character(x[[2]]) %in% states$name, ]
    data.frame(state = as.character(x[[2]]), year = y,
               published = as.numeric(x[[3]]))
  }))
  chk <- merge(rk, data.frame(state = as.character(medicaid$state),
                              year = medicaid$year,
                              computed = medicaid$uninsured))
  chk$diff <- chk$computed - chk$published
  print(chk[order(-abs(chk$diff)), ][1:8, ], row.names = FALSE)
  stopifnot(nrow(chk) == 102, max(abs(chk$diff)) <= 0.055)
}

# (b) Hard-coded published values from the Census report: Smith, J. C. and
#     Medalia, C. (2015), Health Insurance Coverage in the United States:
#     2014, Current Population Reports P60-253, U.S. Census Bureau,
#     Appendix Table A-1 "Population Without Health Insurance Coverage by
#     State: 2013 and 2014" (ACS; percent uninsured, all ages), p. 26,
#     https://www.census.gov/content/dam/Census/library/publications/2015/demo/p60-253.pdf
#     Twelve jurisdictions covering every adoption cohort and the prior
#     coverage states; the pre/post 2014 drop is the core of the design.
published <- data.frame(
  abb = rep(c("AR", "CA", "DC", "FL", "KY", "LA", "MA", "MI", "NV", "PA",
              "TX", "WI"), times = 2),
  year = rep(2013:2014, each = 12),
  value = c(16.0, 17.2, 6.7, 20.0, 14.3, 16.6, 3.7, 11.0, 20.7, 9.7, 22.1, 9.1,
            11.8, 12.4, 5.3, 16.6, 8.5, 14.8, 3.3, 8.5, 15.2, 8.5, 19.1, 7.3),
  stringsAsFactors = FALSE
)
v <- merge(published, data.frame(abb = as.character(medicaid$abb),
                                 year = medicaid$year,
                                 computed = medicaid$uninsured))
v$diff <- round(v$computed - v$value, 3)
print(v, row.names = FALSE)
stopifnot(nrow(v) == 24, max(abs(v$diff)) <= 0.055)

# (c) Low-income adults: the B27016 aggregation (sum of the detailed
#     poverty bins) against the separately published collapsed table
#     C27016 (lines "Under 1.00" and "1.00 to 1.37 of poverty threshold",
#     19 to 64 years, no health insurance coverage), 2019, all states.
lk19 <- read_lookup(2019)
c16 <- table_lines(lk19, "C27016")
bin <- age <- NA_character_
for (i in seq_len(nrow(c16))) {
  lab <- c16$label[i]
  if (grepl("poverty threshold", lab)) bin <- lab else if (grepl("years", lab)) age <- lab
  c16$bin[i] <- bin; c16$age[i] <- age
}
c_tot <- c16$label == "19 to 64 years:" & grepl("^(Under 1.00|1.00 to 1.37)", c16$bin)
c_no <- grepl("^No health", c16$label) & c16$age %in% "19 to 64 years:" &
  grepl("^(Under 1.00|1.00 to 1.37)", c16$bin)
stopifnot(sum(c_tot) == 2, sum(c_no) == 2)
c16chk <- do.call(rbind, lapply(seq_len(nrow(states)), function(s) {
  e <- read_seq(2019, states$name[s], states$abb[s], c16$seq[1])
  data.frame(abb = states$abb[s],
             c27016 = round(100 * sum(as.numeric(e[c16$col[c_no]])) /
                              sum(as.numeric(e[c16$col[c_tot]])), 2))
}))
c16chk <- merge(c16chk, data.frame(abb = as.character(medicaid$abb),
                                   year = medicaid$year,
                                   uninslow = medicaid$uninslow)[medicaid$year == 2019, ])
print(head(c16chk[order(c16chk$abb), ], 6), row.names = FALSE)
stopifnot(nrow(c16chk) == 51, max(abs(c16chk$c27016 - c16chk$uninslow)) < 0.011)

# (d) Expansion coding against KFF counts and the Miller, Johnson and
#     Wherry (2021) cohorts (footnote 14).
st <- unique(medicaid[, c("abb", "adopt", "prior")])
stopifnot(
  sum(!is.na(st$adopt)) == 41,                                 # KFF: 41 incl. DC
  setequal(as.character(st$abb[st$adopt %in% 2014]),
           c("AR", "AZ", "CA", "CO", "CT", "HI", "IL", "IA", "KY", "MD",
             "MI", "MN", "NJ", "NM", "NV", "ND", "OH", "OR", "RI", "WA",
             "WV", "DE", "DC", "MA", "NY", "VT")),               # MJW 2014 + prior states
  setequal(as.character(st$abb[st$adopt %in% 2015]), c("NH", "PA", "IN")),
  setequal(as.character(st$abb[st$adopt %in% 2016]), c("AK", "MT")),
  setequal(as.character(st$abb[st$adopt %in% 2017]), "LA"),
  setequal(as.character(st$abb[st$adopt %in% 2019]), c("ME", "VA")),
  setequal(as.character(st$abb[st$adopt %in% 2020]), c("ID", "UT")),
  setequal(as.character(st$abb[st$adopt %in% 2021]), "NE"),
  setequal(as.character(st$abb[st$adopt %in% 2022]), c("OK", "MO")),
  setequal(as.character(st$abb[st$adopt %in% 2024]), c("SD", "NC"))
)

# ---------------------------------------------------------------------------
# 7. Save
# ---------------------------------------------------------------------------
stopifnot(nrow(medicaid) == 51 * length(years),
          !anyDuplicated(medicaid[, c("state", "year")]))
save(medicaid, file = "data/medicaid.rda", compress = "xz")
