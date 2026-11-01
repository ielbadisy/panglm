# Build the `healthspend` dataset: a country-year panel of health
# expenditure, income and population health, 2000-2022, from the World
# Bank World Development Indicators (WDI).
#
# Source
#   World Bank, World Development Indicators, API v2
#   https://api.worldbank.org/v2/country/all/indicator/<code>?date=2000:2023&format=json&per_page=20000
#   (years then subset to 2000-2022)
#   Country metadata (region, income group):
#   https://api.worldbank.org/v2/country?format=json&per_page=400
#   Vintage used for the shipped data: WDI "lastupdated" 2026-10-08,
#   health expenditure series "Global Health Expenditure Database, updated
#   December 12th, 2025, World Health Organization (WHO)". Downloaded
#   2026-11-01.
#
# Licence (verified 2026-11-01)
#   CC BY 4.0. The WDI catalog page
#   https://datacatalog.worldbank.org/search/dataset/0037712/World-Development-Indicators
#   states: "This dataset is licensed under Creative Commons Attribution
#   4.0." The per-indicator WDI metadata
#   (https://api.worldbank.org/v2/sources/2/indicators/<code>/metadata?format=json)
#   carries "License_Type": "CC BY-4.0" for every one of the nine series
#   used here, including the four WHO GHED health expenditure series.
#   The World Bank licence page https://datacatalog.worldbank.org/public-licenses
#   says CC BY 4.0 "allows users to copy, modify and distribute data in any
#   format for any purpose, including commercial use", with attribution
#   and an indication of changes. Changes made here: country subset,
#   renaming, rounding, one derived column (chereal).
#   Note: WHO's own dataset terms
#   (https://www.who.int/about/policies/publishing/data-policy/terms-and-conditions)
#   are more restrictive than CC BY 4.0; this build uses the copies
#   redistributed by the World Bank under CC BY 4.0.
#
# Raw files are cached outside the package (never committed) in `rawdir`.
# Downloads are treated as data only (parsed with jsonlite, nothing run).

library(jsonlite)

rawdir <- Sys.getenv(
  "HEALTHSPEND_RAW",
  "../archive/raw/datasets-2026-11-01/healthspend"
)
dir.create(rawdir, recursive = TRUE, showWarnings = FALSE)

years <- 2000:2022
codes <- c(
  chepc   = "SH.XPD.CHEX.PP.CD",  # CHE per capita, PPP, current intl $
  chegdp  = "SH.XPD.CHEX.GD.ZS",  # CHE % of GDP
  oop     = "SH.XPD.OOPC.CH.ZS",  # out-of-pocket % of CHE
  gghe    = "SH.XPD.GHED.CH.ZS",  # domestic general government % of CHE
  gdppc   = "NY.GDP.PCAP.PP.KD",  # GDP per capita, PPP, constant 2021 intl $
  lifeexp = "SP.DYN.LE00.IN",     # life expectancy at birth, years
  u5mr    = "SH.DYN.MORT",        # under-5 mortality per 1,000 live births
  pop     = "SP.POP.TOTL",        # population, total
  pop65   = "SP.POP.65UP.TO.ZS"   # population aged 65+, % of total
)

getraw <- function(url, file) {
  path <- file.path(rawdir, file)
  if (!file.exists(path)) {
    utils::download.file(url, path, mode = "wb", quiet = TRUE)
  }
  path
}

## ---- country metadata ------------------------------------------------------
cm <- fromJSON(getraw(
  "https://api.worldbank.org/v2/country?format=json&per_page=400",
  "countries.json"
))[[2]]
cm <- data.frame(
  iso3c   = cm$id,
  country = trimws(cm$name),
  region  = trimws(cm$region$value),
  income  = trimws(cm$incomeLevel$value),
  stringsAsFactors = FALSE
)
cm <- cm[cm$region != "Aggregates", ]          # drop regional/income aggregates

## ---- indicators ------------------------------------------------------------
readind <- function(code) {
  url <- sprintf(
    "https://api.worldbank.org/v2/country/all/indicator/%s?date=2000:2023&format=json&per_page=20000",
    code
  )
  js <- fromJSON(getraw(url, paste0(code, ".json")))
  stopifnot(js[[1]]$pages == 1)                 # everything in one page
  d <- js[[2]]
  data.frame(iso3c = d$countryiso3code, year = as.integer(d$date),
             value = as.numeric(d$value), stringsAsFactors = FALSE)
}

panel <- expand.grid(year = years, iso3c = sort(cm$iso3c),
                     stringsAsFactors = FALSE)[, c("iso3c", "year")]
for (v in names(codes)) {
  d <- readind(codes[[v]])
  d <- d[d$iso3c %in% cm$iso3c & d$year %in% years, ]
  stopifnot(!anyDuplicated(d[, c("iso3c", "year")]))
  names(d)[3] <- v
  panel <- merge(panel, d, by = c("iso3c", "year"), all.x = TRUE)
}

## keep countries with at least one year of health expenditure data
has <- tapply(!is.na(panel$chepc), panel$iso3c, any)
panel <- panel[panel$iso3c %in% names(has)[has], ]
panel <- merge(panel, cm, by = "iso3c", all.x = TRUE)

## ---- types, rounding, derived column ---------------------------------------
inclev <- c("Low income", "Lower middle income", "Upper middle income",
            "High income")
panel$income[!panel$income %in% inclev] <- NA   # "Not classified"
panel$income <- factor(panel$income, levels = inclev,
                       labels = c("low", "lowermiddle", "uppermiddle", "high"))
panel$region <- factor(panel$region)
panel$iso3c  <- factor(panel$iso3c)
panel$year   <- as.integer(panel$year)
panel$pop    <- as.integer(round(panel$pop))

panel$chepc   <- round(panel$chepc, 2)
panel$gdppc   <- round(panel$gdppc, 1)
for (v in c("chegdp", "oop", "gghe", "pop65")) panel[[v]] <- round(panel[[v]], 3)
panel$lifeexp <- round(panel$lifeexp, 3)
panel$u5mr    <- round(panel$u5mr, 1)
## real CHE per capita, constant 2021 PPP dollars: CHE share of GDP times
## real GDP per capita (both from WDI)
panel$chereal <- round(panel$chegdp / 100 * panel$gdppc, 2)

healthspend <- panel[order(panel$iso3c, panel$year),
  c("iso3c", "year", "country", "region", "income",
    "chepc", "chereal", "chegdp", "oop", "gghe",
    "gdppc", "lifeexp", "u5mr", "pop", "pop65")]
rownames(healthspend) <- NULL

## ---- validation ------------------------------------------------------------
## (1) Independent channel: values read on 2026-11-01 from the WDI bulk CSV
## downloads (https://api.worldbank.org/v2/en/indicator/<code>?downloadformat=csv,
## files API_<code>_DS2_en_csv_v2_*.csv, "Last Updated Date" 2026-10-08),
## hard-coded here. They come from a different distribution channel (the
## DataBank bulk files) than the JSON API used above. A later rebuild
## against a newer WDI release may move these by revision; then refresh the
## reference values from the bulk files of that release.
ref <- read.table(header = TRUE, text = "
iso3c year var     ref
USA   2000 chepc   4548.55482519
USA   2019 chepc  10550.22069455
USA   2021 chepc  12080.07362010
DEU   2010 chepc   4301.94597122
IND   2019 chepc    211.16661580
MAR   2000 chepc    143.24442368
MAR   2021 chepc    477.34754814
NGA   2010 chepc    149.82187001
BRA   2019 chepc   1545.03321145
USA   2019 oop       11.30872917
IND   2000 oop       71.70202637
MAR   2010 oop       55.61273193
NGA   2021 oop       76.24337006
DEU   2021 oop       11.36333561
BRA   2010 lifeexp   73.779
NGA   2000 lifeexp   47.143
MAR   2019 lifeexp   74.245
DEU   2019 lifeexp   81.29268293
")
ref$ours <- mapply(function(i, y, v)
  healthspend[[v]][healthspend$iso3c == i & healthspend$year == y],
  ref$iso3c, ref$year, ref$var)
ref$reldiff <- signif(ref$ours / ref$ref - 1, 3)
print(ref, row.names = FALSE)
stopifnot(all(abs(ref$ours / ref$ref - 1) < 1e-3))   # rounding only

## (2) Internal accounting identity: OOP and domestic government shares
## are parts of CHE, so their sum cannot exceed 100 percent (up to rounding).
stopifnot(all(healthspend$oop + healthspend$gghe <= 100.01, na.rm = TRUE))

## (3) Structure.
stopifnot(
  !anyDuplicated(healthspend[, c("iso3c", "year")]),
  nlevels(healthspend$region) == 7,
  all(table(healthspend$iso3c) == length(years))
)
cat("countries:", nlevels(healthspend$iso3c), " rows:", nrow(healthspend), "\n")
print(table(healthspend$income, useNA = "ifany") / length(years))

save(healthspend, file = "data/healthspend.rda", compress = "xz")
