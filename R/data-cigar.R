#' Cigarette demand in US states, 1963-1992 (Baltagi and Levin)
#'
#' The classic real public panel of cigarette consumption, prices and
#' income for 46 US states (45 states plus the District of Columbia;
#' Alaska, Colorado, Hawaii, North Carolina and Oregon are not included)
#' observed every year from 1963 to 1992 (30 years), used by Baltagi and
#' Levin (1986, 1992) and in Baltagi's textbook to estimate a dynamic
#' demand equation for an addictive good. The panel is balanced with no
#' missing values. The first nine columns are the original file; the last
#' four are the log variables of the textbook equation
#' ln C(it) = a + b1 ln C(i,t-1) + b2 ln P(it) + b3 ln Pn(it) + b4 ln Y(it) + u(it).
#' The build script, with the validation against the published estimates,
#' is `data-raw/cigar.R` in the package's GitHub repository.
#'
#' @format A data frame with 1,380 rows (46 states x 30 years) and 13
#'   columns. The panel index pair is (`state`, `year`).
#' \describe{
#'   \item{state}{factor, two-letter US state abbreviation (46 levels,
#'     `DC` for the District of Columbia); the panel individual index}
#'   \item{year}{integer calendar year, 1963-1992; the panel time index}
#'   \item{price}{numeric, average retail price per pack of cigarettes,
#'     cents, nominal (taxes included)}
#'   \item{pop}{numeric, state population, thousands}
#'   \item{pop16}{numeric, population of smoking age, thousands (labelled
#'     "population above the age of 16" in the data archive; the textbook
#'     describes smoking age as 14 and older)}
#'   \item{cpi}{numeric, consumer price index, 1983 = 100}
#'   \item{ndi}{numeric, per capita disposable income, dollars, nominal}
#'   \item{sales}{numeric, cigarette sales, packs per head of total
#'     population}
#'   \item{pimin}{numeric, minimum price per pack in adjoining states,
#'     cents, nominal (proxy for casual smuggling across state borders)}
#'   \item{lnc}{numeric, log real per capita consumption by persons of
#'     smoking age, `log(sales * pop / pop16)`}
#'   \item{lnp}{numeric, log real price, `log(price / cpi * 100)` (1983
#'     cents)}
#'   \item{lnpn}{numeric, log real minimum neighbouring-state price,
#'     `log(pimin / cpi * 100)`}
#'   \item{lny}{numeric, log real per capita disposable income,
#'     `log(ndi / cpi * 100)` (1983 dollars)}
#' }
#'
#' @source File `Cigar.txt` of the data archive of Baltagi, Econometric
#'   Analysis of Panel Data (Wiley companion site),
#'   <https://www.wiley.com/legacy/wileychi/baltagi/supp/Cigar.txt>, listed
#'   at <https://www.wiley.com/legacy/wileychi/baltagi/datasets.html>,
#'   downloaded 2026-11-01. The original variables are identical to the
#'   copies distributed under GPL (>= 2) in the CRAN packages \pkg{plm}
#'   (`plm::Cigar`) and \pkg{Ecdat} (`Ecdat::Cigar`), which trace the data
#'   to the same archive; this package redistributes it on the same terms.
#'   The underlying series are public statistics (Tobacco Institute, The
#'   Tax Burden on Tobacco; US Census population estimates; consumer price
#'   index and personal income series). Cite Baltagi and Levin (1992) and
#'   Baltagi (2005) when using the data.
#'
#' @references
#' Baltagi, B. H. and Levin, D. (1986). Estimating dynamic demand for
#' cigarettes using panel data: the effects of bootlegging, taxation and
#' advertising reconsidered. Review of Economics and Statistics, 68(1),
#' 148-155.
#'
#' Baltagi, B. H. and Levin, D. (1992). Cigarette taxation: raising
#' revenues and reducing consumption. Structural Change and Economic
#' Dynamics, 3(2), 321-335.
#'
#' Baltagi, B. H., Griffin, J. M. and Xiong, W. (2000). To pool or not to
#' pool: homogeneous versus heterogeneous estimators applied to cigarette
#' demand. Review of Economics and Statistics, 82(1), 117-126.
#'
#' Baltagi, B. H. (2005). Econometric Analysis of Panel Data, 3rd edition.
#' Wiley. Section 8.8, Table 8.1 (Pooled estimation results, cigarette
#' demand equation 1963-92) and Table 8.2 (Arellano-Bond estimates).
#'
#' Arellano, M. and Bond, S. (1991). Some tests of specification for panel
#' data: Monte Carlo evidence and an application to employment equations.
#' Review of Economic Studies, 58(2), 277-297.
#'
#' Nickell, S. (1981). Biases in dynamic models with fixed effects.
#' Econometrica, 49(6), 1417-1426.
#'
#' @section Teaching use:
#' \itemize{
#'   \item How habit forming is smoking, and what are the short- and
#'     long-run price elasticities? Dynamic demand `lnc ~ lag(lnc, 1) + lnp
#'     + lnpn + lny`: pooled OLS (`panglm(model = "pooling")`) overstates
#'     persistence because the state effect sits in the lagged outcome;
#'     two-way within (`model = "within", effect = "twoways"`) understates
#'     it (Nickell bias, of order 1 / T); Arellano-Bond difference GMM and
#'     Blundell-Bond system GMM (`pangmm()` with `transformation = "d"` or
#'     `"ld"`) are consistent for fixed T. A credible GMM estimate of the
#'     lag coefficient lies between the within (0.83) and pooled (0.97)
#'     values. The long-run elasticity is b2 / (1 - b1).
#'   \item How much do the GMM results depend on the instrument set? With
#'     N = 46 states the full Arellano-Bond instrument matrix has 437
#'     columns; compare `collapse = TRUE` and `maxlag`, and read the
#'     Hansen test (it is uninformative, p close to 1, when instruments
#'     proliferate) and the AR(2) test.
#'   \item Does cross-border bootlegging matter? The coefficient of the
#'     neighbouring-state price `lnpn` in the static and dynamic models.
#'   \item Static demand: fixed versus random state effects for `lnc ~ lnp
#'     + lnpn + lny`, `panglm()` with `model = "within"` and `"random"`,
#'     and `panhausman()`. Here the classical contrast V_fe - V_re is not
#'     positive definite (the classical statistic is negative), so use the
#'     regression-based `method = "mundlak"` version, which rejects random
#'     effects (F = 7.17 on 3 and 45 df).
#'   \item Initial conditions: the Wooldridge dynamic correlated random
#'     effects model (`panglm(dynamic = TRUE)`), here a contrast with GMM
#'     because the state-effect variance is estimated at zero once the
#'     initial outcome is conditioned on.
#' }
#'
#' @examples
#' data(cigar)
#' ix <- c("state", "year")
#' f <- lnc ~ lag(lnc, 1) + lnp + lnpn + lny
#'
#' ## Baltagi (2005), Table 8.1: OLS, two-way within and one-step
#' ## Arellano-Bond GMM
#' ols <- panglm(f, data = cigar, index = ix, model = "pooling")
#' fe  <- panglm(f, data = cigar, index = ix, model = "within",
#'               effect = "twoways")
#' ab1 <- pangmm(lnc ~ lag(lnc, 1) + lnp + lnpn + lny | lag(lnc, 2:99),
#'               data = cigar, index = ix, transformation = "d",
#'               model = "onestep")
#' round(rbind(OLS = coef(ols)[-1], Within = coef(fe),
#'             "GMM one-step" = coef(ab1)[1:4]), 3)
#' ## published: 0.97 -0.090 0.024 -0.03 / 0.83 -0.299 0.034 0.10 /
#' ##            0.84 -0.377 -0.016 0.14
#'
#' ## long-run price elasticity, b2 / (1 - b1)
#' lr <- function(b) b[["lnp"]] / (1 - b[[1]])
#' c(within = lr(coef(fe)), gmm = lr(coef(ab1)))
#'
#' ## collapsed instruments: system GMM, two-step, Windmeijer errors
#' sys <- pangmm(lnc ~ lag(lnc, 1) + lnp + lnpn + lny | lag(lnc, 2:99),
#'               data = cigar, index = ix, transformation = "ld",
#'               model = "twostep", collapse = TRUE, maxlag = 6)
#' summary(sys)
#'
#' ## static demand: fixed vs random state effects (regression-based
#' ## Mundlak version of the Hausman test, robust to serial correlation)
#' sfe <- panglm(lnc ~ lnp + lnpn + lny, data = cigar, index = ix,
#'               model = "within")
#' panhausman(sfe, method = "mundlak", data = cigar)
#' @keywords datasets
"cigar"
