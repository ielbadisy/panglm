# Counts, costs, bounded and censored outcomes

Health outcomes are rarely Gaussian. This vignette shows the `panglm`
models for the common cases:

| Outcome | Example | Model |
|----|----|----|
| Count | GP visits, exacerbations | Poisson or negative binomial (NB2), `family = "poisson"` or `"negbin"` |
| Count with many zeros | exacerbations | hurdle, [`panhurdle()`](https://elbadisyimad.com/panglm/reference/panhurdle.md) |
| Cost with zeros, skewed when positive | annual expenditure | two-part model, [`twopart()`](https://elbadisyimad.com/panglm/reference/twopart.md) |
| Proportion in (0, 1) | proportion of days covered | beta, `family = "beta"` |
| Ordered category | severity grade | ordered logit or probit, [`orderedfamily()`](https://elbadisyimad.com/panglm/reference/panglm-families.md) |
| Censored at a bound | EQ-5D index at 1 | tobit, [`tobitfamily()`](https://elbadisyimad.com/panglm/reference/panglm-families.md) |

All of them use the same `index`, `model`, `vcov`,
[`partialeffects()`](https://elbadisyimad.com/panglm/reference/partialeffects.md)
and [`predict()`](https://rdrr.io/r/stats/predict.html) interface as the
Gaussian model of
[`vignette("panglm-intro")`](https://elbadisyimad.com/panglm/articles/panglm-intro.md).
The examples use small samples so that the vignette builds in seconds;
the companion book uses the full data and more quadrature nodes (`R`).

``` r

library(panglm)
```

## Counts: Poisson and negative binomial

We simulate a panel of 200 patients over 5 years. The number of GP
visits depends on a chronic condition, age and sex, plus a patient
effect `a` that also raises the chance of having the condition, so
pooled models are biased. The counts are overdispersed (NB2 with size
2).

``` r

set.seed(42)
N <- 200; TT <- 5
d <- data.frame(id = rep(1:N, each = TT), year = rep(1:TT, N))
a <- rnorm(N, sd = 0.6)
d$age <- rep(round(runif(N, 50, 80)), each = TT) + d$year - 1
d$female <- rep(rbinom(N, 1, 0.5), each = TT)
d$chronic <- rbinom(N * TT, 1, plogis(-0.5 + a[d$id]))
mu <- exp(0.5 * d$chronic + 0.03 * (d$age - 65) + 0.2 * d$female + a[d$id])
d$visits <- rnbinom(N * TT, mu = mu, size = 2)
```

**Fixed-effects Poisson** is the conditional Poisson estimator of
Hausman, Hall and Griliches (1984). It is consistent for the slopes
whenever the conditional mean is correct, whatever the distribution of
the counts, so with a cluster covariance it is robust to overdispersion
and serial correlation (Wooldridge 1999). This makes it the default
choice for panel counts.

``` r

fe_pois <- panglm(visits ~ chronic + age, data = d, index = c("id", "year"),
                  model = "within", family = "poisson", vcov = "cluster")
summary(fe_pois)
#> 
#> Call:
#> panglm(formula = visits ~ chronic + age, data = d, index = c("id", 
#>     "year"), model = "within", family = "poisson", vcov = "cluster")
#> 
#> Family: poisson  Link: log  Model: within  Effect: individual 
#> N = 1000  groups = 200  vcov: cluster 
#> 
#> Coefficients:
#>         Estimate Std. Error z value Pr(>|z|)    
#> chronic  0.50872    0.08120   6.265 3.72e-10 ***
#> age      0.01126    0.02643   0.426     0.67    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Log-likelihood: -1224 
#> Dispersion parameter: 1 
#> Residual df: 798  Iterations: 3
```

The true effect is 0.5: a chronic condition raises the visit rate by a
factor of about $`\exp(0.51) = 1.66`$. Sex does not vary within a
patient, so it is not in the model.

The dispersion test of Cameron and Trivedi (1990) checks whether the
variance exceeds the Poisson mean. A ratio well above 1 calls for
cluster standard errors, or for a negative binomial model when
efficiency or predicted probabilities matter.

``` r

pandispersiontest(fe_pois)
#> 
#>  Pearson chi-squared/df overdispersion test
#> 
#> chisq = 1367 , df = 773 , ratio = 1.768 , p-value = < 2.2e-16 
#> alternative hypothesis: overdispersion (ratio > 1)
#> 25 observation(s) with a zero or unavailable fitted value excluded
```

**Negative binomial** fits add the NB2 variance
$`\mu + \mu^2 / \theta`$. The fixed-effects version uses one intercept
per patient (the incidental parameter bias is small for count models);
the random-effects version is the beta-negative-binomial model of
Hausman, Hall and Griliches (1984). Random effects keep time-invariant
regressors such as sex.

``` r

fe_nb <- panglm(visits ~ chronic + age, data = d, index = c("id", "year"),
                model = "within", family = "negbin")
#> 5 group(s) with all-zero outcomes dropped -- their fixed effect is unbounded under the dummy-variable NB2 MLE and carries no information about the shared slope.
coef(fe_nb)
#>    chronic        age 
#> 0.53822486 0.02721607
fe_nb$theta
#> [1] 3.810086
re_nb <- panglm(visits ~ chronic + age + female, data = d,
                index = c("id", "year"), model = "random", family = "negbin")
coef(re_nb)
#> (Intercept)     chronic         age      female 
#> -1.16617986  0.62950329  0.02170266  0.16128873
```

The random-effects Poisson (gamma mixing) and negative binomial models
assume that the patient effect is independent of the regressors. Here it
is not: patients with a high `a` are more often chronic, and the
random-effects NB slope of `chronic` (0.63) is above the fixed-effects
one (0.54) and the true 0.5. `cre = TRUE` removes the bias from the
time-varying regressors, as in the Gaussian case.

``` r

re_pois <- panglm(visits ~ chronic + age + female, data = d,
                  index = c("id", "year"), model = "random",
                  family = "poisson", cre = TRUE, vcov = "cluster")
partialeffects(re_pois)
#> 
#> Average partial effects
#> Delta-method standard errors, vcov: cluster 
#> 
#>     term contrast estimate std.error statistic   p.value conf.low conf.high
#>  chronic    1 - 0  0.95880   0.15190    6.3100 2.786e-10  0.66100    1.2570
#>      age    dY/dX  0.02452   0.04927    0.4977 6.187e-01 -0.07204    0.1211
#>   female    1 - 0  0.29200   0.18370    1.5900 1.119e-01 -0.06803    0.6520
```

The average partial effect is on the count scale: a chronic condition
adds about one visit a year on average.

## Excess zeros: the hurdle model

A hurdle model separates whether any event happens from how many events
happen once the hurdle is crossed.
[`panhurdle()`](https://elbadisyimad.com/panglm/reference/panhurdle.md)
fits a conditional (fixed-effects) logit for $`P(y > 0)`$ and a
zero-truncated Poisson or negative binomial fixed-effects model for the
positive counts. Patients whose outcome never changes between zero and
positive carry no information for the binary part and are dropped from
it.

``` r

data(copd)
hur <- panhurdle(exacerbations ~ crp, data = copd, index = c("id", "visit"))
#> 26 group(s) with no within-group outcome variation dropped (all-0 or all-1) -- they carry no information for the conditional likelihood.
summary(hur)
#> 
#> Formula: exacerbations ~ crp
#> 
#> Estimator: fixed-effects hurdle model (conditional logit and zero-truncated Poisson)
#> N = 300  positive = 238  vcov: model-based (inverse Hessian)
#> 
#> Zero part, 1(y > 0), conditional logit (groups used = 34):
#>     Estimate Std. Error z value Pr(>|z|)
#> crp  0.10889    0.08173   1.332    0.183
#> 
#> Count part, y | y > 0, zero-truncated Poisson (groups used = 59):
#>     Estimate Std. Error z value Pr(>|z|)  
#> crp  0.04246    0.02485   1.709   0.0875 .
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

Within a patient, higher CRP goes with higher odds of any exacerbation
and with more exacerbations among visits with at least one, but with 60
patients neither effect is precisely estimated. Only 34 patients switch
between zero and positive counts, which is typical of short panels.

## Costs: the two-part model

Health care costs have a mass at zero and a long right tail. The
two-part model of Duan et al. (1983) fits a probit or logit for
$`P(y > 0)`$ on all observations and a gamma GLM with a log link for
$`E(y \mid y > 0)`$ on the positive costs, so that
$`E(y) = P(y > 0)\, E(y \mid y > 0)`$ needs no retransformation (Manning
and Mullahy 2001). `positive = "lognormal"` fits a linear model for
$`\log y`$ with smearing instead.

``` r

set.seed(11)
N <- 150; TT <- 4
cst <- data.frame(id = rep(1:N, each = TT), year = rep(1:TT, N))
a <- rnorm(N, sd = 0.5)
cst$age <- rep(round(runif(N, 40, 85)), each = TT)
cst$chronic <- rbinom(N * TT, 1, plogis(-0.5 + a[cst$id]))
p_any <- pnorm(-0.3 + 0.6 * cst$chronic + 0.01 * (cst$age - 60) + a[cst$id])
mu <- exp(7 + 0.5 * cst$chronic + 0.01 * (cst$age - 60) + 0.5 * a[cst$id])
cst$cost <- ifelse(runif(N * TT) < p_any,
                   rgamma(N * TT, shape = 2, rate = 2 / mu), 0)
mean(cst$cost == 0)
#> [1] 0.5433333
```

``` r

tp <- twopart(cost ~ chronic + age, data = cst, index = c("id", "year"),
              model = "pooling", vcov = "cluster")
summary(tp)
#> 
#> Call:
#> twopart(formula = cost ~ chronic + age, data = cst, index = c("id", 
#>     "year"), model = "pooling", vcov = "cluster")
#> 
#> Estimator: two-part model (pooling)
#> N = 600  groups = 150  positive = 274  vcov: cluster 
#> 
#> Part 1: P(y > 0), binomial probit 
#>              Estimate Std. Error z value Pr(>|z|)    
#> (Intercept) -1.256927   0.341659  -3.679 0.000234 ***
#> chronic      0.751339   0.107478   6.991 2.74e-12 ***
#> age          0.013565   0.005223   2.597 0.009393 ** 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Part 2: E[y | y > 0], gamma GLM with log link 
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  6.30418    0.28162  22.385  < 2e-16 ***
#> chronic      0.52261    0.09123   5.729 1.01e-08 ***
#> age          0.01292    0.00426   3.032  0.00243 ** 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> Log-likelihood: -2664
partialeffects(tp)
#> 
#> Average partial effects (total mean)
#> Delta-method standard errors, vcov: cluster 
#> 
#>     term contrast estimate std.error statistic   p.value conf.low conf.high
#>  chronic    1 - 0   903.00   124.000     7.281 3.324e-13   659.90   1146.00
#>      age    dY/dX    18.11     4.666     3.880 1.044e-04     8.96     27.25
```

The two parts are on different scales;
[`partialeffects()`](https://elbadisyimad.com/panglm/reference/partialeffects.md)
combines them into effects on expected cost. A chronic condition adds
about 900 currency units a year: part of it because more patients have
any cost, part because their costs are higher. `model = "random"` adds a
random intercept to each part (and `correlated = TRUE` lets the two
intercepts correlate, for the lognormal part).

## Proportions: the beta model

Outcomes strictly between 0 and 1, such as the proportion of days
covered by a prescription, follow a beta regression with a logit link
and precision $`\phi`$ (Ferrari and Cribari-Neto 2004).
`model = "random"` adds a normal random intercept on the link scale.

``` r

m <- plogis(1 - 0.4 * cst$chronic - 0.01 * (cst$age - 60) + a[cst$id])
cst$pdc <- rbeta(N * TT, m * 20, (1 - m) * 20)
cst$age10 <- (cst$age - 60) / 10
fit_beta <- panglm(pdc ~ chronic + age10, data = cst, model = "random",
                   index = c("id", "year"), family = "beta", R = 5)
summary(fit_beta)
#> 
#> Call:
#> panglm(formula = pdc ~ chronic + age10, data = cst, index = c("id", 
#>     "year"), model = "random", family = "beta", R = 5)
#> 
#> Family: beta  Link: logit  Model: random  Effect: individual 
#> N = 600  groups = 150  vcov: classical 
#> 
#> Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  0.91896    0.04849  18.953  < 2e-16 ***
#> chronic     -0.31919    0.04581  -6.968 3.22e-12 ***
#> age10       -0.08485    0.03516  -2.413   0.0158 *  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Log-likelihood: 428.7 
#> Dispersion parameter: 0.04733 
#> Random-intercept SD: 0.4776 
#> Residual df: 595  Iterations: 39
partialeffects(fit_beta)
#> 
#> Average partial effects
#> Delta-method standard errors, vcov: classical 
#> 
#>     term contrast estimate std.error statistic   p.value conf.low conf.high
#>  chronic    1 - 0 -0.06721  0.009739    -6.902 5.134e-12 -0.08630 -0.048130
#>    age10    dY/dX -0.01764  0.007290    -2.420 1.551e-02 -0.03193 -0.003355
```

Patients with a chronic condition cover about 7 percentage points fewer
days. Age is rescaled to decades from 60: the beta fit is sensitive to
badly scaled regressors, so centring and scaling them is good practice.
Values of exactly 0 or 1 must be moved inside the interval first, or
modelled with a two-part model.

## Ordered outcomes

Severity grades and satisfaction scales are ordinal.
`family = orderedfamily("logit")` (or `"probit"`) fits the cumulative
link model with estimated cut points. The random-effects version adds a
normal random intercept; the fixed-effects ordered logit uses the
“blow-up and cluster” estimator of Baetschmann, Staub and Winkelmann
(2015), which is consistent for the slopes.

``` r

copd$severity <- cut(copd$exacerbations, c(-Inf, 0, 2, Inf),
                     labels = c("none", "mild", "severe"),
                     ordered_result = TRUE)
fit_ord <- panglm(severity ~ treatment + age + crp, data = copd,
                  index = c("id", "visit"), model = "random",
                  family = orderedfamily("logit"), R = 8)
summary(fit_ord)
#> 
#> Call:
#> panglm(formula = severity ~ treatment + age + crp, data = copd, 
#>     index = c("id", "visit"), model = "random", family = orderedfamily("logit"), 
#>     R = 8)
#> 
#> Family: ordered  Link: logit  Model: random  Effect: individual 
#> N = 300  groups = 60  vcov: classical 
#> 
#> Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)    
#> treatment   -0.70419    0.42706  -1.649 0.099165 .  
#> age          0.09392    0.02538   3.700 0.000216 ***
#> crp          0.07963    0.06376   1.249 0.211707    
#> none|mild    4.27286    1.70980   2.499 0.012453 *  
#> mild|severe  7.13186    1.75713   4.059 4.93e-05 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Log-likelihood: -279.5 
#> Random-intercept SD: 1.292 
#> Residual df: 294  Iterations: 42
fe_ord <- panglm(severity ~ crp, data = copd, index = c("id", "visit"),
                 model = "within", family = orderedfamily("logit"))
coef(fe_ord)
#>       crp 
#> 0.1083205
```

[`partialeffects()`](https://elbadisyimad.com/panglm/reference/partialeffects.md)
reports effects on the probability of each category, and
`predict(type = "prob")` gives the predicted probabilities.

``` r

partialeffects(fit_ord)
#> 
#> Average partial effects (category probability mean)
#> Delta-method standard errors, vcov: classical 
#> 
#>  category      term contrast  estimate std.error statistic   p.value  conf.low
#>      none treatment    1 - 0  0.086690  0.053490     1.621 1.051e-01 -0.018140
#>      mild treatment    1 - 0  0.024270  0.017670     1.374 1.695e-01 -0.010360
#>    severe treatment    1 - 0 -0.111000  0.066380    -1.671 9.463e-02 -0.241100
#>      none       age    dY/dX -0.011320  0.002998    -3.775 1.601e-04 -0.017190
#>      mild       age    dY/dX -0.003486  0.001651    -2.112 3.469e-02 -0.006722
#>    severe       age    dY/dX  0.014800  0.003596     4.117 3.837e-05  0.007757
#>      none       crp    dY/dX -0.009596  0.007671    -1.251 2.109e-01 -0.024630
#>      mild       crp    dY/dX -0.002956  0.002674    -1.106 2.689e-01 -0.008196
#>    severe       crp    dY/dX  0.012550  0.009992     1.256 2.090e-01 -0.007032
#>   conf.high
#>   0.1915000
#>   0.0588900
#>   0.0191500
#>  -0.0054420
#>  -0.0002508
#>   0.0218500
#>   0.0054380
#>   0.0022840
#>   0.0321400
head(predict(fit_ord, type = "prob"), 3)
#>           none      mild    severe
#> [1,] 0.1674925 0.4782422 0.3542652
#> [2,] 0.2257297 0.4966803 0.2775900
#> [3,] 0.1843286 0.4861148 0.3295566
```

Active treatment lowers the probability of the “severe” category by
about 11 percentage points and raises that of “none” by about 9
(imprecisely estimated); the effects across categories sum to zero.

## Censored outcomes: the tobit model

The EQ-5D index is capped at 1 (full health), and many patients sit at
the cap. A linear model then underestimates the effects, because the
outcome cannot rise above 1 whatever the latent health. The tobit model
(Tobin 1958) treats the observed value as a censored latent variable;
`tobitfamily(left, right)` sets the bounds. Pooled and random-effects
tobit are available; the fixed-effects tobit is inconsistent with short
panels, so use `cre = TRUE` instead.

``` r

set.seed(7)
N <- 80; TT <- 4
eq <- data.frame(id = rep(1:N, each = TT), wave = rep(1:TT, N))
u <- rnorm(N, sd = 0.15)
eq$age <- rep(round(runif(N, 40, 85)), each = TT)
eq$chronic <- rpois(N * TT, 0.5 + 0.02 * (eq$age - 40))
ystar <- 1.15 - 0.006 * (eq$age - 40) - 0.08 * eq$chronic + u[eq$id] +
  rnorm(N * TT, sd = 0.12)
eq$eq5d <- pmin(ystar, 1)
mean(eq$eq5d == 1)
#> [1] 0.428125
```

``` r

fit_tob <- panglm(eq5d ~ age + chronic, data = eq, index = c("id", "wave"),
                  model = "random", family = tobitfamily(left = -Inf, right = 1),
                  R = 5)
coef(fit_tob)
#>  (Intercept)          age      chronic 
#>  1.567750520 -0.008401564 -0.082276325
c(sigma_e = fit_tob$sigma_e, sigma_u = fit_tob$sigma)
#>   sigma_e   sigma_u 
#> 0.1177384 0.1346040
fit_lin <- panglm(eq5d ~ age + chronic, data = eq, index = c("id", "wave"),
                  model = "random")
coef(fit_lin)
#>  (Intercept)          age      chronic 
#>  1.254339319 -0.005045475 -0.053889729
```

The tobit slopes refer to latent health. The linear random-effects
slopes are attenuated by the ceiling, by about 40% here, since 43% of
the observations sit at 1. Use
[`partialeffects()`](https://elbadisyimad.com/panglm/reference/partialeffects.md)
for effects on the observed, censored outcome.

## Validated against

Fixed-effects Poisson reproduces `pglm` and
[`fixest::fepois()`](https://lrberge.github.io/fixest/reference/feglm.html);
negative binomial fits reproduce
[`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html) and `pglm`;
random-effects Poisson and negative binomial reproduce `pglm`; the
hurdle count part reproduces
[`pscl::hurdle()`](https://rdrr.io/pkg/pscl/man/hurdle.html) on the same
data; the two-part gamma part reproduces
[`stats::glm()`](https://rdrr.io/r/stats/glm.html) and `MASS`; beta
models reproduce `betareg` and `glmmTMB`; ordered models reproduce
[`MASS::polr()`](https://rdrr.io/pkg/MASS/man/polr.html),
[`ordinal::clmm()`](https://rdrr.io/pkg/ordinal/man/clmm.html) and, for
the fixed-effects logit,
[`survival::clogit()`](https://rdrr.io/pkg/survival/man/clogit.html);
tobit models reproduce
[`AER::tobit()`](https://rdrr.io/pkg/AER/man/tobit.html) and `censReg`.
The scripts are in `tests/testthat` and `validation/`.

## References

Baetschmann, G., Staub, K. E., and Winkelmann, R. (2015). Consistent
estimation of the fixed effects ordered logit model. *Journal of the
Royal Statistical Society A*, 178, 685-703.

Cameron, A. C., and Trivedi, P. K. (1990). Regression-based tests for
overdispersion in the Poisson model. *Journal of Econometrics*, 46,
347-364.

Duan, N., Manning, W. G., Morris, C. N., and Newhouse, J. P. (1983). A
comparison of alternative models for the demand for medical care.
*Journal of Business and Economic Statistics*, 1, 115-126.

Ferrari, S., and Cribari-Neto, F. (2004). Beta regression for modelling
rates and proportions. *Journal of Applied Statistics*, 31, 799-815.

Hausman, J., Hall, B. H., and Griliches, Z. (1984). Econometric models
for count data with an application to the patents-R&D relationship.
*Econometrica*, 52, 909-938.

Manning, W. G., and Mullahy, J. (2001). Estimating log models: to
transform or not to transform? *Journal of Health Economics*, 20,
461-494.

Tobin, J. (1958). Estimation of relationships for limited dependent
variables. *Econometrica*, 26, 24-36.

Wooldridge, J. M. (1999). Distribution-free estimation of some nonlinear
panel data models. *Journal of Econometrics*, 90, 77-97.
