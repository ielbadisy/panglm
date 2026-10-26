# Family specifications for panglm

Lightweight family objects used by
[`panglm()`](https://elbadisyimad.com/panglm/reference/panglm.md),
analogous in spirit to
[stats::family](https://rdrr.io/r/stats/family.html) objects but keyed
to the integer codes the C++ backend expects (see `src/families.h`).

## Usage

``` r
gaussianfamily(link = "identity")

poissonfamily(link = "log")

binomialfamily(link = c("logit", "probit"))

negbinfamily(link = "log")

gammafamily(link = "log")

orderedfamily(link = c("logit", "probit"))

betafamily(link = c("logit", "probit"))

tobitfamily(left = 0, right = Inf)
```

## Arguments

- link:

  canonical link name

- left, right:

  censoring limits of `tobitfamily()`

## Value

a list with `family`, `link`, `family_id`, `link_id`

## Details

`gammafamily()` is the gamma GLM with a log link, the usual choice for
skewed positive outcomes such as health care costs: the variance is
proportional to the squared mean and the coefficients act
multiplicatively on the mean, with no retransformation needed.

`orderedfamily()` is the ordered logit or ordered probit model for
ordinal outcomes such as self-assessed health:
`P(y <= j) = F(c_j - x'beta - u_i)`. The response may be a factor, an
ordered factor or integer-valued; categories follow the factor levels or
the sorted values. The intercept is absorbed in the cutpoints `c_j`,
which are reported after the slopes. `model = "pooling"` and
`model = "random"` (Gaussian random intercept by adaptive Gauss-Hermite
quadrature) are available for both links, `model = "within"` for the
logit link only (the "blow-up and cluster" fixed-effects ordered logit
of Baetschmann, Staub and Winkelmann 2015). The fit has class
`c("panordered", "panglm")`; see
[`predict.panordered()`](https://elbadisyimad.com/panglm/reference/predict.panordered.md).

`betafamily()` is beta regression for continuous outcomes in the open
interval (0, 1), such as utility index scores rescaled to (0, 1),
proportions of days in hospital or adherence fractions:
`y ~ Beta(mu * phi, (1 - mu) * phi)` with `g(mu) = x'beta`, so
`E(y) = mu` and `Var(y) = mu (1 - mu) / (1 + phi)`. The precision `phi`
is reported as `fit$precision`. `model = "pooling"` and
`model = "random"` (Gaussian random intercept on the link scale,
adaptive Gauss-Hermite quadrature) are available; there is no
fixed-effects estimator. Outcomes equal to 0 or 1 must be transformed
first, e.g. with `(y * (n - 1) + 0.5) / n` (Smithson and Verkuilen
2006).

`tobitfamily()` is the censored normal (Tobit) model for outcomes piled
up at a floor or a ceiling, such as health scores at the ceiling of a
scale (motor function scores, EQ-5D utility at 1) or spending at zero:
the latent outcome is `y* = x'beta + u_i + e_it` with
`e_it ~ N(0, sigma_e^2)`, and `y = left` when `y* <= left`, `y = right`
when `y* >= right`, `y = y*` otherwise. Use `left = -Inf` or
`right = Inf` for one-sided censoring. `model = "pooling"` and
`model = "random"` (Gaussian random intercept `u_i ~ N(0, sigma_u^2)`,
adaptive Gauss-Hermite quadrature) are available. There is no
fixed-effects Tobit (incidental parameters); use `model = "random"` with
`cre = TRUE`. `predict(type = "link")` returns the latent mean `x'beta`
and `predict(type = "response")` the mean of the observed, censored
outcome (averaged over the random intercept), which is the scale of
[`partialeffects()`](https://elbadisyimad.com/panglm/reference/partialeffects.md).
`fit$sigma_e` and `fit$sigma` (random intercept SD) are stored on the
fit.

## References

Baetschmann, G., Staub, K. E. and Winkelmann, R. (2015). Consistent
estimation of the fixed effects ordered logit model. *Journal of the
Royal Statistical Society A*, 178, 685-703.
[doi:10.1111/rssa.12090](https://doi.org/10.1111/rssa.12090)

Ferrari, S. and Cribari-Neto, F. (2004). Beta regression for modelling
rates and proportions. *Journal of Applied Statistics*, 31, 799-815.
[doi:10.1080/0266476042000214501](https://doi.org/10.1080/0266476042000214501)

Smithson, M. and Verkuilen, J. (2006). A better lemon squeezer?
Maximum-likelihood regression with beta-distributed dependent variables.
*Psychological Methods*, 11, 54-71.
[doi:10.1037/1082-989X.11.1.54](https://doi.org/10.1037/1082-989X.11.1.54)

Tobin, J. (1958). Estimation of relationships for limited dependent
variables. *Econometrica*, 26, 24-36.
[doi:10.2307/1907382](https://doi.org/10.2307/1907382)

## See also

[`panglm()`](https://elbadisyimad.com/panglm/reference/panglm.md), whose
`family` argument takes these objects, and
[panglm-deprecated](https://elbadisyimad.com/panglm/reference/panglm-deprecated.md)
for the former names.

## Examples

``` r
gaussianfamily()
#> $family
#> [1] "gaussian"
#> 
#> $link
#> [1] "identity"
#> 
#> $family_id
#> [1] 0
#> 
#> $link_id
#> [1] 0
#> 
poissonfamily()
#> $family
#> [1] "poisson"
#> 
#> $link
#> [1] "log"
#> 
#> $family_id
#> [1] 1
#> 
#> $link_id
#> [1] 1
#> 
binomialfamily("probit")
#> $family
#> [1] "binomial"
#> 
#> $link
#> [1] "probit"
#> 
#> $family_id
#> [1] 2
#> 
#> $link_id
#> [1] 3
#> 
negbinfamily()
#> $family
#> [1] "negbin"
#> 
#> $link
#> [1] "log"
#> 
#> $family_id
#> [1] 1
#> 
#> $link_id
#> [1] 1
#> 
gammafamily()
#> $family
#> [1] "gamma"
#> 
#> $link
#> [1] "log"
#> 
#> $family_id
#> [1] 3
#> 
#> $link_id
#> [1] 1
#> 
orderedfamily("probit")
#> $family
#> [1] "ordered"
#> 
#> $link
#> [1] "probit"
#> 
#> $family_id
#> [1] NA
#> 
#> $link_id
#> [1] 3
#> 
betafamily()
#> $family
#> [1] "beta"
#> 
#> $link
#> [1] "logit"
#> 
#> $family_id
#> [1] NA
#> 
#> $link_id
#> [1] 2
#> 
tobitfamily(left = 0)
#> $family
#> [1] "tobit"
#> 
#> $link
#> [1] "identity"
#> 
#> $family_id
#> [1] NA
#> 
#> $link_id
#> [1] 0
#> 
#> $left
#> [1] 0
#> 
#> $right
#> [1] Inf
#> 

data(copd)
fit <- panglm(hospitalized ~ crp, data = copd, index = c("id", "visit"),
              model = "within", family = binomialfamily("logit"))
#> 17 group(s) with no within-group outcome variation dropped (all-0 or all-1) -- they carry no information for the conditional likelihood.
coef(fit)
#>        crp 
#> 0.09667869 
```
