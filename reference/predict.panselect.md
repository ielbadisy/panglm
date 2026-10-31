# Predictions from a panel sample-selection fit

Predictions from a
[`panselect()`](https://elbadisyimad.com/panglm/reference/panselect.md)
fit on `newdata`, a long data frame with the panel index columns and the
variables of both equations (the outcome and the selection indicator
themselves are not needed).

## Usage

``` r
# S3 method for class 'panselect'
predict(
  object,
  newdata,
  type = c("unconditional", "conditional", "selection"),
  ...
)
```

## Arguments

- object:

  a `"panselect"` fit

- newdata:

  a data frame in long format (see Details); the fit does not store its
  data, so `newdata` is required (pass the fitting data for in-sample
  predictions)

- type:

  `"unconditional"`, `"conditional"` or `"selection"`

- ...:

  unused

## Value

a numeric vector with one prediction per row of `newdata`

## Details

- `type = "unconditional"`: the linear index of the outcome equation,
  `x'b`, the mean of the outcome in the whole population (selected or
  not) given the regressors. With `model = "pooling"` it includes the
  intercept and, with `mundlak = TRUE`, the Mundlak means; with
  `model = "within"` it excludes the individual effect (and the
  intercept, which the within transformation removes).

- `type = "conditional"`: the mean in the selected population,
  `x'b + sum_t g_t lambda_it`, adding the inverse Mills ratio terms with
  their estimated coefficients (one per period with `lambda = "period"`,
  a single one with `lambda = "common"`).

- `type = "selection"`: the selection probability `pnorm(z'd_t)` of the
  first-step probit of each row's period.

The inverse Mills ratios and probabilities use the first-step probit of
the row's period; the Mundlak means are the means of the time-varying
selection regressors over the rows of each individual in `newdata`, so
`newdata` should hold every period of an individual (as in the fit). A
period without a probit in the fit (every observation selected) has an
inverse Mills ratio of zero and a missing selection probability; a
period that is not in the fit gives a missing prediction for the
conditional and selection types.

## Validation

Numerical accuracy is checked in the package tests and in the validation
suite; see the validation article on the package website.

## See also

Other panselect methods:
[`summary.panselect()`](https://elbadisyimad.com/panglm/reference/summary.panselect.md),
[`vcov.panselect()`](https://elbadisyimad.com/panglm/reference/vcov.panselect.md)

## Examples

``` r
set.seed(1)
d <- data.frame(id = rep(1:150, each = 3), t = rep(1:3, 150))
d$x <- rnorm(450); d$z <- rnorm(450); v <- rnorm(450)
d$s <- as.numeric(0.5 + d$x + d$z + v > 0)
d$y <- ifelse(d$s == 1, 1 + d$x + 0.5 * v + rnorm(450), NA)
fit <- panselect(s ~ x + z, y ~ x, data = d, index = c("id", "t"),
                 model = "pooling", boot = 0)
head(predict(fit, newdata = d, type = "conditional"))
#> [1] 0.4389674 1.5181401 2.2310140 3.0856300 1.7378664 0.7881404
head(predict(fit, newdata = d, type = "selection"))
#> [1] 0.78550320 0.48020531 0.02973796 0.33258371 0.32511926 0.39075972
```
