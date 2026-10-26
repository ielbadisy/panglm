# Predictions from a latent class panel model

Predictions from a latent class panel model

## Usage

``` r
# S3 method for class 'panfmm'
predict(
  object,
  newdata = NULL,
  type = c("response", "class", "prior", "posterior", "components"),
  weights = c("prior", "posterior"),
  ...
)
```

## Arguments

- object:

  a `"panfmm"` fit

- newdata:

  optional data.frame; by default the estimation sample (in its original
  row order). Offsets are included in the class means: the
  estimation-sample offsets without `newdata`, the
  [`offset()`](https://rdrr.io/r/stats/offset.html) terms of the formula
  evaluated in `newdata` otherwise (an `offset` argument of
  [`panfmm()`](https://ielbadisy.github.io/panglm/reference/panfmm.md)
  cannot be evaluated in `newdata` and is not applied, with a warning)

- type:

  `"response"` for the marginal mean `sum_c w_ic mu_c(x_it)`; `"class"`
  for the modal class of each observation's individual (an integer, the
  class with the largest posterior probability in the estimation sample,
  or the largest prior probability, from the concomitant variables when
  there are some, with `newdata`); `"posterior"` / `"prior"` for the
  matrix of class probabilities of each observation's individual;
  `"components"` for the matrix of class-specific means `mu_c(x_it)`
  (one column per class). In panglm 1.5.0 and earlier, `type = "class"`
  returned the class-specific means now given by `"components"`; a
  message recalls this once per session.

- weights:

  for `type = "response"`: `"prior"` (the class probabilities, from the
  concomitant variables when there are some) or `"posterior"` (the
  in-sample posterior probabilities; only without `newdata`)

- ...:

  unused

## Value

a numeric vector (`"response"`), an integer vector of class labels
`1, ..., k` (`"class"`), or a matrix with one column per class
(`"posterior"`, `"prior"`, `"components"`)

## See also

Other panfmm methods:
[`summary.panfmm()`](https://ielbadisy.github.io/panglm/reference/summary.panfmm.md),
[`vcov.panfmm()`](https://ielbadisy.github.io/panglm/reference/vcov.panfmm.md)

## Examples

``` r
set.seed(3)
d <- data.frame(id = rep(1:150, each = 3), t = rep(1:3, 150), x = rnorm(450))
d$y <- rpois(450, exp(ifelse(rep(rbinom(150, 1, 0.4), each = 3) == 1, 1.5, 0) + 0.3 * d$x))
fit <- panfmm(y ~ x, data = d, index = c("id", "t"), nstart = 2, seed = 1)
head(predict(fit))
#> [1] 1.690137 2.067845 2.442469 1.596154 2.396391 2.279392
table(predict(fit, type = "class"))              # modal class per row
#> 
#>   1   2 
#> 282 168 
head(predict(fit, type = "posterior"))           # class probabilities
#>         class1     class2
#> [1,] 0.9112466 0.08875336
#> [2,] 0.9112466 0.08875336
#> [3,] 0.9112466 0.08875336
#> [4,] 0.9229608 0.07703924
#> [5,] 0.9229608 0.07703924
#> [6,] 0.9229608 0.07703924
predict(fit, newdata = data.frame(x = c(-1, 0, 1)), type = "components")
#>         class1   class2
#> [1,] 0.7443682 3.187864
#> [2,] 0.9513565 4.399313
#> [3,] 1.2159025 6.071134
```
