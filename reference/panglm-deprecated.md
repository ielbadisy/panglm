# Deprecated functions in panglm

These functions still work but will be removed in a future release. Each
one emits a deprecation warning once per session and calls its
replacement, so it returns exactly the same result.

## Usage

``` r
gaussian_family(link = "identity")

binomial_family(link = c("logit", "probit"))

poisson_family(link = "log")

negbin_family(link = "log")

gamma_family(link = "log")

beta_family(link = c("logit", "probit"))

ordered_family(link = c("logit", "probit"))

tobit_family(left = 0, right = Inf)

panglm_dispersiontest(object)

panglm_hausman(
  fe,
  re = NULL,
  method = c("classical", "mundlak"),
  data = NULL,
  vcov = c("cluster", "HC1", "classical")
)

panglm_hurdle(
  formula,
  data,
  index,
  count_family = c("poisson", "negbin"),
  maxit = 100,
  tol = 1e-10
)
```

## Arguments

- link, left, right:

  passed to the replacement

- object, fe, re, method, data, vcov, formula, index, maxit, tol:

  passed to the replacement

- count_family:

  passed as `family` to
  [`panhurdle()`](https://elbadisyimad.com/panglm/reference/panhurdle.md)

## Value

the value of the replacement function

## Details

Exported function names in panglm are glued lowercase, without
underscores:

|  |  |
|----|----|
| **Deprecated** | **Replacement** |
| `gaussian_family()` | [`gaussianfamily()`](https://elbadisyimad.com/panglm/reference/panglm-families.md) |
| `binomial_family()` | [`binomialfamily()`](https://elbadisyimad.com/panglm/reference/panglm-families.md) |
| `poisson_family()` | [`poissonfamily()`](https://elbadisyimad.com/panglm/reference/panglm-families.md) |
| `negbin_family()` | [`negbinfamily()`](https://elbadisyimad.com/panglm/reference/panglm-families.md) |
| `gamma_family()` | [`gammafamily()`](https://elbadisyimad.com/panglm/reference/panglm-families.md) |
| `beta_family()` | [`betafamily()`](https://elbadisyimad.com/panglm/reference/panglm-families.md) |
| `ordered_family()` | [`orderedfamily()`](https://elbadisyimad.com/panglm/reference/panglm-families.md) |
| `tobit_family()` | [`tobitfamily()`](https://elbadisyimad.com/panglm/reference/panglm-families.md) |
| `panglm_dispersiontest()` | [`pandispersiontest()`](https://elbadisyimad.com/panglm/reference/pandispersiontest.md) |
| `panglm_hausman()` | [`panhausman()`](https://elbadisyimad.com/panglm/reference/panhausman.md) |
| `panglm_hurdle()` | [`panhurdle()`](https://elbadisyimad.com/panglm/reference/panhurdle.md) |

`panglm_hurdle(count_family = )` is `panhurdle(family = )`. The objects
returned keep their classes (`"panglm_dispersiontest"`,
`"panglm_hausman"`, `"panglm_hurdle"`), so code that tests or dispatches
on them is unaffected.
