# Simulate data with an unmeasured confounder

A small helper that simulates a dataset with one measured confounder
`c1` and one unmeasured confounder `u` that jointly confound the
exposure-outcome relationship. The naive model `y ~ x + c1` (which omits
`u`) is therefore biased for the exposure effect, making the data useful
for examples and tests of
[`qbaconfound()`](https://remlapmot.github.io/qbaconfound/reference/qbaconfound.md).

## Usage

``` r
sim_confounding(
  n = 1000L,
  beta_x = 0,
  family = c("gaussian", "binomial"),
  seed = NULL
)
```

## Arguments

- n:

  Number of observations.

- beta_x:

  True exposure effect (on the linear-predictor scale).

- family:

  Outcome type: `"gaussian"` (continuous outcome) or `"binomial"`
  (binary outcome).

- seed:

  Optional integer seed for reproducibility.

## Value

A data frame with columns `y`, `x`, `c1`, and `u` (the unmeasured
confounder, included so it can be removed to mimic the unmeasured case).
The true exposure effect is stored in `attr(, "beta_x")`.

## Examples

``` r
df <- sim_confounding(n = 1000, beta_x = 0.5, seed = 42)
# Naive (biased) versus full (unbiased) model:
coef(lm(y ~ x + c1, df))["x"]
#>         x 
#> 0.8357274 
coef(lm(y ~ x + c1 + u, df))["x"]
#>        x 
#> 0.480904 
```
