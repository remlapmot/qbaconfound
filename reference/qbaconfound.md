# Monte Carlo quantitative bias analysis for unmeasured confounding

Conducts the flexible Monte Carlo quantitative bias analysis (QBA) for
unmeasured confounding of Hughes et al. Given a naive analysis model
(which omits one or more unmeasured confounders) and informative priors
for a small number of bias parameters, the function returns a
bias-adjusted estimate of the exposure effect together with an interval
that accounts for both the unmeasured confounding and sampling
variability.

## Usage

``` r
qbaconfound(
  formula,
  data,
  exposure = NULL,
  confounders,
  family = stats::gaussian(),
  reps = 1000L,
  sampling_error = TRUE,
  seed = NULL
)
```

## Arguments

- formula:

  The naive analysis model, e.g. `y ~ x + c1 + c2` for a GLM or
  `Surv(time, status) ~ x + c1` for a Cox model. The unmeasured
  confounders are *not* included in `formula`.

- data:

  A data frame containing the outcome, exposure, and measured
  confounders.

- exposure:

  Character vector naming the exposure term(s) in `formula` whose effect
  is of interest. Defaults to the first term on the right-hand side.

- confounders:

  A single
  [`u_continuous()`](https://remlapmot.github.io/qbaconfound/reference/u_continuous.md)/[`u_binary()`](https://remlapmot.github.io/qbaconfound/reference/u_binary.md)
  object, or a list of them, describing the unmeasured confounder(s).

- family:

  The outcome model family: a
  [`stats::glm()`](https://rdrr.io/r/stats/glm.html) family object, a
  family name, or `"cox"` for a Cox proportional hazards model. Ignored
  (and inferred as Cox) when the response is a
  [`survival::Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) call.

- reps:

  Number of Monte Carlo replications.

- sampling_error:

  Logical; if `TRUE` (the default) Monte Carlo sampling error is
  incorporated at step 4 above. Set to `FALSE` to obtain the
  distribution of bias-adjusted point estimates without sampling error.

- seed:

  Optional integer seed for reproducibility.

## Value

An object of class `qbaconfound`, a list with elements including
`estimates` (a data frame of naive and bias-adjusted estimates per
exposure term), `draws` (the matrix of bias-adjusted estimates across
replications), and `n_failed` (the number of replications whose model
fit failed).

## Details

The substantive analysis may be a generalised linear model (any
[`stats::glm()`](https://rdrr.io/r/stats/glm.html) `family`) or a Cox
proportional hazards model. Survival outcomes are detected automatically
when the left-hand side of `formula` is a
[`survival::Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) call,
or `family` can be set to `"cox"`.

Each Monte Carlo replication (see Hughes et al., section 2.4):

1.  draws a value for every bias parameter from its prior (the priors
    are specified through
    [`u_continuous()`](https://remlapmot.github.io/qbaconfound/reference/u_continuous.md)
    /
    [`u_binary()`](https://remlapmot.github.io/qbaconfound/reference/u_binary.md));

2.  simulates a proxy for each unmeasured confounder as a function of
    the exposure only (a continuous proxy from a normal model, a binary
    proxy from a Bernoulli model whose intercept reproduces the drawn
    prevalence);

3.  refits the outcome model including the simulated proxies, with their
    coefficients fixed to the drawn values (implemented as a model
    offset), and reads off the exposure coefficient and its standard
    error;

4.  adds Monte Carlo sampling error by drawing the bias-adjusted
    estimate from a normal distribution centred on that coefficient.

The point estimate is the median and the interval the 2.5th and 97.5th
percentiles of the resulting distribution of bias-adjusted estimates.

## References

Hughes RA, Kawabata E, Palmer TM, et al. A flexible Monte Carlo
quantitative bias analysis for unmeasured confounding. *Statistical
Methods in Medical Research* (under review).

## See also

[`u_continuous()`](https://remlapmot.github.io/qbaconfound/reference/u_continuous.md),
[`u_binary()`](https://remlapmot.github.io/qbaconfound/reference/u_binary.md),
[`sim_confounding()`](https://remlapmot.github.io/qbaconfound/reference/sim_confounding.md)

## Examples

``` r
df <- sim_confounding(n = 500, beta_x = 0, seed = 1)

# Naive model y ~ x + c1 omits the confounder u; adjust for one continuous U.
fit <- qbaconfound(
  y ~ x + c1, data = df, exposure = "x",
  confounders = u_continuous(coef_out = c(0.8, 0.1),
                             coef_exp = c(0.6, 0.05),
                             resid_sd = c(0.9, 1.1)),
  reps = 200, seed = 1
)
fit
#> Monte Carlo QBA for unmeasured confounding
#> Outcome family: gaussian | replications: 200 | observations: 500
#> Unmeasured confounder(s): 1 (continuous)
#> 
#>  term naive bias-adjusted          95% CI
#>     x 0.371        -0.099 (-0.305, 0.074)
```
