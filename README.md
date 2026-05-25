# qbaconfound

<!-- badges: start -->
<!-- badges: end -->

A flexible **Monte Carlo quantitative bias analysis (QBA)** for unmeasured
confounding in observational studies, as described in:

> Hughes RA, Kawabata E, Palmer TM, et al. *A flexible Monte Carlo quantitative
> bias analysis for unmeasured confounding.* Statistical Methods in Medical
> Research (forthcoming; manuscript SMM-25-0703).

When unmeasured confounding is suspected, this method produces a bias-adjusted
estimate of the exposure effect together with an interval that reflects both the
unmeasured confounding and sampling variability. It uses informative priors for
a **small number of bias parameters** that encode external information about the
unmeasured confounder(s) — and, unlike most existing software, applies to:

* a generalised linear model (any `glm` family) **or** a Cox proportional
  hazards model;
* a binary, continuous, or categorical exposure with measured confounders;
* one or more binary or continuous unmeasured confounders, which may be
  correlated with the measured confounders.

The number of bias parameters does **not** grow with the number of measured
confounders.

> **Status:** early development. The Monte Carlo QBA is implemented; the
> Bayesian QBA from the paper is not yet included. A companion Stata command of
> the same name is described in the paper.

## Installation

```r
# install.packages("remotes")
remotes::install_github("remlapmot/qbaconfound")
```

## How it works

Each Monte Carlo replication (paper, section 2.4):

1. draws a value for every bias parameter from its prior;
2. simulates a proxy for each unmeasured confounder *as a function of the
   exposure only*;
3. refits the outcome model including the simulated proxies, with their
   coefficients **fixed** to the drawn values (a model offset);
4. reads off the exposure coefficient and adds Monte Carlo sampling error.

The bias-adjusted point estimate is the median, and the interval the 2.5th and
97.5th percentiles, of the resulting distribution.

## Example

```r
library(qbaconfound)

# Simulate data where `u` confounds the x -> y relationship.
df <- sim_confounding(n = 1000, beta_x = 0.5, seed = 1)

# The naive model omits the unmeasured confounder u and is biased:
coef(lm(y ~ x + c1, df))["x"]        # ~0.88 (true value is 0.5)

# Adjust for a single continuous unmeasured confounder, supplying informative
# priors for its three bias parameters.
fit <- qbaconfound(
  y ~ x + c1, data = df, exposure = "x",
  confounders = u_continuous(
    coef_out = c(0.8, 0.1),   # N(mean, sd) prior for the U -> Y coefficient
    coef_exp = c(0.6, 0.05),  # N(mean, sd) prior for the X -> U coefficient
    resid_sd = c(0.9, 1.1)    # Uniform(min, max) prior for U's residual SD
  ),
  reps = 1000, seed = 1
)
fit
#> Monte Carlo QBA for unmeasured confounding
#> Outcome family: gaussian | replications: 1000 | observations: 1000
#> Unmeasured confounder(s): 1 (continuous)
#>
#>  term naive bias-adjusted       95% CI
#>     x 0.878         0.401 (0.22, 0.56)
```

The naive estimate (0.88) is badly biased away from the true value of 0.5; the
bias-adjusted estimate (0.40) is much closer and its interval covers the truth.

### Binary outcome, multiple unmeasured confounders

```r
qbaconfound(
  y ~ x + c1, data = df, exposure = "x", family = binomial(),
  confounders = list(
    u_continuous(coef_out = c(0.8, 0.1), coef_exp = c(0.6, 0.05),
                 resid_sd = c(0.9, 1.1)),
    u_binary(coef_out = c(0.4, 0.1), coef_exp = c(0.3, 0.05),
             prevalence = c(0.15, 0.25))
  ),
  reps = 1000, seed = 1
)
```

### Survival (Cox) outcome

```r
qbaconfound(
  survival::Surv(time, status) ~ x + c1, data = surv_df, exposure = "x",
  confounders = u_continuous(coef_out = c(0.5, 0.1), coef_exp = c(0.6, 0.05),
                             resid_sd = c(0.9, 1.1)),
  reps = 1000, seed = 1
)
```

## Specifying the bias parameters

| Constructor      | `coef_out` (β_U) | `coef_exp` (α_X) | third parameter            |
|------------------|------------------|------------------|----------------------------|
| `u_continuous()` | `c(mean, sd)`    | `c(mean, sd)`*   | `resid_sd` — `Uniform(min, max)` or `Gamma(shape, scale)` on the precision |
| `u_binary()`     | `c(mean, sd)`    | `c(mean, sd)`*   | `prevalence` — `Uniform(min, max)` or `Beta(a, b)` |

\* For a categorical exposure with several terms, pass `coef_exp` as a matrix
(one `c(mean, sd)` row per term) or a list of such vectors. Set a prior `sd` to
`0` for a fixed (point-mass) value.

## License

MIT
