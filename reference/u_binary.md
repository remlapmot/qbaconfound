# Specify a binary unmeasured confounder

Describes the prior distributions of the bias parameters for a single
binary unmeasured confounder, for use in
[`qbaconfound()`](https://remlapmot.github.io/qbaconfound/reference/qbaconfound.md).

## Usage

``` r
u_binary(coef_out, coef_exp, prevalence, prev_dist = c("uniform", "beta"))
```

## Arguments

- coef_out:

  Length-2 numeric `c(mean, sd)` giving the normal prior for the
  coefficient of the unmeasured confounder in the outcome model
  (\\\beta_U\\). Set `sd = 0` for a point-mass (fixed value) prior.

- coef_exp:

  Normal prior for the coefficient(s) of the exposure in the model for
  the unmeasured confounder (\\\alpha_X\\). For a single exposure term,
  a length-2 numeric `c(mean, sd)`. For a categorical exposure with
  several terms, a matrix with one row `c(mean, sd)` per exposure term,
  or a list of such length-2 vectors.

- prevalence:

  Length-2 numeric giving the prior for the marginal prevalence \\\pi\\
  of the binary unmeasured confounder. If `prev_dist = "uniform"` (the
  default) this is `c(min, max)` of a uniform prior; if
  `prev_dist = "beta"` this is `c(a, b)` of a beta prior.

- prev_dist:

  Prior family for the prevalence: `"uniform"` or `"beta"`.

## Value

An object of class `qba_u` describing a binary unmeasured confounder.

## Details

As for
[`u_continuous()`](https://remlapmot.github.io/qbaconfound/reference/u_continuous.md)
the coefficient of `U` in the outcome model (`coef_out`) and the
coefficient(s) of the exposure in the model for `U` (`coef_exp`) are
bias parameters. For a binary confounder the third bias parameter is the
marginal prevalence of `U` (`prevalence`, i.e. \\\pi\\) rather than a
residual standard deviation. A prevalence is usually easier to elicit
and more readily reported in the literature than the intercept of a
logistic model; the intercept needed to reproduce the drawn prevalence
is derived internally.

## See also

[`u_continuous()`](https://remlapmot.github.io/qbaconfound/reference/u_continuous.md),
[`qbaconfound()`](https://remlapmot.github.io/qbaconfound/reference/qbaconfound.md)

## Examples

``` r
u_binary(coef_out = c(0.7, 0.1), coef_exp = c(0.4, 0.05),
         prevalence = c(0.15, 0.25))
#> <binary unmeasured confounder>
#>   coef_out (beta_U): N(mean = 0.7, sd = 0.1)
#>   coef_exp (alpha_X)[1]: N(mean = 0.4, sd = 0.05)
#>   prevalence (pi): Uniform(0.15, 0.25)
```
