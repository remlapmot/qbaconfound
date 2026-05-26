# Specify a continuous unmeasured confounder

Describes the prior distributions of the bias parameters for a single
continuous unmeasured confounder, for use in
[`qbaconfound()`](https://remlapmot.github.io/qbaconfound/reference/qbaconfound.md).

## Usage

``` r
u_continuous(coef_out, coef_exp, resid_sd, resid_dist = c("uniform", "gamma"))
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

- resid_sd:

  Length-2 numeric giving the prior for the residual standard deviation
  \\\eta\\. If `resid_dist = "uniform"` (the default) this is
  `c(min, max)` of a uniform prior on \\\eta\\. If
  `resid_dist = "gamma"` this is `c(shape, scale)` of a gamma prior on
  the residual *precision* \\1/\eta^2\\.

- resid_dist:

  Prior family for the residual standard deviation: either `"uniform"`
  (on \\\eta\\) or `"gamma"` (on the precision \\1/\eta^2\\).

## Value

An object of class `qba_u` describing a continuous unmeasured
confounder.

## Details

The bias model relates the unmeasured confounder `U` to the study data
through three bias parameters: the coefficient of `U` in the outcome
model (`coef_out`, i.e. \\\beta_U\\), the coefficient(s) of the exposure
in the model for `U` (`coef_exp`, i.e. \\\alpha_X\\), and the residual
standard deviation of `U` given the exposure (`resid_sd`, i.e.
\\\eta\\). Values for these parameters cannot be estimated from the data
and so are drawn from the prior distributions specified here.

## See also

[`u_binary()`](https://remlapmot.github.io/qbaconfound/reference/u_binary.md),
[`qbaconfound()`](https://remlapmot.github.io/qbaconfound/reference/qbaconfound.md)

## Examples

``` r
u_continuous(coef_out = c(0.8, 0.1), coef_exp = c(0.3, 0.05),
             resid_sd = c(0.9, 1.1))
#> <continuous unmeasured confounder>
#>   coef_out (beta_U): N(mean = 0.8, sd = 0.1)
#>   coef_exp (alpha_X)[1]: N(mean = 0.3, sd = 0.05)
#>   resid_sd (eta): Uniform(0.9, 1.1)
```
