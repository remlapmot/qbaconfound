# Summarise a Monte Carlo QBA

Summarise a Monte Carlo QBA

## Usage

``` r
# S3 method for class 'qbaconfound'
summary(object, ...)
```

## Arguments

- object:

  A `qbaconfound` object returned by
  [`qbaconfound()`](https://remlapmot.github.io/qbaconfound/reference/qbaconfound.md).

- ...:

  Unused.

## Value

A data frame of the naive and bias-adjusted estimates for each exposure
term, with columns `term`, `naive`, `naive_se`, `estimate`, `conf.low`,
and `conf.high`.
