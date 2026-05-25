#' Simulate data with an unmeasured confounder
#'
#' A small helper that simulates a dataset with one measured confounder `c1`
#' and one unmeasured confounder `u` that jointly confound the exposure-outcome
#' relationship. The naive model `y ~ x + c1` (which omits `u`) is therefore
#' biased for the exposure effect, making the data useful for examples and
#' tests of [qbaconfound()].
#'
#' @param n Number of observations.
#' @param beta_x True exposure effect (on the linear-predictor scale).
#' @param family Outcome type: `"gaussian"` (continuous outcome) or
#'   `"binomial"` (binary outcome).
#' @param seed Optional integer seed for reproducibility.
#'
#' @return A data frame with columns `y`, `x`, `c1`, and `u` (the unmeasured
#'   confounder, included so it can be removed to mimic the unmeasured case).
#'   The true exposure effect is stored in `attr(, "beta_x")`.
#' @export
#' @examples
#' df <- sim_confounding(n = 1000, beta_x = 0.5, seed = 42)
#' # Naive (biased) versus full (unbiased) model:
#' coef(lm(y ~ x + c1, df))["x"]
#' coef(lm(y ~ x + c1 + u, df))["x"]
sim_confounding <- function(n = 1000L, beta_x = 0,
                            family = c("gaussian", "binomial"),
                            seed = NULL) {
  family <- match.arg(family)
  if (!is.null(seed)) set.seed(seed)
  c1 <- stats::rnorm(n)
  u <- 0.5 * c1 + stats::rnorm(n)
  x <- 0.4 * c1 + 0.6 * u + stats::rnorm(n)
  lp <- beta_x * x + 0.5 * c1 + 0.8 * u
  y <- if (family == "gaussian") {
    lp + stats::rnorm(n)
  } else {
    stats::rbinom(n, 1L, stats::plogis(lp))
  }
  out <- data.frame(y = y, x = x, c1 = c1, u = u)
  attr(out, "beta_x") <- beta_x
  out
}
