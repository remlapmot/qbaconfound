test_that("confounder constructors validate and normalise inputs", {
  u <- u_continuous(coef_out = c(0.8, 0.1), coef_exp = c(0.3, 0.05),
                    resid_sd = c(0.9, 1.1))
  expect_s3_class(u, "qba_u")
  expect_equal(dim(u$coef_exp), c(1L, 2L))

  # A list/matrix coef_exp expands to one row per exposure term.
  u2 <- u_binary(coef_out = c(0.5, 0.1),
                 coef_exp = list(c(0.3, 0.05), c(-0.2, 0.05)),
                 prevalence = c(0.1, 0.3))
  expect_equal(nrow(u2$coef_exp), 2L)

  expect_error(u_continuous(c(1), c(0.3, 0.05), c(0.9, 1.1)), "length-2")
  expect_error(u_binary(c(0.5, 0.1), c(0.3, 0.05), c(0.1, 1.5)), "\\[0, 1\\]")
})

test_that("solve_intercept reproduces the target prevalence", {
  set.seed(1)
  lp <- rnorm(5000, sd = 0.7)
  for (pi in c(0.1, 0.3, 0.6)) {
    a0 <- qbaconfound:::solve_intercept(pi, lp)
    expect_equal(mean(plogis(a0 + lp)), pi, tolerance = 1e-6)
  }
})

test_that("a near-zero-uncertainty prior at the truth recovers the full model (gaussian)", {
  df <- sim_confounding(n = 2000, beta_x = 0.5, seed = 3)
  full <- coef(lm(y ~ x + c1 + u, df))["x"]
  naive <- coef(lm(y ~ x + c1, df))["x"]
  expect_gt(abs(naive - 0.5), abs(full - 0.5))  # naive is more biased

  # Point-mass priors at values close to the data-generating mechanism.
  fit <- qbaconfound(
    y ~ x + c1, data = df, exposure = "x",
    confounders = u_continuous(coef_out = c(0.8, 0),
                               coef_exp = c(0.6, 0),
                               resid_sd = c(1, 1)),
    reps = 300, sampling_error = FALSE, seed = 7
  )
  est <- fit$estimates$estimate
  expect_lt(abs(est - 0.5), abs(naive - 0.5))   # adjustment reduces bias
})

test_that("qbaconfound runs for a binary outcome (glm) and returns structure", {
  df <- sim_confounding(n = 1000, beta_x = 0, family = "binomial", seed = 5)
  fit <- qbaconfound(
    y ~ x + c1, data = df, exposure = "x", family = binomial(),
    confounders = u_continuous(coef_out = c(0.8, 0.1),
                               coef_exp = c(0.6, 0.05),
                               resid_sd = c(0.9, 1.1)),
    reps = 100, seed = 11
  )
  expect_s3_class(fit, "qbaconfound")
  expect_equal(nrow(fit$estimates), 1L)
  expect_true(all(c("estimate", "conf.low", "conf.high") %in%
                    names(fit$estimates)))
  expect_lt(fit$estimates$conf.low, fit$estimates$conf.high)
})

test_that("multiple confounders and binary U are supported", {
  df <- sim_confounding(n = 800, beta_x = 0.3, seed = 9)
  fit <- qbaconfound(
    y ~ x + c1, data = df, exposure = "x",
    confounders = list(
      u_continuous(coef_out = c(0.8, 0.1), coef_exp = c(0.6, 0.05),
                   resid_sd = c(0.9, 1.1)),
      u_binary(coef_out = c(0.4, 0.1), coef_exp = c(0.3, 0.05),
               prevalence = c(0.15, 0.25))
    ),
    reps = 100, seed = 13
  )
  expect_equal(fit$n_failed, 0L)
  expect_length(fit$confounders, 2L)
})

test_that("Cox proportional hazards outcomes are detected and fitted", {
  set.seed(21)
  n <- 600
  c1 <- rnorm(n); u <- 0.5 * c1 + rnorm(n)
  x <- 0.4 * c1 + 0.6 * u + rnorm(n)
  time <- rexp(n, rate = exp(0.3 * x + 0.5 * c1 + 0.8 * u))
  status <- rbinom(n, 1L, 0.8)
  df <- data.frame(time, status, x, c1)

  fit <- qbaconfound(
    survival::Surv(time, status) ~ x + c1, data = df, exposure = "x",
    confounders = u_continuous(coef_out = c(0.8, 0.1),
                               coef_exp = c(0.6, 0.05),
                               resid_sd = c(0.9, 1.1)),
    reps = 60, seed = 23
  )
  expect_equal(fit$family, "cox")
  expect_equal(fit$n_failed, 0L)
})
