#' Monte Carlo quantitative bias analysis for unmeasured confounding
#'
#' Conducts the flexible Monte Carlo quantitative bias analysis (QBA) for
#' unmeasured confounding of Hughes et al. Given a naive analysis model (which
#' omits one or more unmeasured confounders) and informative priors for a small
#' number of bias parameters, the function returns a bias-adjusted estimate of
#' the exposure effect together with an interval that accounts for both the
#' unmeasured confounding and sampling variability.
#'
#' @details
#' The substantive analysis may be a generalised linear model (any [stats::glm()]
#' `family`) or a Cox proportional hazards model. Survival outcomes are detected
#' automatically when the left-hand side of `formula` is a [survival::Surv()]
#' call, or `family` can be set to `"cox"`.
#'
#' Each Monte Carlo replication (see Hughes et al., section 2.4):
#' \enumerate{
#'   \item draws a value for every bias parameter from its prior (the priors are
#'     specified through [u_continuous()] / [u_binary()]);
#'   \item simulates a proxy for each unmeasured confounder as a function of the
#'     exposure only (a continuous proxy from a normal model, a binary proxy from
#'     a Bernoulli model whose intercept reproduces the drawn prevalence);
#'   \item refits the outcome model including the simulated proxies, with their
#'     coefficients fixed to the drawn values (implemented as a model offset),
#'     and reads off the exposure coefficient and its standard error;
#'   \item adds Monte Carlo sampling error by drawing the bias-adjusted estimate
#'     from a normal distribution centred on that coefficient.
#' }
#' The point estimate is the median and the interval the 2.5th and 97.5th
#' percentiles of the resulting distribution of bias-adjusted estimates.
#'
#' @param formula The naive analysis model, e.g. `y ~ x + c1 + c2` for a GLM or
#'   `Surv(time, status) ~ x + c1` for a Cox model. The unmeasured confounders
#'   are *not* included in `formula`.
#' @param data A data frame containing the outcome, exposure, and measured
#'   confounders.
#' @param exposure Character vector naming the exposure term(s) in `formula`
#'   whose effect is of interest. Defaults to the first term on the right-hand
#'   side.
#' @param confounders A single [u_continuous()]/[u_binary()] object, or a list
#'   of them, describing the unmeasured confounder(s).
#' @param family The outcome model family: a [stats::glm()] family object,
#'   a family name, or `"cox"` for a Cox proportional hazards model. Ignored
#'   (and inferred as Cox) when the response is a [survival::Surv()] call.
#' @param reps Number of Monte Carlo replications.
#' @param sampling_error Logical; if `TRUE` (the default) Monte Carlo sampling
#'   error is incorporated at step 4 above. Set to `FALSE` to obtain the
#'   distribution of bias-adjusted point estimates without sampling error.
#' @param seed Optional integer seed for reproducibility.
#'
#' @return An object of class `qbaconfound`, a list with elements including
#'   `estimates` (a data frame of naive and bias-adjusted estimates per exposure
#'   term), `draws` (the matrix of bias-adjusted estimates across replications),
#'   and `n_failed` (the number of replications whose model fit failed).
#'
#' @references
#' Hughes RA, Kawabata E, Palmer TM, et al. A flexible Monte Carlo quantitative
#' bias analysis for unmeasured confounding. *Statistical Methods in Medical
#' Research* (forthcoming).
#'
#' @seealso [u_continuous()], [u_binary()], [sim_confounding()]
#' @export
#' @examples
#' df <- sim_confounding(n = 500, beta_x = 0, seed = 1)
#'
#' # Naive model y ~ x + c1 omits the confounder u; adjust for one continuous U.
#' fit <- qbaconfound(
#'   y ~ x + c1, data = df, exposure = "x",
#'   confounders = u_continuous(coef_out = c(0.8, 0.1),
#'                              coef_exp = c(0.6, 0.05),
#'                              resid_sd = c(0.9, 1.1)),
#'   reps = 200, seed = 1
#' )
#' fit
qbaconfound <- function(formula, data, exposure = NULL, confounders,
                        family = stats::gaussian(), reps = 1000L,
                        sampling_error = TRUE, seed = NULL) {
  if (!inherits(formula, "formula") || length(formula) != 3L) {
    stop("`formula` must be a two-sided formula, e.g. y ~ x + c.",
         call. = FALSE)
  }
  if (missing(confounders)) {
    stop("`confounders` must be supplied (see u_continuous(), u_binary()).",
         call. = FALSE)
  }
  if (inherits(confounders, "qba_u")) confounders <- list(confounders)
  if (!is.list(confounders) ||
      !all(vapply(confounders, inherits, logical(1L), "qba_u"))) {
    stop("`confounders` must be a u_continuous()/u_binary() object or a list ",
         "of them.", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)

  is_cox <- is_surv_lhs(formula) ||
    (is.character(family) && family[1L] %in% c("cox", "coxph"))
  if (!is_cox) family <- resolve_family(family)

  # Default exposure to the first right-hand-side term.
  rhs_terms <- attr(stats::terms(formula, data = data), "term.labels")
  if (is.null(exposure)) exposure <- rhs_terms[1L]

  # Work on complete cases for the model variables so offsets align with rows.
  model_vars <- all.vars(formula)
  data <- data[stats::complete.cases(data[model_vars]), , drop = FALSE]
  n <- nrow(data)
  if (n == 0L) stop("No complete cases for the model variables.", call. = FALSE)

  ed <- exposure_design(formula, data, exposure)
  Xexp <- ed$X
  exp_names <- ed$names
  n_terms <- length(exp_names)

  # Validate that each confounder's coef_exp matches the number of exposure
  # terms.
  for (k in seq_along(confounders)) {
    ce <- confounders[[k]]$coef_exp
    if (nrow(ce) != n_terms) {
      stop(sprintf(paste0("Confounder %d: `coef_exp` has %d row(s) but the ",
                          "exposure expands to %d term(s) (%s)."),
                   k, nrow(ce), n_terms,
                   paste(exp_names, collapse = ", ")), call. = FALSE)
    }
  }

  # Naive estimate (outcome model with no offset).
  naive_fit <- fit_outcome(formula, data, family, offset = rep(0, n), is_cox)
  naive_est <- stats::coef(naive_fit)[exp_names]
  naive_se <- sqrt(diag(stats::vcov(naive_fit))[exp_names])

  # Monte Carlo replications.
  draws <- matrix(NA_real_, nrow = reps, ncol = n_terms,
                  dimnames = list(NULL, exp_names))
  n_failed <- 0L
  for (m in seq_len(reps)) {
    offset <- numeric(n)
    for (u in confounders) {
      alpha_x <- stats::rnorm(n_terms, mean = u$coef_exp[, "mean"],
                              sd = u$coef_exp[, "sd"])
      lp <- as.numeric(Xexp %*% alpha_x)
      proxy <- simulate_proxy(u, lp, n)
      beta_u <- stats::rnorm(1L, u$coef_out[1L], u$coef_out[2L])
      offset <- offset + beta_u * proxy$value
    }

    fit <- tryCatch(
      fit_outcome(formula, data, family, offset = offset, is_cox),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      n_failed <- n_failed + 1L
      next
    }
    bhat <- stats::coef(fit)[exp_names]
    if (sampling_error) {
      se <- sqrt(diag(stats::vcov(fit))[exp_names])
      draws[m, ] <- stats::rnorm(n_terms, mean = bhat, sd = se)
    } else {
      draws[m, ] <- bhat
    }
  }

  qs <- apply(draws, 2L, stats::quantile,
              probs = c(0.5, 0.025, 0.975), na.rm = TRUE)
  estimates <- data.frame(
    term = exp_names,
    naive = unname(naive_est),
    naive_se = unname(naive_se),
    estimate = qs["50%", ],
    conf.low = qs["2.5%", ],
    conf.high = qs["97.5%", ],
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  structure(
    list(
      estimates = estimates,
      draws = draws,
      exposure = exp_names,
      confounders = confounders,
      family = if (is_cox) "cox" else family$family,
      reps = reps,
      n_failed = n_failed,
      n_obs = n,
      sampling_error = sampling_error,
      call = match.call()
    ),
    class = "qbaconfound"
  )
}
