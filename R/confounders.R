#' Specify a continuous unmeasured confounder
#'
#' Describes the prior distributions of the bias parameters for a single
#' continuous unmeasured confounder, for use in [qbaconfound()].
#'
#' The bias model relates the unmeasured confounder `U` to the study data
#' through three bias parameters: the coefficient of `U` in the outcome model
#' (`coef_out`, i.e. \eqn{\beta_U}), the coefficient(s) of the exposure in the
#' model for `U` (`coef_exp`, i.e. \eqn{\alpha_X}), and the residual standard
#' deviation of `U` given the exposure (`resid_sd`, i.e. \eqn{\eta}). Values
#' for these parameters cannot be estimated from the data and so are drawn from
#' the prior distributions specified here.
#'
#' @param coef_out Length-2 numeric `c(mean, sd)` giving the normal prior for
#'   the coefficient of the unmeasured confounder in the outcome model
#'   (\eqn{\beta_U}). Set `sd = 0` for a point-mass (fixed value) prior.
#' @param coef_exp Normal prior for the coefficient(s) of the exposure in the
#'   model for the unmeasured confounder (\eqn{\alpha_X}). For a single
#'   exposure term, a length-2 numeric `c(mean, sd)`. For a categorical
#'   exposure with several terms, a matrix with one row `c(mean, sd)` per
#'   exposure term, or a list of such length-2 vectors.
#' @param resid_sd Length-2 numeric giving the prior for the residual standard
#'   deviation \eqn{\eta}. If `resid_dist = "uniform"` (the default) this is
#'   `c(min, max)` of a uniform prior on \eqn{\eta}. If `resid_dist = "gamma"`
#'   this is `c(shape, scale)` of a gamma prior on the residual *precision*
#'   \eqn{1/\eta^2}.
#' @param resid_dist Prior family for the residual standard deviation: either
#'   `"uniform"` (on \eqn{\eta}) or `"gamma"` (on the precision \eqn{1/\eta^2}).
#'
#' @return An object of class `qba_u` describing a continuous unmeasured
#'   confounder.
#' @seealso [u_binary()], [qbaconfound()]
#' @export
#' @examples
#' u_continuous(coef_out = c(0.8, 0.1), coef_exp = c(0.3, 0.05),
#'              resid_sd = c(0.9, 1.1))
u_continuous <- function(coef_out, coef_exp, resid_sd,
                         resid_dist = c("uniform", "gamma")) {
  resid_dist <- match.arg(resid_dist)
  coef_out <- check_pair(coef_out, "coef_out")
  coef_exp <- normalise_coef_exp(coef_exp)
  resid_sd <- check_pair(resid_sd, "resid_sd")
  if (resid_dist == "uniform" && resid_sd[1] < 0) {
    stop("For a uniform prior, `resid_sd` lower limit must be non-negative.",
         call. = FALSE)
  }
  structure(
    list(type = "continuous",
         coef_out = coef_out,
         coef_exp = coef_exp,
         resid = resid_sd,
         resid_dist = resid_dist),
    class = c("qba_u_continuous", "qba_u")
  )
}

#' Specify a binary unmeasured confounder
#'
#' Describes the prior distributions of the bias parameters for a single binary
#' unmeasured confounder, for use in [qbaconfound()].
#'
#' As for [u_continuous()] the coefficient of `U` in the outcome model
#' (`coef_out`) and the coefficient(s) of the exposure in the model for `U`
#' (`coef_exp`) are bias parameters. For a binary confounder the third bias
#' parameter is the marginal prevalence of `U` (`prevalence`, i.e.
#' \eqn{\pi}) rather than a residual standard deviation. A prevalence is
#' usually easier to elicit and more readily reported in the literature than
#' the intercept of a logistic model; the intercept needed to reproduce the
#' drawn prevalence is derived internally.
#'
#' @inheritParams u_continuous
#' @param prevalence Length-2 numeric giving the prior for the marginal
#'   prevalence \eqn{\pi} of the binary unmeasured confounder. If
#'   `prev_dist = "uniform"` (the default) this is `c(min, max)` of a uniform
#'   prior; if `prev_dist = "beta"` this is `c(a, b)` of a beta prior.
#' @param prev_dist Prior family for the prevalence: `"uniform"` or `"beta"`.
#'
#' @return An object of class `qba_u` describing a binary unmeasured
#'   confounder.
#' @seealso [u_continuous()], [qbaconfound()]
#' @export
#' @examples
#' u_binary(coef_out = c(0.7, 0.1), coef_exp = c(0.4, 0.05),
#'          prevalence = c(0.15, 0.25))
u_binary <- function(coef_out, coef_exp, prevalence,
                     prev_dist = c("uniform", "beta")) {
  prev_dist <- match.arg(prev_dist)
  coef_out <- check_pair(coef_out, "coef_out")
  coef_exp <- normalise_coef_exp(coef_exp)
  prevalence <- check_pair(prevalence, "prevalence")
  if (prev_dist == "uniform" &&
      (prevalence[1] < 0 || prevalence[2] > 1)) {
    stop("For a uniform prior, `prevalence` must lie within [0, 1].",
         call. = FALSE)
  }
  structure(
    list(type = "binary",
         coef_out = coef_out,
         coef_exp = coef_exp,
         prev = prevalence,
         prev_dist = prev_dist),
    class = c("qba_u_binary", "qba_u")
  )
}

#' @export
print.qba_u <- function(x, ...) {
  cat(sprintf("<%s unmeasured confounder>\n", x$type))
  cat(sprintf("  coef_out (beta_U): N(mean = %g, sd = %g)\n",
              x$coef_out[1], x$coef_out[2]))
  ce <- x$coef_exp
  for (i in seq_len(nrow(ce))) {
    cat(sprintf("  coef_exp (alpha_X)[%d]: N(mean = %g, sd = %g)\n",
                i, ce[i, 1], ce[i, 2]))
  }
  if (x$type == "continuous") {
    if (x$resid_dist == "uniform") {
      cat(sprintf("  resid_sd (eta): Uniform(%g, %g)\n",
                  x$resid[1], x$resid[2]))
    } else {
      cat(sprintf("  1/eta^2: Gamma(shape = %g, scale = %g)\n",
                  x$resid[1], x$resid[2]))
    }
  } else {
    if (x$prev_dist == "uniform") {
      cat(sprintf("  prevalence (pi): Uniform(%g, %g)\n",
                  x$prev[1], x$prev[2]))
    } else {
      cat(sprintf("  prevalence (pi): Beta(%g, %g)\n",
                  x$prev[1], x$prev[2]))
    }
  }
  invisible(x)
}
