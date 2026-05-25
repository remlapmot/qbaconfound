#' @export
print.qbaconfound <- function(x, digits = 3L, ...) {
  cat("Monte Carlo QBA for unmeasured confounding\n")
  cat(sprintf("Outcome family: %s | replications: %d | observations: %d\n",
              x$family, x$reps, x$n_obs))
  cat(sprintf("Unmeasured confounder(s): %d (%s)\n",
              length(x$confounders),
              paste(vapply(x$confounders, `[[`, character(1L), "type"),
                    collapse = ", ")))
  if (x$n_failed > 0L) {
    cat(sprintf("Note: %d of %d replications failed to fit.\n",
                x$n_failed, x$reps))
  }
  cat("\n")

  est <- x$estimates
  out <- data.frame(
    term = est$term,
    naive = round(est$naive, digits),
    `bias-adjusted` = round(est$estimate, digits),
    `95% CI` = sprintf("(%s, %s)",
                       round(est$conf.low, digits),
                       round(est$conf.high, digits)),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  print(out, row.names = FALSE)
  invisible(x)
}

#' Summarise a Monte Carlo QBA
#'
#' @param object A `qbaconfound` object returned by [qbaconfound()].
#' @param ... Unused.
#' @return A data frame of the naive and bias-adjusted estimates for each
#'   exposure term, with columns `term`, `naive`, `naive_se`, `estimate`,
#'   `conf.low`, and `conf.high`.
#' @export
summary.qbaconfound <- function(object, ...) {
  object$estimates
}
