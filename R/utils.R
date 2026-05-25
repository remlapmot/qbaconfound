# Internal helpers (not exported) -------------------------------------------

# Validate a length-2 numeric prior specification.
check_pair <- function(x, name) {
  x <- as.numeric(x)
  if (length(x) != 2L || anyNA(x)) {
    stop(sprintf("`%s` must be a length-2 numeric vector.", name),
         call. = FALSE)
  }
  x
}

# Normalise a `coef_exp` specification to a matrix with one row c(mean, sd)
# per exposure term.
normalise_coef_exp <- function(coef_exp) {
  if (is.list(coef_exp)) {
    coef_exp <- do.call(rbind, lapply(coef_exp, function(p) {
      p <- as.numeric(p)
      if (length(p) != 2L) {
        stop("Each `coef_exp` entry must be a length-2 numeric vector.",
             call. = FALSE)
      }
      p
    }))
  } else if (is.matrix(coef_exp)) {
    if (ncol(coef_exp) != 2L) {
      stop("A `coef_exp` matrix must have 2 columns (mean, sd).", call. = FALSE)
    }
    storage.mode(coef_exp) <- "double"
  } else {
    coef_exp <- matrix(as.numeric(coef_exp), ncol = 2L)
  }
  colnames(coef_exp) <- c("mean", "sd")
  coef_exp
}

# Resolve a `family` argument to a stats family object, mirroring glm().
resolve_family <- function(family) {
  if (is.character(family)) {
    family <- get(family, mode = "function", envir = parent.frame())
  }
  if (is.function(family)) {
    family <- family()
  }
  if (is.null(family$family)) {
    stop("`family` not recognised.", call. = FALSE)
  }
  family
}

# Is the left-hand side of `formula` a Surv() call? Handles both `Surv(...)`
# and the namespace-qualified `survival::Surv(...)`.
is_surv_lhs <- function(formula) {
  lhs <- formula[[2L]]
  if (!is.call(lhs)) return(FALSE)
  fn <- lhs[[1L]]
  nm <- if (is.call(fn)) as.character(fn[[length(fn)]]) else as.character(fn)
  identical(nm, "Surv")
}

# Build the exposure design (model-matrix columns belonging to `exposure`)
# and return the columns together with their coefficient names.
exposure_design <- function(formula, data, exposure) {
  tt <- stats::terms(formula, data = data)
  labs <- attr(tt, "term.labels")
  if (!all(exposure %in% labs)) {
    missing <- setdiff(exposure, labs)
    stop(sprintf("Exposure term(s) not found on the right-hand side: %s",
                 paste(missing, collapse = ", ")), call. = FALSE)
  }
  mf <- stats::model.frame(tt, data = data)
  mm <- stats::model.matrix(tt, mf)
  assign <- attr(mm, "assign")
  cols <- which(assign %in% match(exposure, labs))
  list(X = mm[, cols, drop = FALSE], names = colnames(mm)[cols])
}

# Solve for the logistic intercept a0 so that the marginal prevalence
# mean(plogis(a0 + lp)) equals the target prevalence `pi`.
solve_intercept <- function(pi, lp) {
  if (pi <= 0) return(-Inf)
  if (pi >= 1) return(Inf)
  f <- function(a0) mean(stats::plogis(a0 + lp)) - pi
  stats::uniroot(f, interval = c(-10, 10), extendInt = "upX",
                 tol = .Machine$double.eps^0.5)$root
}

# Draw a residual standard deviation eta from its prior.
draw_eta <- function(u) {
  if (u$resid_dist == "uniform") {
    stats::runif(1L, u$resid[1L], u$resid[2L])
  } else {
    sqrt(1 / stats::rgamma(1L, shape = u$resid[1L], scale = u$resid[2L]))
  }
}

# Draw a prevalence pi from its prior.
draw_pi <- function(u) {
  if (u$prev_dist == "uniform") {
    stats::runif(1L, u$prev[1L], u$prev[2L])
  } else {
    stats::rbeta(1L, u$prev[1L], u$prev[2L])
  }
}

# Simulate one proxy unmeasured confounder given drawn bias parameters and
# the exposure linear predictor `lp` = Xexp %*% alphaX.
simulate_proxy <- function(u, lp, n) {
  if (u$type == "continuous") {
    eta <- draw_eta(u)
    list(value = lp + eta * stats::rnorm(n), nuisance = c(eta = eta))
  } else {
    pi <- draw_pi(u)
    a0 <- solve_intercept(pi, lp)
    pr <- stats::plogis(a0 + lp)
    list(value = stats::rbinom(n, 1L, pr), nuisance = c(pi = pi, a0 = a0))
  }
}

# Fit the outcome model with the proxy contribution as a fixed offset.
# The offset is added to `data` as a column and referenced through an offset()
# term in the formula. This works uniformly for glm and coxph and avoids glm
# resolving a bare `offset` symbol in the formula's environment.
fit_outcome <- function(formula, data, family, offset, is_cox) {
  data[[".qba_offset"]] <- offset
  f <- stats::update(formula, . ~ . + offset(.qba_offset))
  if (is_cox) {
    survival::coxph(f, data = data)
  } else {
    stats::glm(f, data = data, family = family)
  }
}
