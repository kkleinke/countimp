## Shared data generators.
##
## Every generator is seeded so that a failure is reproducible from the test
## name alone. Sizes are deliberately small: the suite runs inside R CMD check,
## where glmmTMB fits dominate the runtime, so tests assert on structure and on
## invariants that hold at n = 150-400 rather than on asymptotic precision.

sim_count <- function(n = 300, nmis = 80, theta = 2, seed = 11) {
  set.seed(seed)
  x <- stats::rnorm(n)
  y <- MASS::rnegbin(n, mu = exp(1.2 + 0.6 * x), theta = theta)
  d <- data.frame(y = y, x = x)
  d$y[sample.int(n, nmis)] <- NA
  d
}

sim_zi_2l <- function(n = 300, ngrp = 30, nmis = 60, seed = 11) {
  set.seed(seed)
  id <- rep(seq_len(ngrp), each = n / ngrp)
  d <- data.frame(id = id, x1 = stats::rnorm(n), z1 = stats::rnorm(n))
  mu <- exp(-0.3 + 0.5 * d$x1 + rep(stats::rnorm(ngrp, 0, 0.3), each = n / ngrp))
  d$y <- MASS::rnegbin(n, mu, theta = 2)
  d$y[stats::rbinom(n, 1, stats::plogis(-0.5 + 0.6 * d$z1)) == 1] <- 0
  d$y[sample.int(n, nmis)] <- NA
  d
}

## Underdispersed counts: binomial draws have variance < mean, which is the
## configuration that produced NaN imputations before B04 was fixed.
sim_underdispersed <- function(n = 300, nmis = 90, seed = 3) {
  set.seed(seed)
  x <- stats::rnorm(n)
  d <- data.frame(y = stats::rbinom(n, size = 6, prob = stats::plogis(0.3 + 0.5 * x)),
                  x = x)
  d$y[sample.int(n, nmis)] <- NA
  d
}

## Two-level negative binomial fit, used by the tests that inspect internals
## directly rather than going through countimp().
fit_2l_nb <- function(n = 300, ngrp = 15, seed = 11) {
  set.seed(seed)
  g <- rep(seq_len(ngrp), each = n / ngrp)
  x <- stats::rnorm(n)
  y <- MASS::rnegbin(n, mu = exp(1 + 0.5 * x + rep(stats::rnorm(ngrp), each = n / ngrp)),
                     theta = 2)
  glmmTMB::glmmTMB(y ~ x + (1 | g), data = data.frame(y, x, g),
                   family = glmmTMB::nbinom2)
}

## Collect imputed values as a plain vector, without requiring mice.
imputed_values <- function(imp, var = "y") as.numeric(unlist(imp$imp[[var]]))

quiet_impute <- function(...) suppressWarnings(countimp(..., printFlag = FALSE))

## Two-level count data plus the fit on it. The tests for the joint draw need
## both: the data to build newdata from, the fit to draw from. `g` is a factor
## because .countimp_re_design() matches against levels(flist[[gf]]).
sim_count_2l <- function(n = 300, ngrp = 15, theta = 2, seed = 11) {
  set.seed(seed)
  g <- rep(seq_len(ngrp), each = n / ngrp)
  x <- stats::rnorm(n)
  y <- MASS::rnegbin(n, mu = exp(1 + 0.5 * x + rep(stats::rnorm(ngrp), each = n / ngrp)),
                     theta = theta)
  data.frame(y = y, x = x, g = factor(g))
}

## Fit on data produced by sim_count_2l(). Returns NULL when glmmTMB is absent
## or the fit fails, so tests can skip cleanly.
fit_2l_on <- function(d, form = y ~ x + (1 | g)) {
  if (!requireNamespace("glmmTMB", quietly = TRUE)) return(NULL)
  tryCatch(suppressWarnings(glmmTMB::glmmTMB(form, data = d,
           family = glmmTMB::nbinom2)), error = function(e) NULL)
}
