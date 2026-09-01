## B47 -- the pooled-inference contract for every consolidated family.
##
## skripte/k08_ueberdeckung_alle.R measures coverage of nominal 95 % intervals
## for all twelve representative cells at R = 1000 (one-level) / 300 (two-level).
## That takes ~20 minutes, so it cannot live in R CMD check. What CAN live here
## are the properties a coverage run depends on -- if any of them breaks, the
## coverage number becomes meaningless, and it breaks silently.
##
## The three properties, in the order a coverage run needs them:
##   (1) the method returns exactly one finite draw per missing cell, so the
##       completed data set is well defined;
##   (2) the between-imputation variance is strictly positive -- a method that
##       returns the same values every call would pool to zero B and produce
##       intervals that are too narrow no matter how correct beta is;
##   (3) the drawn values carry the predictor signal, so the analysis model
##       fitted on completed data can recover the slope at all.
##
## Property (2) is what B09 (fixed theta) and B40 (unscaled covariance) both
## damaged, in different ways: coverage failures in this package have so far
## always been failures of the between-imputation component, never of the point
## estimate. Asserting B > 0 per family is therefore the cheapest guard that
## would have caught them.

zi_type_1l <- c(x1 = 2L, x2 = 2L)   # both predictors in the count part

b47_daten_1l <- function(prozess, n = 250L, seed = 47L) {
  set.seed(seed)
  x1 <- stats::rnorm(n); x2 <- stats::rnorm(n)
  mu <- exp(1.0 + 0.6 * x1 - 0.4 * x2)
  y <- switch(prozess,
    poisson = stats::rpois(n, mu),
    nb      = MASS::rnegbin(n, mu = mu, theta = 1.5),
    zip     = ifelse(stats::rbinom(n, 1L, 0.25) == 1L, 0L, stats::rpois(n, mu)),
    hnb     = { z <- stats::rbinom(n, 1L, 0.65)
                ifelse(z == 0L, 0L, ci(".countimp_rktnb")(n, size = 1.5, k = 0,
                                                          mu = mu)) })
  ry <- stats::rbinom(n, 1L, stats::plogis(-0.8 + 0.7 * x2)) == 0L
  list(y = y, ry = ry, x = data.frame(x1 = x1, x2 = x2))
}

## M draws from one method on one data set.
b47_draws <- function(methode, d, M = 5L, type = NULL) {
  f <- get(paste0("mice.impute.", methode), envir = ci_home())
  arg <- list(y = d$y, ry = d$ry, x = d$x)
  if ("type" %in% names(formals(f))) arg$type <- type
  out <- matrix(NA_real_, sum(!d$ry), M)
  for (k in seq_len(M)) out[, k] <- suppressWarnings(do.call(f, arg))
  out
}

testthat::test_that("B47: every family returns one finite draw per missing cell", {
  ## pscl is in Suggests, not in Imports: the one-level zero-inflation and
  ## hurdle methods need it, the rest of the package does not. Without this
  ## skip the test fails in an environment without Suggests -- and it fails on
  ## the missing dependency, not on the property under test.
  skip_if_not_installed("pscl")
  faelle <- list(c("poisson", "poisson"), c("poisson.boot", "poisson"),
                 c("nb", "nb"), c("quasipoisson", "poisson"),
                 c("zip", "zip"), c("zinb", "zip"),
                 c("hp", "hnb"), c("hnb", "hnb"))
  for (f in faelle) {
    d <- b47_daten_1l(f[2])
    z <- b47_draws(f[1], d, M = 3L, type = zi_type_1l)
    testthat::expect_equal(nrow(z), sum(!d$ry), info = f[1])
    testthat::expect_true(all(is.finite(z)), info = f[1])
    testthat::expect_true(all(z >= 0), info = f[1])
    testthat::expect_true(all(z == round(z)), info = paste(f[1], "counts"))
  }
})

testthat::test_that("B47: between-imputation variance is strictly positive", {
  ## pscl is in Suggests, not in Imports: the one-level zero-inflation and
  ## hurdle methods need it, the rest of the package does not. Without this
  ## skip the test fails in an environment without Suggests -- and it fails on
  ## the missing dependency, not on the property under test.
  skip_if_not_installed("pscl")
  ## A method whose draws do not vary across m pools to B = 0 and yields
  ## intervals that are too narrow regardless of the point estimate. This is
  ## the failure mode of B09 and B40; it is invisible in a single-draw test.
  for (f in list(c("poisson", "poisson"), c("poisson.boot", "poisson"),
                 c("nb", "nb"), c("quasipoisson", "poisson"),
                 c("zip", "zip"), c("hnb", "hnb"))) {
    d <- b47_daten_1l(f[2])
    z <- b47_draws(f[1], d, M = 5L, type = zi_type_1l)
    ## mean of each completed data set must differ across m
    mittel <- colMeans(z)
    testthat::expect_gt(stats::var(mittel), 0, label = paste("B(", f[1], ")"))
    ## and the variation must not be a single outlying draw
    testthat::expect_gt(length(unique(round(mittel, 6))), 2L,
                        label = paste("distinct draws", f[1]))
  }
})

testthat::test_that("B47: drawn values carry the predictor signal", {
  ## pscl is in Suggests, not in Imports: the one-level zero-inflation and
  ## hurdle methods need it, the rest of the package does not. Without this
  ## skip the test fails in an environment without Suggests -- and it fails on
  ## the missing dependency, not on the property under test.
  skip_if_not_installed("pscl")
  ## If the draws were unconditional the analysis model could not recover any
  ## slope, and coverage would be an artefact of wide intervals. A positive
  ## correlation between the drawn count and the true linear predictor is the
  ## weakest form of "the imputation model conditions on x".
  for (f in list(c("poisson", "poisson"), c("nb", "nb"), c("zip", "zip"))) {
    d <- b47_daten_1l(f[2])
    z <- b47_draws(f[1], d, M = 3L, type = zi_type_1l)
    eta <- 1.0 + 0.6 * d$x[!d$ry, "x1"] - 0.4 * d$x[!d$ry, "x2"]
    r <- stats::cor(rowMeans(z), eta, method = "spearman")
    testthat::expect_gt(r, 0.3, label = paste("cor(draw, eta)", f[1]))
  }
})

testthat::test_that("B47: pooled interval covers the true slope on a fixed seed", {
  ## One seeded end-to-end pool per one-level family: not a coverage estimate
  ## (that is k08), but a guard that the whole chain -- draw, complete, refit,
  ## Rubin -- produces a usable interval rather than NaN or a degenerate width.
  M <- 10L
  for (f in list(c("poisson", "poisson"), c("nb", "nb"))) {
    d <- b47_daten_1l(f[2], seed = 471L)
    z <- b47_draws(f[1], d, M = M, type = zi_type_1l)
    B <- SE <- matrix(NA_real_, M, 3L)
    for (k in seq_len(M)) {
      dk <- data.frame(y = d$y, d$x); dk$y[!d$ry] <- z[, k]
      fk <- stats::glm(y ~ x1 + x2, family = stats::poisson, data = dk)
      B[k, ]  <- stats::coef(fk)
      SE[k, ] <- sqrt(diag(stats::vcov(fk)))
    }
    qbar <- colMeans(B); ubar <- colMeans(SE^2)
    bvar <- apply(B, 2, stats::var)
    tot  <- ubar + (1 + 1 / M) * bvar
    testthat::expect_true(all(is.finite(tot)), info = f[1])
    testthat::expect_true(all(tot > ubar), info = paste(f[1], "B contributes"))
    lam <- (1 + 1 / M) * bvar / tot
    testthat::expect_true(all(lam > 0 & lam < 1), info = paste(f[1], "lambda"))
    ## the x1 slope: true value 0.6 must lie in the pooled interval
    nu <- (M - 1) / pmax(lam^2, 1e-12)
    crit <- stats::qt(0.975, df = nu)
    lo <- qbar - crit * sqrt(tot); hi <- qbar + crit * sqrt(tot)
    testthat::expect_true(lo[2] <= 0.6 && 0.6 <= hi[2],
                          info = sprintf("%s: [%.3f, %.3f]", f[1], lo[2], hi[2]))
  }
})
