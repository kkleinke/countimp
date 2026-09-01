## test-B69-ztrunc.R -- zero-truncated single-level count imputation
##
## Bugs found while building R/ztrunc.R, each pinned below:
##
## B69  The theta component of the ztnb gradient had the wrong sign on the
##      truncation term. Non-obvious failure mode: the beta components stayed
##      exact and the fitted coefficients looked right, so nothing downstream
##      complained -- only the finite-difference check caught it (off by a
##      factor of -90). A gradient error costs convergence, not correctness,
##      which is why it hides.
## B70  log(1 - exp(-mu)) computed naively loses all precision for small mu
##      (1 - exp(-1e-8) rounds to 0 in double) and underflows for large mu.
##      Small mu is not exotic: a Bayesian beta draw with an extreme predictor
##      produces it routinely.
## B71  A timing claim measured against a cold glmmTMB call suggested a 46x
##      speed advantage that was pure TMB template loading. Pinned here as a
##      warm-start comparison so the documented factors stay honest.
## B72  Imputing zero-truncated data with "poisson" writes structurally
##      impossible zeros (12.9% of imputations measured) and then the correct
##      analysis model cannot be fitted at all -- 200 of 200 runs aborted.

zt_dat <- function(seed = 7, n = 600, dist = "ztp", theta = 2,
                   b = c(0.8, 0.5, -0.3)) {
  set.seed(seed)
  x1 <- rnorm(n); x2 <- rnorm(n)
  X  <- cbind(1, x1, x2)
  mu <- exp(drop(X %*% b))
  y  <- if (dist == "ztp") stats::rpois(n, mu) else
        MASS::rnegbin(n, mu = mu, theta = theta)
  while (any(y == 0)) {
    i <- y == 0
    y[i] <- if (dist == "ztp") stats::rpois(sum(i), mu[i]) else
            MASS::rnegbin(sum(i), mu = mu[i], theta = theta)
  }
  list(y = y, X = X, x = data.frame(x1 = x1, x2 = x2), b = b, mu = mu)
}

## ---- B69: analytic gradients agree with finite differences ---------------
## The check that found the sign error. Central differences at a point away
## from the optimum, where the gradient is large and a wrong term shows.

test_that("B69: ztp gradient matches finite differences", {
  d  <- zt_dat(dist = "ztp")
  nll <- ci(".countimp_ztp_nll"); gr <- ci(".countimp_ztp_gr")
  b0 <- c(0.5, 0.3, -0.1)
  fd <- vapply(seq_along(b0), function(j) {
    e <- numeric(length(b0)); e[j] <- 1e-6
    (nll(b0 + e, d$X, d$y) - nll(b0 - e, d$X, d$y)) / 2e-6
  }, numeric(1))
  an <- gr(b0, d$X, d$y)
  expect_lt(max(abs(an - fd) / abs(fd)), 1e-6)
})

test_that("B69: ztnb gradient matches finite differences in every component", {
  d   <- zt_dat(dist = "ztnb")
  nll <- ci(".countimp_ztnb_nll"); gr <- ci(".countimp_ztnb_gr")
  p0  <- c(0.5, 0.3, -0.1, log(1.5))
  fd  <- vapply(seq_along(p0), function(j) {
    e <- numeric(length(p0)); e[j] <- 1e-6
    (nll(p0 + e, d$X, d$y) - nll(p0 - e, d$X, d$y)) / 2e-6
  }, numeric(1))
  an <- gr(p0, d$X, d$y)
  ## component-wise, not just the maximum: the bug left beta exact and only
  ## theta (the last component) wrong, so a norm over all four could pass
  for (j in seq_along(p0))
    expect_lt(abs(an[j] - fd[j]) / abs(fd[j]), 1e-6)
})

## ---- B70: log1mexp is accurate at both ends ------------------------------

test_that("B70: log1mexp is accurate where the naive form fails", {
  f <- ci(".countimp_log1mexp")
  ## Measured breakdown of log(1 - exp(-a)), not assumed: relative error
  ## 8e-7 at a = 1e-12, 2.8e-3 at a = 1e-16, and -Inf from a = 1e-17 down.
  ## The split at log 2 holds all of them to double precision.
  for (a in c(1e-12, 1e-16, 1e-17, 1e-20))
    expect_lt(abs(f(a) - (log(a) - a / 2)), 1e-9, label = paste("a =", a))
  ## the naive form: already 3 digits off at 1e-16, gone at 1e-17
  expect_gt(abs(log(1 - exp(-1e-16)) - log(1e-16)) / abs(log(1e-16)), 1e-3)
  expect_true(is.infinite(log(1 - exp(-1e-17))))
  ## upper end: exp(-a) underflows at a = 750, so the result is -0 exactly
  expect_equal(f(800), -exp(-800))
  expect_lt(abs(f(700) - (-exp(-700))), 1e-310)
  expect_true(all(is.finite(f(c(1e-12, 1e-3, 0.5, log(2), 1, 50, 700, 800)))))
})

## ---- fit accuracy against glmmTMB ---------------------------------------
## glmmTMB is the reference here, not a dependency: R/ztrunc.R must reproduce
## it without needing it at run time.

test_that("ztp fit reproduces glmmTMB::truncated_poisson", {
  skip_if_not_installed("glmmTMB")
  d <- zt_dat(dist = "ztp", n = 800)
  f <- ci(".countimp_zt_fit")(d$X, d$y, "ztp")
  g <- glmmTMB::glmmTMB(y ~ x1 + x2, cbind(y = d$y, d$x),
                        family = glmmTMB::truncated_poisson())
  expect_lt(max(abs(f$beta - glmmTMB::fixef(g)$cond)), 1e-4)
  expect_lt(max(abs(sqrt(diag(f$cov)) - sqrt(diag(vcov(g)$cond)))), 1e-5)
  expect_lt(abs(f$ll - as.numeric(stats::logLik(g))), 1e-5)
})

test_that("ztnb fit reproduces glmmTMB::truncated_nbinom2 including theta", {
  skip_if_not_installed("glmmTMB")
  d <- zt_dat(dist = "ztnb", n = 800, seed = 9)
  f <- ci(".countimp_zt_fit")(d$X, d$y, "ztnb")
  g <- glmmTMB::glmmTMB(y ~ x1 + x2, cbind(y = d$y, d$x),
                        family = glmmTMB::truncated_nbinom2())
  expect_lt(max(abs(f$beta - glmmTMB::fixef(g)$cond)), 1e-3)
  ## ci(), not a bare call: the bare name resolved to nothing under
  ## R CMD check (tests do not see package internals) and the test errored
  ## instead of failing -- B38 did not catch it because there is no
  ## countimp::: prefix to grep for.
  expect_lt(abs(f$theta - ci(".countimp_theta_num")(g)) / f$theta, 1e-3)
  expect_lt(abs(f$ll - as.numeric(stats::logLik(g))), 1e-4)
})

## ---- B72: imputations are never zero ------------------------------------
## The defining property. If this fails the method is pointless.

test_that("B72: no method in the ztp family ever imputes a zero", {
  d <- zt_dat(dist = "ztp", n = 500)
  set.seed(11)
  miss <- stats::runif(500) < stats::plogis(-0.4 + 0.9 * d$x$x2)
  ym <- d$y; ym[miss] <- NA
  for (meth in c("ztp", "ztp.boot", "ztnb", "ztnb.boot")) {
    im <- suppressWarnings(get(paste0("mice.impute.", meth))(ym, !miss, d$x))
    expect_true(all(im >= 1), info = meth)
    expect_true(all(is.finite(im)), info = meth)
    expect_equal(length(im), sum(miss), info = meth)
  }
})

test_that("B72: poisson imputation of the same data does produce zeros", {
  ## The contrast that motivates the method. Not a test of countimp's
  ## correctness -- a test that the problem being solved is real.
  d <- zt_dat(dist = "ztp", n = 500)
  set.seed(11)
  miss <- stats::runif(500) < stats::plogis(-0.4 + 0.9 * d$x$x2)
  ym <- d$y; ym[miss] <- NA
  ip <- suppressWarnings(mice.impute.poisson(ym, !miss, d$x))
  expect_gt(sum(ip == 0), 0)
})

## ---- observed zeros are an error, with a message that helps --------------

test_that("B73: observed zeros are refused with the cause, not a retry failure", {
  ## pscl is in Suggests, not in Imports: the one-level zero-inflation and
  ## hurdle methods need it, the rest of the package does not. Without this
  ## skip the test fails in an environment without Suggests -- and it fails on
  ## the missing dependency, not on the property under test.
  skip_if_not_installed("pscl")
  ## The check has to run before .countimp_draw_retry(). Raised from inside the
  ## draw closure it was retried three times and then reported as "failed on 3
  ## successive draws ... extreme predictor values, a separated model or too few
  ## observed cases" -- a convergence diagnosis for what is a wrong-method
  ## error, with the real cause in a parenthesis.
  d <- zt_dat(dist = "ztp", n = 200)
  y0 <- d$y; y0[c(3, 17)] <- 0
  ry <- rep(TRUE, 200); ry[150:200] <- FALSE
  for (m in c("ztp", "ztp.boot", "ztnb", "ztnb.boot")) {
    f <- get(paste0("mice.impute.", m))
    msg <- tryCatch({ f(y0, ry, d$x); "" }, error = conditionMessage)
    expect_match(msg, "zero-truncated", info = m)
    expect_match(msg, "2 observed value", info = m)
    ## names the alternatives that do model zeros
    expect_match(msg, "hp", info = m)
    expect_match(msg, "zip", info = m)
    ## and does NOT masquerade as a convergence problem
    expect_false(grepl("successive draws", msg), info = m)
    expect_false(grepl("collinear", msg), info = m)
  }
})

test_that("B73: the zero check counts only observed cases", {
  ## NA is not a zero: a variable whose gap is being imputed must not trip the
  ## check, or the method could never be used at all.
  d <- zt_dat(dist = "ztp", n = 200)
  ry <- rep(TRUE, 200); ry[150:200] <- FALSE
  ym <- d$y; ym[!ry] <- NA
  expect_silent(ci(".countimp_zt_check_zeros")(ym[ry], "ztp"))
  im <- suppressWarnings(mice.impute.ztp(ym, ry, d$x))
  expect_true(all(im >= 1))
})

## ---- recovery of the data-generating parameters -------------------------

test_that("ztp imputation recovers beta under MAR", {
  ## Measured: bias -0.0024, coverage 94.5% over 200 replications (R = 200,
  ## n = 300, M = 5). Tolerance set against that measurement -- loose enough
  ## not to flutter on one seed, tight enough to fail if the draw stops being
  ## zero-truncated or the fit loses the truncation term.
  b <- c(0.8, 0.5, -0.3); n <- 400; M <- 5
  est <- numeric(M)
  set.seed(4242)
  x1 <- rnorm(n); x2 <- rnorm(n); X <- cbind(1, x1, x2)
  mu <- exp(drop(X %*% b))
  y  <- stats::rpois(n, mu)
  while (any(y == 0)) { i <- y == 0; y[i] <- stats::rpois(sum(i), mu[i]) }
  miss <- stats::runif(n) < stats::plogis(-0.4 + 0.9 * x2)
  ym <- y; ym[miss] <- NA
  for (m in seq_len(M)) {
    im <- suppressWarnings(mice.impute.ztp(ym, !miss, data.frame(x1 = x1, x2 = x2)))
    yc <- ym; yc[miss] <- im
    ## the correct analysis model, fitted without repair: it only works
    ## because the imputations respect the truncation
    est[m] <- ci(".countimp_zt_fit")(X, yc, "ztp")$beta[2]
  }
  expect_lt(abs(mean(est) - b[2]), 0.12)
})

test_that("ztnb recovers theta and does not collapse to Poisson", {
  d <- zt_dat(dist = "ztnb", n = 900, seed = 21, theta = 1.5)
  f <- ci(".countimp_zt_fit")(d$X, d$y, "ztnb")
  expect_gt(f$theta, 0.8)
  expect_lt(f$theta, 3.5)
  ## and the ztp fit on the same data must be clearly worse
  fp <- ci(".countimp_zt_fit")(d$X, d$y, "ztp")
  expect_gt(f$ll, fp$ll + 5)
})

## ---- B71: the documented speed factors are warm-start factors -----------

test_that("B71: ztp fit is faster than glmmTMB with the template loaded", {
  skip_if_not_installed("glmmTMB")
  skip_on_cran()
  d <- zt_dat(dist = "ztp", n = 800)
  dd <- cbind(y = d$y, d$x)
  ## warm both sides; a cold glmmTMB call measures template loading
  invisible(glmmTMB::glmmTMB(y ~ x1, dd[1:100, ],
                            family = glmmTMB::truncated_poisson()))
  invisible(ci(".countimp_zt_fit")(d$X[1:100, ], d$y[1:100], "ztp"))
  t1 <- system.time(for (i in 1:5) ci(".countimp_zt_fit")(d$X, d$y, "ztp"))[["elapsed"]]
  t2 <- system.time(for (i in 1:5) glmmTMB::glmmTMB(y ~ x1 + x2, dd,
          family = glmmTMB::truncated_poisson()))[["elapsed"]]
  expect_lt(t1, t2)
})

## ---- input contract ------------------------------------------------------

test_that("ztrunc methods honour wy and the method contract", {
  d <- zt_dat(dist = "ztp", n = 300)
  ry <- rep(TRUE, 300); ry[200:300] <- FALSE
  ym <- d$y; ym[!ry] <- NA
  wy <- !ry; wy[250:300] <- FALSE          # impute only part of the gap
  im <- suppressWarnings(mice.impute.ztp(ym, ry, d$x, wy = wy))
  expect_equal(length(im), sum(wy))
  expect_true(all(im >= 1))
})

test_that("every ztrunc method is registered and callable", {
  for (m in c("ztp", "ztp.boot", "ztnb", "ztnb.boot")) {
    fn <- paste0("mice.impute.", m)
    expect_true(exists(fn, envir = asNamespace("countimp")), info = fn)
    expect_true(is.function(get(fn, envir = asNamespace("countimp"))), info = fn)
  }
})

## ---- integration with countimp_fit_diag() -------------------------------

test_that("B74: a converged truncated candidate is ranked, not listed as failed", {
  ## .countimp_fit_candidate() returned status "fitted" where the table tests
  ## for "ok", so ztp/ztnb appeared under "not ranked" with status "fitted" --
  ## a line contradicting itself. The status strings are a contract.
  d <- zt_dat(dist = "ztp", n = 400)
  r <- countimp_fit_diag(y ~ x1 + x2, cbind(y = d$y, d$x))
  for (f in c("ztp", "ztnb")) {
    row <- r$table[r$table$family == f, ]
    expect_equal(nrow(row), 1L, info = f)
    expect_true(row$fitted, info = f)
    expect_equal(row$status, "ok", info = f)
    expect_true(is.finite(row$aic), info = f)
    ## P(Y = 0) is 0 by construction, not NA
    expect_equal(row$p0_model, 0, info = f)
  }
})

test_that("truncated candidates win on truncated data", {
  d <- zt_dat(dist = "ztp", n = 500)
  r <- countimp_fit_diag(y ~ x1 + x2, cbind(y = d$y, d$x))
  expect_equal(r$recommendation, "ztp")
  ## and by a wide margin over the untruncated families -- measured at 123 AIC
  ## units on this design
  expect_gt(r$table$aic[r$table$family == "poisson"] -
            r$table$aic[r$table$family == "ztp"], 50)
})

test_that("truncated candidates are dropped when the data contain zeros", {
  set.seed(31)
  n <- 400; x1 <- rnorm(n); x2 <- rnorm(n)
  y <- stats::rpois(n, exp(0.3 + 0.5 * x1))
  stopifnot(any(y == 0))
  r <- countimp_fit_diag(y ~ x1 + x2, data.frame(y = y, x1 = x1, x2 = x2))
  expect_false(any(c("ztp", "ztnb") %in% r$table$family))
  ## dropped silently: not reported as a failure the user has to interpret
  expect_false(any(grepl("truncated", r$table$label)))
})

test_that("zero-free data get the substantive-choice note", {
  d <- zt_dat(dist = "ztp", n = 400)
  r <- countimp_fit_diag(y ~ x1 + x2, cbind(y = d$y, d$x))
  txt <- paste(utils::capture.output(print(r)), collapse = " ")
  expect_match(txt, "no zeros")
  expect_match(txt, "impossible by design")
  ## the note must say the AIC cannot settle it
  expect_match(txt, "not one the AIC can answer")
})

test_that("at a high rate Poisson stays in the tie band on zero-free data", {
  ## Measured: on Poisson samples that are zero-free of their own accord,
  ## Poisson lands in the tie band 15/15 times at mu = 6 and mu = 8 (median
  ## dAIC 1.03 and 0.15). At mu = 4 it does not, and that is correct -- see
  ## the selection-effect note in R/fitdiag.R.
  found <- FALSE
  for (s in 1:400) {
    set.seed(s * 7)
    x1 <- rnorm(200, 0, 0.25)
    y  <- stats::rpois(200, exp(log(8) + 0.2 * x1))
    if (any(y == 0)) next
    r <- countimp_fit_diag(y ~ x1, data.frame(y = y, x1 = x1))
    expect_true(any(grepl("^Poisson$", r$ties)))
    found <- TRUE
    break
  }
  expect_true(found)
})

test_that("logLik/AIC/npar work on the truncated fit through the generic", {
  ## One S3 method instead of four inherits() branches in R/fitdiag.R.
  d <- zt_dat(dist = "ztnb", n = 400, seed = 13)
  f <- ci(".countimp_zt_fit")(d$X, d$y, "ztnb")
  class(f) <- "countimp_zt_fit"
  ll <- stats::logLik(f)
  expect_equal(as.numeric(ll), f$ll)
  expect_equal(attr(ll, "df"), 4L)          # 3 beta + theta
  expect_equal(stats::AIC(f), -2 * f$ll + 2 * 4)
  expect_equal(ci(".countimp_npar")(f), 4)
  expect_true(ci(".countimp_fit_usable")(f))
})
