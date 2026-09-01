## B75-B78: interval-truncated counts (bounds)
##
## What these tests defend, in order of how easily each was got wrong:
##
##  B75  the analytic gradient, component by component against finite
##       differences. The estimates look right even with a sign error in the
##       truncation term (B69), so a coefficient comparison would not catch it.
##  B76  countimp() accepting the arguments of its own methods. The
##       hand-written allowlist rejected EV, bounds, donors and theta.
##  B77  the internal-fit contract: field names beta/cov/scale/theta/ll/conv/
##       nobs/npar and NO class attribute.
##  B78  the extreme-mu regime, where the naive log pmf loses every
##       k-dependent digit and the weights come out uniform.

## `ci()` is the helper for reaching package internals (helper-internals.R).
## The imputation shortcut therefore needs a different name.
imp_q <- function(...) suppressWarnings(countimp(..., printFlag = FALSE))

## --- B75: gradient ----------------------------------------------------------

test_that("B75: the interval-truncated gradient matches finite differences", {
  nll <- ci(".countimp_bd_nll")
  gr  <- ci(".countimp_bd_gr")
  fd <- function(par, ...) vapply(seq_along(par), function(j) {
    h <- max(1e-7, abs(par[j]) * 1e-7)
    p1 <- p2 <- par; p1[j] <- p1[j] + h; p2[j] <- p2[j] - h
    (nll(p1, ...) - nll(p2, ...)) / (2 * h)
  }, numeric(1))

  set.seed(21); n <- 500
  X <- cbind(1, rnorm(n), runif(n))
  mu <- exp(drop(X %*% c(0.9, 0.4, -0.3)))

  y <- ci(".countimp_rint")(mu, 1, 8, "poisson")
  par <- c(0.8, 0.35, -0.2)
  ga <- gr(par, X, y, 1, 8, "poisson")
  gn <- fd(par, x = X, y = y, lo = 1, hi = 8, dist = "poisson")
  expect_equal(length(ga), 3L)
  ## component by component: a sign error in one term must not be averaged away
  for (j in seq_along(ga))
    expect_lt(abs(ga[j] - gn[j]) / max(abs(gn[j]), 1e-8), 1e-5)

  yn <- ci(".countimp_rint")(mu, 2, 12, "negbin", theta = 1.8)
  parn <- c(0.85, 0.38, -0.25, log(1.6))
  ga <- gr(parn, X, yn, 2, 12, "negbin")
  gn <- fd(parn, x = X, y = yn, lo = 2, hi = 12, dist = "negbin")
  expect_equal(length(ga), 4L)
  for (j in seq_along(ga))
    expect_lt(abs(ga[j] - gn[j]) / max(abs(gn[j]), 1e-8), 1e-4)
})

test_that("B75: the fit recovers the generating parameters", {
  set.seed(22); n <- 4000
  X <- cbind(1, rnorm(n)); b <- c(1, 0.5)
  mu <- exp(drop(X %*% b))
  y <- ci(".countimp_rint")(mu, 1, 9, "poisson")
  f <- ci(".countimp_bd_fit")(X, y, "poisson", 1, 9)
  se <- sqrt(diag(f$cov))
  expect_lt(abs(f$beta[2] - b[2]), 3 * se[2])
  expect_equal(f$conv, 0L)

  yn <- ci(".countimp_rint")(mu, 0, 9, "negbin", theta = 2)
  fn <- ci(".countimp_bd_fit")(X, yn, "negbin", 0, 9)
  expect_lt(abs(fn$beta[2] - b[2]), 0.1)
  expect_lt(abs(fn$theta - 2), 0.6)
})

test_that("B75: an unbounded fit is biased where the bound binds", {
  ## The reason bounds enter the likelihood and not only the draw. If this
  ## stops holding, the extra machinery is not earning its place.
  set.seed(23); n <- 3000
  X <- cbind(1, rnorm(n)); b <- c(0.7, 0.45)
  mu <- exp(drop(X %*% b))
  y <- ci(".countimp_rint")(mu, 0, 4, "poisson")
  bu <- stats::glm.fit(X, y, family = stats::poisson())$coefficients
  bt <- ci(".countimp_bd_fit")(X, y, "poisson", 0, 4)$beta
  expect_lt(abs(bt[2] - b[2]), abs(bu[2] - b[2]) / 2)
})

## --- B76: countimp() accepts its own methods' arguments ---------------------

test_that("B76: countimp() accepts arguments declared by its own methods", {
  ## EV = TRUE tops up outlying draws through mice.impute.midastouch, and mice
  ## is in Suggests. Without mice countimp reports that clearly -- but what is
  ## under test here is something else, so the skip belongs here rather than an
  ## expectation on the error.
  skip_if_not_installed("mice")
  set.seed(24); n <- 200
  d <- data.frame(x1 = rnorm(n))
  d$y <- rpois(n, exp(0.5 + 0.4 * d$x1))
  d$y[sample(n, 40)] <- NA
  m <- c(y = "poisson", x1 = "")

  ## EV is documented on every count method; it used to be rejected as unknown
  expect_no_error(imp_q(d, method = m, m = 1, EV = TRUE))
  expect_no_error(imp_q(d, method = c(y = "nb", x1 = ""), m = 1, quiet = TRUE))
  ## and a genuine typo must still be caught, with the near miss named
  expect_error(imp_q(d, method = m, m = 1, maxti = 3), "maxit")
  expect_error(imp_q(d, method = m, m = 1, nonsense_arg = 1), "Unknown argument")
})

test_that("B76: the allowlist is derived, not hand-written", {
  ## A new method's arguments must be accepted without editing families.R.
  known <- ci(".countimp_known_args")()
  ns <- ci_home()
  ms <- grep("^mice\\.impute\\.", ls(ns), value = TRUE)
  decl <- setdiff(unique(unlist(lapply(ms, function(m)
    names(formals(get(m, envir = ns)))))),
    c("y", "ry", "x", "wy", "...", "type"))
  expect_true(length(decl) > 0L)
  expect_true(all(decl %in% known))
  expect_true(all(c("EV", "bounds", "m", "maxit", "quiet") %in% known))
})

## --- B77: internal fit contract --------------------------------------------

test_that("B77: the bounded fit follows the internal fit contract", {
  set.seed(25); n <- 400
  X <- cbind(1, rnorm(n))
  y <- ci(".countimp_rint")(exp(drop(X %*% c(1, 0.4))), 0, 7, "poisson")
  f <- ci(".countimp_bd_fit")(X, y, "poisson", 0, 7)
  ## same field names as .countimp_1l_fit() and .countimp_zt_fit()
  expect_true(all(c("beta", "cov", "scale", "theta", "ll", "conv",
                    "nobs", "npar") %in% names(f)))
  ## no class: countimp_check() recognises an internal fit by its fields, and a
  ## class attribute sent it down the coef()/vcov() path, which returns NULL
  expect_null(attr(f, "class"))
  ## and the diagnostic must accept it
  expect_length(ci(".countimp_check_fit")(f), 0L)
})

test_that("B77: a lower bound of 1 with no ceiling delegates to ztrunc", {
  set.seed(26); n <- 300
  X <- cbind(1, rnorm(n))
  y <- ci(".countimp_rktp")(n, k = 0, mu = exp(drop(X %*% c(1, 0.4))))
  f1 <- ci(".countimp_bd_fit")(X, y, "poisson", 1, Inf)
  f2 <- ci(".countimp_zt_fit")(X, y, "ztp")
  expect_equal(f1$beta, f2$beta, tolerance = 1e-8)
  ## an infinite ceiling with any other floor has no grid and must say so
  expect_error(ci(".countimp_bd_fit")(X, y, "poisson", 2, Inf), "infinite")
})

## --- B78: extreme mu -------------------------------------------------------

test_that("B78: the log kernel keeps its k-dependence at extreme mu", {
  lg <- ci(".countimp_bd_lgrid")
  lse <- ci(".countimp_bd_lse")
  ## exact against the analytic truncated pmf
  for (cfg in list(list(4, 0L, 8L, "poisson", NULL),
                   list(5, 2L, 12L, "negbin", 1.5))) {
    k <- cfg[[2]]:cfg[[3]]
    w <- exp(lg(log(cfg[[1]]), cfg[[2]], cfg[[3]], cfg[[4]], cfg[[5]]) -
               lse(lg(log(cfg[[1]]), cfg[[2]], cfg[[3]], cfg[[4]], cfg[[5]])))
    ex <- if (cfg[[4]] == "negbin")
      stats::dnbinom(k, size = cfg[[5]], mu = cfg[[1]]) else
        stats::dpois(k, cfg[[1]])
    expect_lt(max(abs(as.numeric(w) - ex / sum(ex))), 1e-12)
  }
  ## the regime that broke: mass must sit at the ceiling, not spread uniformly
  for (mu in c(1e2, 1e10, 1e300)) {
    l <- lg(log(mu), 1L, 7L, "poisson")
    w <- as.numeric(exp(l - lse(l)))
    expect_equal(which.max(w), 7L)
    if (mu >= 1e10) expect_gt(w[7], 1 - 1e-8)
  }
  ## and at the floor for tiny mu
  l <- lg(log(1e-300), 2L, 9L, "poisson")
  w <- as.numeric(exp(l - lse(l)))
  expect_equal(which.max(w), 1L)
})

test_that("B78: draws stay inside the bounds at extreme mu", {
  set.seed(27)
  z <- ci(".countimp_rint")(rep(1e50, 200), 1, 7, "poisson")
  expect_true(all(z == 7))
  z <- ci(".countimp_rint")(rep(1e-50, 200), 3, 9, "poisson")
  expect_true(all(z == 3))
  ## a degenerate one-value interval
  expect_true(all(ci(".countimp_rint")(rep(4, 50), 2, 2, "poisson") == 2))
})

test_that("B78: clamping is not equivalent to a truncated draw", {
  ## The reason .countimp_rint() exists. Clamping errs in both directions, so
  ## the test checks that it is wrong, not the direction of the error.
  set.seed(28)
  vr <- vapply(c(1, 4, 20), function(mu) {
    tr <- ci(".countimp_rint")(rep(mu, 20000), 1, 7, "poisson")
    cl <- pmin(pmax(stats::rpois(20000, mu), 1), 7)
    stats::var(cl) / stats::var(tr)
  }, numeric(1))
  ## measured values: 0.74, 1.17, 0.001 (analyse/k21_bounds_klemmen.csv)
  expect_lt(vr[1], 0.95)   # near the floor: clipping DEFLATES the variance
  expect_gt(vr[2], 1.05)   # in the interior: it INFLATES it
  expect_lt(vr[3], 0.1)    # far above the ceiling: variance destroyed
})

## --- end to end -------------------------------------------------------------

test_that("bounded methods respect the bounds through countimp()", {
  set.seed(29); n <- 400
  d <- data.frame(x1 = rnorm(n))
  d$y <- rbinom(n, 8, plogis(0.2 + 0.6 * d$x1))
  d$y[sample(n, 90)] <- NA
  for (meth in c("bp", "bp.boot", "bnb", "bnb.boot")) {
    imp <- imp_q(d, method = c(y = meth, x1 = ""), m = 2, bounds = c(0, 8))
    z <- unlist(imp$imp$y)
    expect_true(all(z >= 0 & z <= 8), info = meth)
    expect_equal(sum(is.na(z)), 0L, info = meth)
  }
})

test_that("bounds may be given per variable", {
  set.seed(30); n <- 400
  d <- data.frame(x1 = rnorm(n))
  d$y <- rbinom(n, 8, plogis(0.2 + 0.6 * d$x1))
  d$days <- rbinom(n, 5, 0.4)
  d$y[sample(n, 80)] <- NA
  d$days[sample(n, 70)] <- NA
  m <- c(y = "bp", x1 = "", days = "bp")
  imp <- imp_q(d, method = m, m = 2, bounds = list(y = c(0, 8), days = c(0, 5)))
  expect_true(all(unlist(imp$imp$y) <= 8))
  expect_true(all(unlist(imp$imp$days) <= 5))
  ## a variable left out of the list is an error, not silent unbounded imputation
  expect_error(imp_q(d, method = m, m = 1, bounds = list(y = c(0, 8))),
               "does not cover variable 'days'")
})

test_that("impossible bounds specifications are rejected", {
  set.seed(31); n <- 200
  d <- data.frame(x1 = rnorm(n))
  d$y <- rbinom(n, 8, 0.5); d$y[sample(n, 40)] <- NA
  m <- c(y = "bp", x1 = "")
  expect_error(imp_q(d, method = m, m = 1, bounds = c(0, 3)), "outside bounds")
  expect_error(imp_q(d, method = m, m = 1, bounds = c(5, 5)), "lower < upper")
  expect_error(imp_q(d, method = m, m = 1, bounds = c(0, 8.4)), "whole numbers")
  expect_error(imp_q(d, method = m, m = 1, bounds = c(-2, 8)), "cannot be negative")
  expect_error(imp_q(d, method = m, m = 1, bounds = 5), "length 2")
  expect_error(imp_q(d, method = m, m = 1), "required")
  ## a grid wide enough to be pointless is refused rather than made slow
  expect_error(imp_q(d, method = m, m = 1, bounds = c(0, 99999)), "grid limit")
})
