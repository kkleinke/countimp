## B82: Conway-Maxwell-Poisson for underdispersed counts.
##
## The internal kernel is checked against an INDEPENDENT full grid from 0, not
## against itself: the windowed grid is the thing under test, so a reference
## that shares its window would prove nothing. Where glmmTMB is installed the
## fit is additionally compared against it, but the suite does not require it --
## countimp must not depend on a C++ toolchain for this family.

test_that("moments on the window match a full grid from zero", {
  ## Independent reference: plain sum over 0:top in R arithmetic.
  ## top is per case, not fixed: at nu = 0.2, lam = 8 the mean is 32770, so a
  ## grid of 2e4 truncates the distribution and the reference itself is wrong
  ## (measured: mu = 19990 against the true 32770). expect_lt(rand) below is
  ## what caught that -- a reference has to prove its own adequacy.
  voll <- function(lam, nu, top) {
    j <- 0:top
    lp <- j * log(lam) - nu * lgamma(j + 1)
    mx <- max(lp); p <- exp(lp - mx); s <- sum(p); p <- p / s
    m1 <- sum(p * j); lf <- lgamma(j + 1); e_lf <- sum(p * lf)
    list(logZ = mx + log(s), mu = m1, va = sum(p * j * j) - m1^2,
         elf = e_lf, cov = sum(p * j * lf) - m1 * e_lf,
         rand = p[top + 1L] / max(p))
  }
  for (nu in c(0.2, 0.5, 1, 2, 5)) {
    for (lam in c(0.5, 2, 8)) {
      ## Grid width from the asymptotic mean lam^(1/nu) plus a wide margin.
      top <- max(4e3, ceiling(3 * exp(log(lam) / nu) + 40 * sqrt(exp(log(lam) / nu) / nu)))
      r <- voll(lam, nu, top)
      expect_lt(r$rand, 1e-30)   # reference grid really is wide enough
      m <- ci(".countimp_cmp_mom")(lam, nu)
      expect_false(m$capped)
      expect_equal(m$logZ, r$logZ, tolerance = 1e-10)
      expect_equal(m$mu,   r$mu,   tolerance = 1e-10)
      expect_equal(m$va,   r$va,   tolerance = 1e-9)
      expect_equal(m$elf,  r$elf,  tolerance = 1e-10)
      expect_equal(m$cov,  r$cov,  tolerance = 1e-9)
    }
  }
})

test_that("nu = 1 reproduces the Poisson exactly", {
  for (mu in c(0.5, 3, 40)) {
    m <- ci(".countimp_cmp_mom")(mu, 1)
    expect_equal(m$mu, mu, tolerance = 1e-10)
    expect_equal(m$va, mu, tolerance = 1e-10)
    expect_equal(m$logZ, mu, tolerance = 1e-10)  # Z = exp(lam)
  }
})

test_that("nu > 1 underdisperses and nu < 1 overdisperses", {
  vm <- vapply(c(0.25, 0.5, 1, 2, 4), function(nu) {
    s <- ci(".countimp_cmp_lam")(8, nu)
    s$mom$va / s$mom$mu
  }, numeric(1))
  expect_true(all(diff(vm) < 0))          # monotone in nu
  expect_gt(vm[1], 1)                      # nu = 0.25 overdispersed
  expect_equal(vm[3], 1, tolerance = 1e-6) # nu = 1 is Poisson
  expect_lt(vm[5], 1)                      # nu = 4 underdispersed
})

test_that("the mean solve hits the target at small nu", {
  ## Regression test. The Newton step was limited to +-5 in log lam, but the
  ## mode is lam^(1/nu), so that is a step of exp(5/nu) in the mode. Measured
  ## before the fix: solving for mu = 6 at nu = 0.05 returned lam = 59, whose
  ## true mean is 3e5 -- silently, with no warning and no capped flag.
  for (nu in c(0.02, 0.05, 0.1, 0.5, 1, 3)) {
    for (mu in c(1, 6, 30)) {
      s <- ci(".countimp_cmp_lam")(mu, nu)
      expect_false(is.na(s$mom$mu[1]))
      expect_equal(s$mom$mu[1], mu, tolerance = 1e-6)
    }
  }
})

test_that("the analytic gradient matches finite differences", {
  set.seed(4)
  n <- 400; x <- rnorm(n); X <- cbind(1, x)
  for (nu in c(0.7, 1, 2.2)) {
    y <- ci(".countimp_rcmp")(exp(drop(X %*% c(1.2, 0.35))), nu)
    par <- c(1.2, 0.35, log(nu))
    g <- ci(".countimp_cmp_nll")(par, X, y, gr = TRUE)
    h <- 1e-5
    fd <- vapply(seq_along(par), function(k) {
      p1 <- par; p1[k] <- p1[k] + h
      p0 <- par; p0[k] <- p0[k] - h
      (ci(".countimp_cmp_nll")(p1, X, y) -
       ci(".countimp_cmp_nll")(p0, X, y)) / (2 * h)
    }, numeric(1))
    ## unname(): the gradient carries the design matrix column names, like
    ## every other fit in the package. The values are what is under test.
    expect_equal(unname(g), unname(fd), tolerance = 1e-5)
  }
})

test_that("the gradient in beta equals the Poisson score at nu = 1", {
  ## Independent analytic check: at nu = 1 the COM-Poisson IS the Poisson, so
  ## d l / d beta must be sum (y - mu) x. This does not rely on the grid.
  set.seed(9)
  n <- 300; x <- rnorm(n); X <- cbind(1, x)
  y <- rpois(n, exp(drop(X %*% c(1, 0.4))))
  par <- c(1, 0.4, 0)
  g <- ci(".countimp_cmp_nll")(par, X, y, gr = TRUE)
  mu <- exp(drop(X %*% par[1:2]))
  expect_equal(unname(g[1:2]), unname(-drop(crossprod(X, y - mu))),
               tolerance = 1e-8)
})

test_that("draws follow the target distribution", {
  set.seed(21)
  for (nu in c(0.6, 1, 3)) {
    lam <- ci(".countimp_cmp_lam")(5, nu)$lam
    z <- ci(".countimp_rcmp")(rep(5, 20000), nu)
    j <- 0:40
    lp <- j * log(lam) - nu * lgamma(j + 1)
    p <- exp(lp - max(lp)); p <- p / sum(p)
    obs <- tabulate(z + 1L, nbins = 41L) / length(z)
    expect_lt(max(abs(obs - p)), 0.012)
    expect_equal(mean(z), 5, tolerance = 0.1)
  }
})

test_that("the fit recovers known parameters", {
  set.seed(17)
  n <- 2000; x <- rnorm(n); X <- cbind(1, x); b <- c(1.4, 0.35)
  for (nu in c(0.8, 1, 2.5)) {
    y <- ci(".countimp_rcmp")(exp(drop(X %*% b)), nu)
    f <- ci(".countimp_cmp_fit")(X, y)
    expect_equal(f$conv, 0L)
    expect_equal(unname(f$beta), b, tolerance = 0.08)
    ## unname(): .countimp_zt_fit() and .countimp_bd_fit() both return theta
    ## with an empty name, so this follows the house convention.
    expect_equal(unname(f$theta), nu, tolerance = 0.15 * nu)
    ## Field contract shared with .countimp_1l_fit()/.countimp_zt_fit().
    expect_true(all(c("beta", "cov", "scale", "theta", "ll", "conv",
                      "nobs", "npar") %in% names(f)))
    expect_null(attr(f, "class"))
  }
})

test_that("non-representable parameters are reported, not swallowed", {
  ## Each of these silently returned nonsense at some point during development.
  expect_true(ci(".countimp_cmp_mom")(Inf, 2)$capped)
  expect_true(ci(".countimp_cmp_mom")(10, Inf)$capped)
  expect_true(ci(".countimp_cmp_mom")(10, 0)$capped)
  expect_true(ci(".countimp_cmp_mom")(-1, 2)$capped)
  ## A window wider than the cell budget must be refused BEFORE allocation:
  ## clamping and building the matrix anyway hung the fit for minutes.
  expect_true(ci(".countimp_cmp_mom")(1e5, 1e-3)$capped)
  ## The likelihood turns that into a finite penalty so optim() walks back.
  expect_equal(ci(".countimp_cmp_nll")(c(50, 0, 0), cbind(1, 0), 1), 1e300)
  expect_equal(ci(".countimp_cmp_nll")(c(1, 0, 500), cbind(1, 0), 1), 1e300)
})

test_that("cmp imputes underdispersed counts and keeps the dispersion", {
  set.seed(31)
  n <- 400; x <- rnorm(n)
  y <- rbinom(n, 14, plogis(-0.2 + 0.45 * x))   # conditionally underdispersed
  d <- data.frame(y = y, x = x)
  d$y[runif(n) < 0.3] <- NA
  for (meth in c("cmp", "cmp.boot")) {
    imp <- countimp(d, method = c(y = meth, x = ""), m = 3, printFlag = FALSE)
    z <- unlist(imp$imp$y)
    expect_true(all(is.finite(z)))
    expect_true(all(z >= 0))
    expect_true(all(z == round(z)))
    ## Poisson imputations inflate the dispersion of the completed data; the
    ## point of this family is that these do not.
    vm <- vapply(1:3, function(k) {
      v <- countimp_complete(imp, k)$y; stats::var(v) / mean(v)
    }, numeric(1))
    expect_lt(mean(vm), 1.05)
  }
})

test_that("the diagnostic reports conditional, not marginal, dispersion", {
  ## Regression test. The marginal variance/mean ratio also contains the spread
  ## the covariates induce in the mean, so with a strong predictor it reads as
  ## OVERdispersion on data that are underdispersed given x -- the opposite
  ## verdict from the one the AIC ranking gives.
  set.seed(11)
  n <- 600; x <- rnorm(n)
  y <- rbinom(n, 14, plogis(-0.2 + 0.8 * x))
  d <- data.frame(y = y, x = x)
  fd <- countimp_fit_diag(y ~ x, data = d)
  expect_gt(fd$var_mean_ratio, 1.2)        # marginal says overdispersed
  expect_lt(fd$disp_conditional, 0.8)      # conditional says under
  expect_match(paste(capture.output(print(fd)), collapse = " "),
               "UNDERdispersed")
  expect_equal(fd$recommendation, "cmp")
})

test_that("cmp ranks in the diagnostic without glmmTMB", {
  ## The candidate used to be fitted by glmmTMB, so the one model that detects
  ## underdispersion disappeared on machines without a C++ toolchain.
  set.seed(13)
  n <- 400; x <- rnorm(n)
  y <- rbinom(n, 12, plogis(0.1 + 0.4 * x))
  fd <- countimp_fit_diag(y ~ x, data = data.frame(y = y, x = x))
  tb <- fd$table
  row <- tb[tb$method == "cmp", ]
  expect_equal(nrow(row), 1L)
  expect_equal(row$status, "ok")
  expect_true(is.finite(row$aic))
  ## P(Y = 0) is exact from the normalising constant, not NA as it was under
  ## the glmmTMB route.
  expect_true(is.finite(row$p0_model))
  expect_gte(row$p0_model, 0)
  expect_lte(row$p0_model, 1)
})

test_that("the recommended method string is executable", {
  ## An unusable recommendation is a bug, not a hint (same rule as B61).
  expect_true(exists("mice.impute.cmp", where = asNamespace("countimp")))
  expect_true(exists("mice.impute.cmp.boot", where = asNamespace("countimp")))
})
