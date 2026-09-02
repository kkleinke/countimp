## B61: countimp_fit_diag() -- which count family fits the observed data?
##
## Four bugs found while building this, each recorded because the failure mode
## is not obvious from reading the code:
##
##   1. glmmTMB does not export poisson(); glmmTMB::poisson() errors. The ZIP
##      candidate was silently dropped from every ranking.
##   2. glmmTMB reports fit$convergence == 0 even when the Hessian is not
##      positive definite, so a model that has not identified its parameters
##      looks converged. Finiteness of logLik() catches it; that flag does not.
##   3. Fitting the zero-inflated candidates with ziformula = ~1 while the
##      hurdle candidates get all predictors in their zero model handicaps the
##      ZI family by however much signal sits in those predictors. Measured:
##      ZIP data ranked hurdle Poisson first at dAIC 22 in favour of the wrong
##      family; with a comparable zero part the ZI model wins.
##   4. dpois(0, mean(y)) is not the zero share a Poisson model predicts when
##      the covariates spread the rate: 4.7% marginal against 9.6% from the
##      fitted model on the same data, which reads as excess zeros where there
##      are none.

## Draw from a known process. Kept local rather than in a helper so the
## generating truth is next to the assertion about it.
b61_sim <- function(kind, n = 800, seed = 1) {
  set.seed(seed)
  x1 <- rnorm(n); x2 <- rnorm(n)
  mu <- exp(1 + 0.4 * x1 + 0.3 * x2)
  y <- switch(kind,
    poisson = rpois(n, mu),
    nb      = MASS::rnegbin(n, mu = mu, theta = 1.5),
    zip     = ifelse(runif(n) < plogis(-0.4 + 1.5 * x2), 0, rpois(n, mu)),
    zinb    = ifelse(runif(n) < plogis(-0.4 + 1.5 * x2), 0,
                     MASS::rnegbin(n, mu = mu, theta = 1.5)),
    hp      = { z <- rbinom(n, 1, plogis(0.3 + 1.2 * x2)); v <- integer(n)
                v[z == 1] <- ci(".countimp_rktp")(sum(z == 1), 0, mu[z == 1]); v },
    hnb     = { z <- rbinom(n, 1, plogis(0.3 + 1.2 * x2)); v <- integer(n)
                v[z == 1] <- ci(".countimp_rktnb")(sum(z == 1), size = 1.5, k = 0,
                                             mu = mu[z == 1]); v },
    stop("unknown kind"))
  data.frame(y = y, x1 = x1, x2 = x2)
}

test_that("the generating family is recommended for six known processes", {
  ## Six full countimp_fit_diag() runs at n = 800, each fitting every
  ## candidate family: 76 s. The size is what the claim needs -- at a
  ## smaller n the families stop being distinguishable and the test would
  ## assert less than it says. Kept out of CRAN's budget, not weakened.
  skip_on_cran()
  ## pscl is in Suggests, not in Imports: the one-level zero-inflation and
  ## hurdle methods need it, the rest of the package does not. Without this
  ## skip the test fails in an environment without Suggests -- and it fails on
  ## the missing dependency, not on the property under test.
  skip_if_not_installed("pscl")
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("MASS")
  ## One seed per process. The full 25-run measurement lives in the docs (141
  ## of 150); here we pin one run each so a regression shows up immediately.
  for (k in c("poisson", "nb", "zip", "zinb", "hp", "hnb")) {
    r <- countimp_fit_diag(y ~ x1 + x2, b61_sim(k, seed = 501))
    expect_equal(r$recommendation, k, info = paste("process:", k))
  }
})

test_that("every candidate is either ranked or given a reason", {
  skip_on_cran()          # 64 s: fits every candidate family at n = 800
  ## pscl is in Suggests, not in Imports: the one-level zero-inflation and
  ## hurdle methods need it, the rest of the package does not. Without this
  ## skip the test fails in an environment without Suggests -- and it fails on
  ## the missing dependency, not on the property under test.
  skip_if_not_installed("pscl")
  skip_if_not_installed("glmmTMB")
  r <- countimp_fit_diag(y ~ x1 + x2, b61_sim("zip", seed = 11))
  ## Bug 1: a candidate must never vanish without explanation.
  expect_true(all(r$table$status != ""))
  expect_true(all(is.finite(r$table$aic) | !r$table$fitted |
                    r$table$family == "quasipoisson"))
  ## Quasi-Poisson converges but has no likelihood: fitted, no AIC.
  qp <- r$table[r$table$family == "quasipoisson", ]
  expect_true(qp$fitted)
  expect_true(is.na(qp$aic))
})

test_that("the ZIP candidate is fitted, not silently dropped", {
  skip_on_cran()          # 50 s: full candidate set at n = 800
  ## pscl is in Suggests, not in Imports: the one-level zero-inflation and
  ## hurdle methods need it, the rest of the package does not. Without this
  ## skip the test fails in an environment without Suggests -- and it fails on
  ## the missing dependency, not on the property under test.
  skip_if_not_installed("pscl")
  ## Bug 1 directly: glmmTMB::poisson() does not exist. If the family call
  ## regresses, this candidate fails and the test catches it.
  skip_if_not_installed("glmmTMB")
  r <- countimp_fit_diag(y ~ x1 + x2, b61_sim("zip", seed = 11))
  expect_true(r$table$fitted[r$table$family == "zip"])
  expect_true(is.finite(r$table$aic[r$table$family == "zip"]))
})

test_that("zero-inflated candidates get the same zero-part predictors as hurdle", {
  skip_on_cran()          # 61 s: two full candidate sets at n = 800
  ## pscl is in Suggests, not in Imports: the one-level zero-inflation and
  ## hurdle methods need it, the rest of the package does not. Without this
  ## skip the test fails in an environment without Suggests -- and it fails on
  ## the missing dependency, not on the property under test.
  skip_if_not_installed("pscl")
  ## Bug 3. Comparing families on AIC requires each candidate be given its best
  ## shot. With ziformula = ~1 the ZI candidate loses to hurdle on ZIP data.
  skip_if_not_installed("glmmTMB")
  d <- b61_sim("zip", seed = 11)
  r <- countimp_fit_diag(y ~ x1 + x2, d)
  aic.zip <- r$table$aic[r$table$family == "zip"]
  aic.hp  <- r$table$aic[r$table$family == "hp"]
  expect_true(is.finite(aic.zip))
  expect_lt(aic.zip, aic.hp)
})

test_that("a separated zero-inflation fit is rejected, not ranked", {
  skip_on_cran()          # 10.5 s here, ~85 s on Windows, for one assertion
  ## pscl is in Suggests, not in Imports: the one-level zero-inflation and
  ## hurdle methods need it, the rest of the package does not. Without this
  ## skip the test fails in an environment without Suggests -- and it fails on
  ## the missing dependency, not on the property under test.
  skip_if_not_installed("pscl")
  ## Bug 2 plus the degeneracy rule. On genuine Poisson data the ZI candidate
  ## can separate (fitted zprob spanning the unit interval) and still report a
  ## good likelihood. Such a fit must not win.
  skip_if_not_installed("glmmTMB")
  bad <- 0L
  for (i in 1:8) {
    r <- countimp_fit_diag(y ~ x1 + x2, b61_sim("poisson", seed = 100 + i))
    if (r$recommendation != "poisson") bad <- bad + 1L
    win <- r$table$family[1]
    if (win %in% c("zip", "zinb") && r$table$fitted[1]) {
      zp <- stats::predict(r$fits[[win]], type = "zprob")
      expect_true(min(zp) >= 1e-6 && max(zp) <= 1 - 1e-6)
    }
  }
  ## Measured: 1 of 40 runs misses. Allowing 2 of 8 keeps the test from being
  ## flaky while still failing if the rules are removed (then it was 11 of 40).
  expect_lte(bad, 2L)
})

test_that("the parsimony rule prefers the simpler model inside 2 AIC units", {
  skip_if_not_installed("glmmTMB")
  ## Seed 109 is chosen because it produces the situation the rule exists for:
  ## on genuine Poisson data the negative binomial takes the lowest AIC by 0.41
  ## units -- inside the noise, since AIC charges 2 per parameter -- and the
  ## Poisson must still be recommended. Written unconditionally: a version of
  ## this test that only asserted when the band happened to be wide enough
  ## passed against a build with the rule removed.
  r <- countimp_fit_diag(y ~ x1 + x2, b61_sim("poisson", seed = 109))
  expect_equal(r$recommendation, "poisson")
  expect_true(isTRUE(r$parsimony_applied))
  ## The lowest-AIC candidate really is the richer one, or the seed has drifted
  ## and this test no longer exercises the rule.
  expect_false(r$table$family[1] == "poisson")
  npar.rec <- ci(".countimp_npar")(r$fits$poisson)
  npar.top <- ci(".countimp_npar")(r$fits[[r$table$family[1]]])
  expect_lt(npar.rec, npar.top)
  ## And the print method must disclose that the top-ranked model was passed over.
  out <- paste(capture.output(print(r)), collapse = " ")
  expect_match(out, "not the lowest-AIC candidate")
})

test_that("the Poisson-implied zero share comes from the fitted model", {
  ## Bug 4: the marginal dpois(0, mean(y)) understates it whenever the
  ## covariates spread the rate.
  d <- b61_sim("poisson", seed = 11)
  r <- countimp_fit_diag(y ~ x1 + x2, d)
  marginal <- mean(dpois(0, mean(d$y)))
  expect_gt(r$p0_poisson, marginal)
  ## And it should be close to what was observed, since the data ARE Poisson.
  expect_lt(abs(r$p0_poisson - r$p0_observed), 0.05)
})

test_that("underdispersion is reported even though no method handles it yet", {
  ## The moment evidence must carry the signal when no candidate family can.
  set.seed(4); n <- 600; x <- rnorm(n)
  d <- data.frame(y = rbinom(n, 14, plogis(-0.3 + 0.35 * x)), x = x)
  r <- countimp_fit_diag(y ~ x, d)
  expect_lt(r$var_mean_ratio, 0.9)
  out <- paste(capture.output(print(r)), collapse = " ")
  expect_match(out, "UNDERdispersed")
})

test_that("input is validated", {
  d <- data.frame(y = c(1.5, 2.5, 3.5), x = 1:3)
  expect_error(countimp_fit_diag(y ~ x, d), "integer counts")
  d2 <- data.frame(y = c(-1L, 2L, 3L), x = 1:3)
  expect_error(countimp_fit_diag(y ~ x, d2), "integer counts")
  d3 <- data.frame(y = rpois(10, 2), x = rnorm(10))
  expect_error(countimp_fit_diag(y ~ x, d3), "at least 20")
  expect_error(countimp_fit_diag("y ~ x", d3), "must be a formula")
  expect_error(countimp_fit_diag(y ~ x, list(y = 1)), "data frame")
  d4 <- data.frame(y = rpois(50, 2), x = rnorm(50))
  expect_error(countimp_fit_diag(y ~ x, d4, families = "gaussian"), "unknown famil")
})

test_that("every recommendation is a method countimp() accepts", {
  ## A recommendation the user cannot paste into method= is a bug, not a hint.
  skip_if_not_installed("glmmTMB")
  cand <- ci(".countimp_candidates")()
  for (nm in names(cand)) {
    fn <- paste0("mice.impute.", cand[[nm]]$method)
    expect_true(exists(fn, where = asNamespace("countimp")),
                info = paste("no imputation function for candidate", nm))
  }
})

test_that("plot() produces a rootogram panel per candidate", {
  skip_on_cran()          # 65 s: fits every candidate, then draws it
  ## pscl is in Suggests, not in Imports: the one-level zero-inflation and
  ## hurdle methods need it, the rest of the package does not. Without this
  ## skip the test fails in an environment without Suggests -- and it fails on
  ## the missing dependency, not on the property under test.
  skip_if_not_installed("pscl")
  skip_if_not_installed("glmmTMB")
  r <- countimp_fit_diag(y ~ x1 + x2, b61_sim("zip", seed = 11))
  f <- tempfile(fileext = ".png")
  grDevices::png(f, width = 800, height = 600)
  ## The device must be closed before the file is complete on disk, so the
  ## existence check goes after dev.off(), not before it.
  ok <- tryCatch({ plot(r); TRUE }, error = function(e) conditionMessage(e))
  grDevices::dev.off()
  on.exit(unlink(f), add = TRUE)
  expect_true(isTRUE(ok), info = if (!isTRUE(ok)) ok else "")
  expect_true(file.exists(f) && file.size(f) > 1000)
})

test_that("the hurdle rootogram does not recycle mismatched predictions", {
  ## pscl is in Suggests, not in Imports: the one-level zero-inflation and
  ## hurdle methods need it, the rest of the package does not. Without this
  ## skip the test fails in an environment without Suggests -- and it fails on
  ## the missing dependency, not on the property under test.
  skip_if_not_installed("pscl")
  ## The zero component is fitted on all n cases, the count component on the
  ## positives only. Multiplying the two vectors directly recycles 403 rates
  ## against 800 zero probabilities -- no error, just a wrong rootogram and a
  ## recycling warning per count level.
  skip_if_not_installed("glmmTMB")
  r <- countimp_fit_diag(y ~ x1 + x2, b61_sim("hp", seed = 11))
  h <- r$fits$hp
  skip_if(is.null(h))
  expect_equal(length(stats::fitted(h$zero)), nrow(h$frame))
  e <- ci(".countimp_expected_counts")(h, "hp", 10)
  expect_false(is.null(e))
  expect_true(all(is.finite(e)))
  ## Expected counts over all levels cannot exceed the sample size.
  expect_lte(sum(e), nrow(h$frame) * 1.001)
  n_w <- 0L
  withCallingHandlers(ci(".countimp_expected_counts")(h, "hp", 10),
                      warning = function(w) { n_w <<- n_w + 1L
                                              invokeRestart("muffleWarning") })
  expect_equal(n_w, 0L)
})
