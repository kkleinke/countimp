## B83: right-censored counts (cp, cnb, cp.boot, cnb.boot)
##
## What these tests defend, in the order in which each was easy to get wrong:
##
##  1. the analytic gradient, component by component. The NB censored term
##     carries a factor (theta + c - 1)/(theta + mu) that reduces to 1 as
##     theta -> Inf; putting it the other way up leaves the fitted coefficients
##     almost unchanged on mild censoring, so only a finite-difference check
##     finds it. Mutation probe below.
##  2. the two grid-free identities: censor = Inf must reproduce glm()/glm.nb()
##     exactly, and the censored likelihood must equal an independent grid sum
##     of the tail. Both are references OUTSIDE this file's own code.
##  3. the scale rule. Which scale the imputations live on is decided from
##     `where`, and getting it wrong is silent: the column simply contains
##     values it should not.
##  4. a resolved censored cell must be >= its limit. That is the whole content
##     of "the observed value is a lower bound"; drawing it from the untruncated
##     distribution would throw the information away and look perfectly normal.
##  5. the data-error checks (values above the limit, missing censor, all cases
##     censored), which must fire BEFORE the retry loop for the reason recorded
##     as B73: inside it, a data problem is reported as a convergence failure.

imp_c <- function(...) suppressWarnings(countimp(..., printFlag = FALSE))

## One censored data set, so every block measures the same thing.
.b83_daten <- function(n = 400L, cens = 6, seed = 83L, nmis = 100L) {
  set.seed(seed)
  x1 <- stats::rnorm(n)
  ystar <- stats::rpois(n, exp(1.4 + 0.5 * x1))
  d <- data.frame(y = pmin(ystar, cens), x1 = x1)
  d$ystar <- ystar
  d$y[sample.int(n, nmis)] <- NA
  d
}

## --- 1: gradient ------------------------------------------------------------

test_that("B83: the censored gradient matches finite differences", {
  nll <- ci(".countimp_cs_nll")
  gr  <- ci(".countimp_cs_gr")
  fd <- function(par, ...) vapply(seq_along(par), function(j) {
    h <- max(1e-6, abs(par[j]) * 1e-6)
    p1 <- p2 <- par; p1[j] <- p1[j] + h; p2[j] <- p2[j] - h
    (nll(p1, ...) - nll(p2, ...)) / (2 * h)
  }, numeric(1))

  set.seed(83); n <- 500
  X <- cbind(1, rnorm(n), runif(n))
  mu <- exp(drop(X %*% c(1.4, 0.5, -0.3)))
  cc <- rep(6, n)

  y <- pmin(rpois(n, mu), 6)
  expect_gt(mean(y >= 6), 0.15)          # the censored term must actually bind
  par <- c(1.3, 0.45, -0.25)
  ga <- gr(par, X, y, cc, "poisson")
  gn <- fd(par, x = X, y = y, cc = cc, dist = "poisson")
  for (j in seq_along(ga))
    expect_equal(ga[j], gn[j], tolerance = 1e-5)

  ## NB, including the log-theta component, which uses the finite sum over
  ## 0..c-1 rather than a closed form
  y2 <- pmin(MASS::rnegbin(n, mu = mu, theta = 1.8), 6)
  par2 <- c(1.3, 0.45, -0.25, log(1.6))
  ga2 <- gr(par2, X, y2, cc, "negbin")
  gn2 <- fd(par2, x = X, y = y2, cc = cc, dist = "negbin")
  expect_equal(length(ga2), 4L)
  for (j in seq_along(ga2))
    expect_equal(ga2[j], gn2[j], tolerance = 1e-5)

  ## Mutation probe: the NB censored term with the hazard factor inverted. The
  ## check must fail on it -- otherwise it proves nothing about the factor.
  gr_mut <- function(par, x, y, cc, dist) {
    p <- ncol(x); th <- exp(par[p + 1L])
    eta <- pmin(pmax(drop(x %*% par[seq_len(p)]), -700), 700); mu <- exp(eta)
    cen <- y >= cc
    d <- numeric(length(y))
    d[!cen] <- (y[!cen] - mu[!cen]) * th / (th + mu[!cen])
    k <- cc[cen] - 1; m <- mu[cen]
    lS <- stats::pnbinom(k, size = th, mu = m, lower.tail = FALSE, log.p = TRUE)
    lf <- stats::dnbinom(k, size = th, mu = m, log = TRUE)
    d[cen] <- exp(lf - lS + log(m) + log(th + m) - log(th + k))   # inverted
    -drop(crossprod(x, d))
  }
  gm <- gr_mut(par2, X, y2, cc, "negbin")
  expect_false(isTRUE(all.equal(gm[2], gn2[2], tolerance = 1e-5)))
})


## --- 2: identities against references outside this code ---------------------

test_that("B83: censor = Inf reproduces the uncensored ML fit exactly", {
  set.seed(831); n <- 400
  x1 <- rnorm(n)
  y <- rpois(n, exp(1.2 + 0.4 * x1))
  X <- cbind(1, x1)

  f <- ci(".countimp_cs_fit")(X, y, rep(Inf, n), "poisson")
  g <- stats::glm(y ~ x1, family = stats::poisson())
  expect_equal(as.numeric(f$beta), as.numeric(coef(g)), tolerance = 1e-6)
  expect_equal(f$ll, as.numeric(logLik(g)), tolerance = 1e-8)
  ## the covariance too: the Hessian route must agree with the IRLS one
  expect_equal(sqrt(diag(f$cov)), as.numeric(summary(g)$coefficients[, 2]),
               tolerance = 1e-4, ignore_attr = TRUE)

  y2 <- MASS::rnegbin(n, mu = exp(1.2 + 0.4 * x1), theta = 2)
  f2 <- ci(".countimp_cs_fit")(X, y2, rep(Inf, n), "negbin")
  g2 <- MASS::glm.nb(y2 ~ x1)
  expect_equal(as.numeric(f2$beta), as.numeric(coef(g2)), tolerance = 1e-4)
  expect_equal(as.numeric(f2$theta), as.numeric(g2$theta), tolerance = 1e-3)
})

test_that("B83: the censored log-likelihood equals an independent tail sum", {
  ## Reference: P(Y >= c) summed term by term from c upwards, with a check that
  ## the reference itself runs far enough (the lesson of test-B82: a reference
  ## that stops too early can fail a correct implementation).
  nll <- ci(".countimp_cs_nll")
  set.seed(832); n <- 60
  X <- cbind(1, rnorm(n))
  par <- c(1.1, 0.3)
  mu <- exp(drop(X %*% par))
  cc <- rep(5, n)
  y <- pmin(rpois(n, mu), 5)
  cen <- y >= cc

  obergrenze <- max(cc) + 400L
  rest <- stats::ppois(obergrenze, max(mu), lower.tail = FALSE)
  expect_lt(rest, 1e-30)                       # the reference is wide enough
  ref <- sum(stats::dpois(y[!cen], mu[!cen], log = TRUE)) +
    sum(vapply(which(cen), function(i)
      log(sum(stats::dpois(cc[i]:obergrenze, mu[i]))), numeric(1)))
  expect_equal(nll(par, X, y, cc, "poisson"), -ref, tolerance = 1e-10)
})


## --- 3: the fit is what removes the bias ------------------------------------

test_that("B83: the censored fit recovers the slope the naive fit loses", {
  set.seed(833); n <- 3000
  x1 <- rnorm(n)
  ystar <- rpois(n, exp(1.4 + 0.55 * x1))
  y <- pmin(ystar, 5)
  X <- cbind(1, x1)
  expect_gt(mean(y >= 5), 0.2)

  b_cens <- ci(".countimp_cs_fit")(X, y, rep(5, n), "poisson")$beta[2]
  b_naiv <- stats::glm.fit(X, y, family = stats::poisson())$coefficients[2]
  ## the censored fit is close to the truth, the naive one is not
  expect_lt(abs(b_cens - 0.55), 0.05)
  expect_gt(abs(b_naiv - 0.55), 0.15)
})


## --- 4: the draw, and the scale rule ----------------------------------------

test_that("B83: without censored cells in `where` imputations stay at the limit", {
  d <- .b83_daten()[, c("y", "x1")]
  imp <- imp_c(d, method = c(y = "cp", x1 = ""), m = 2, censor = 6)
  vals <- unlist(imp$imp$y)
  expect_true(all(vals <= 6))
  expect_true(all(vals >= 0))
  expect_true(any(vals == 6))               # the cap is reached, not avoided
})

test_that("B83: censored cells marked in `where` are resolved above their limit", {
  d <- .b83_daten()[, c("y", "x1")]
  wm <- matrix(FALSE, nrow(d), ncol(d), dimnames = dimnames(d))
  wm[, "y"] <- is.na(d$y) | (!is.na(d$y) & d$y >= 6)
  imp <- imp_c(d, method = c(y = "cp", x1 = ""), m = 2, where = wm, censor = 6)

  ## every resolved censored observation must be at least its limit -- this is
  ## the left truncation, and it is the point of the whole method
  zeilen <- rownames(imp$imp$y)
  cen_zeilen <- zeilen %in% rownames(d)[!is.na(d$y) & d$y >= 6]
  expect_gt(sum(cen_zeilen), 10L)
  expect_true(all(imp$imp$y[cen_zeilen, ] >= 6))
  ## and the missing cells are now on the same (latent) scale, so they may pass
  ## the limit; if they could not, the column would carry two scales
  expect_true(any(imp$imp$y[!cen_zeilen, ] > 6))
})

test_that("B83: resolving the censoring restores the latent mean", {
  ## The substantive claim of the method: completing a top-coded column brings
  ## its mean back towards the latent one instead of leaving it understated.
  d <- .b83_daten(n = 800L, nmis = 0L)
  dd <- d[, c("y", "x1")]
  wm <- matrix(FALSE, nrow(dd), ncol(dd), dimnames = dimnames(dd))
  wm[, "y"] <- dd$y >= 6
  imp <- imp_c(dd, method = c(y = "cp", x1 = ""), m = 5,
               where = wm, censor = 6)
  mittel <- mean(vapply(1:5, function(i)
    mean(countimp_complete(imp, i)$y), numeric(1)))
  expect_lt(abs(mittel - mean(d$ystar)), 0.35)
  expect_gt(mittel - mean(d$y), 0.35)        # and it is a real correction
})

test_that("B83: the bootstrap and NB variants obey the same scale rule", {
  d <- .b83_daten(n = 300L, nmis = 70L)[, c("y", "x1")]
  for (meth in c("cp.boot", "cnb", "cnb.boot")) {
    imp <- imp_c(d, method = c(y = meth, x1 = ""), m = 1, censor = 6)
    expect_true(all(unlist(imp$imp$y) <= 6), info = meth)
    expect_false(anyNA(unlist(imp$imp$y)), info = meth)
  }
})


## --- 5: refusals, all of them before the retry loop -------------------------

test_that("B83: contradictory or absent censoring specifications are errors", {
  d <- .b83_daten(n = 200L, nmis = 40L)[, c("y", "x1")]

  expect_error(imp_c(d, method = c(y = "cp", x1 = ""), m = 1),
               "'censor' is required")
  ## a limit below the observed maximum: the data contradict it
  expect_error(imp_c(d, method = c(y = "cp", x1 = ""), m = 1, censor = 3),
               "exceed their censoring limit")
  expect_error(imp_c(d, method = c(y = "cp", x1 = ""), m = 1, censor = 6.5),
               "whole numbers")
  expect_error(imp_c(d, method = c(y = "cp", x1 = ""), m = 1,
                     censor = rep(6, 3)), "length 3")
  expect_error(imp_c(d, method = c(y = "cp", x1 = ""), m = 1, censor = 0),
               "below 1")
  expect_error(imp_c(d, method = c(y = "cp", x1 = ""), m = 1,
                     censor = c(6, NA, rep(6, nrow(d) - 2L))), "must not contain NA")
})

test_that("B83: latent = FALSE with censored cells in `where` is refused", {
  d <- .b83_daten(n = 200L, nmis = 40L)[, c("y", "x1")]
  wm <- matrix(FALSE, nrow(d), ncol(d), dimnames = dimnames(d))
  wm[, "y"] <- is.na(d$y) | (!is.na(d$y) & d$y >= 6)
  expect_error(imp_c(d, method = c(y = "cp", x1 = ""), m = 1, where = wm,
                     censor = 6, latent = FALSE),
               "censored observation")
})

test_that("B83: data with nothing but censored observations are refused", {
  set.seed(834); n <- 60
  X <- cbind(1, rnorm(n))
  expect_error(ci(".countimp_cs_fit")(X, rep(4, n), rep(4, n), "poisson"),
               "no information")
})

test_that("B83: `censor` may be given per variable, and a gap is an error", {
  set.seed(835); n <- 250
  x1 <- rnorm(n)
  d <- data.frame(a = pmin(rpois(n, exp(1.2 + 0.4 * x1)), 5),
                  b = pmin(rpois(n, exp(1.0 - 0.3 * x1)), 4),
                  x1 = x1)
  d$a[sample.int(n, 40)] <- NA
  d$b[sample.int(n, 40)] <- NA
  imp <- imp_c(d, method = c(a = "cp", b = "cp", x1 = ""), m = 1,
               censor = list(a = 5, b = 4))
  expect_true(all(unlist(imp$imp$a) <= 5))
  expect_true(all(unlist(imp$imp$b) <= 4))
  expect_error(imp_c(d, method = c(a = "cp", b = "cp", x1 = ""), m = 1,
                     censor = list(a = 5)),
               "does not cover variable 'b'")
})


test_that("B83: two scale families side by side in one call", {
  ## A call may mix them: one variable top-coded, another bounded by its scale.
  ## Both arguments then reach BOTH methods, and the per-variable check used to
  ## fire on the method that does not read the argument at all -- "'bounds'
  ## does not cover variable 'y', which is imputed by a bounded method", on a
  ## `y` imputed by `cp`. Complete specifications, aborted run.
  set.seed(837); n <- 300
  x1 <- rnorm(n)
  d <- data.frame(y = pmin(rpois(n, exp(1.4 + 0.5 * x1)), 8),
                  z = pmin(rpois(n, exp(1.0 - 0.3 * x1)), 5), x1 = x1)
  d$y[sample.int(n, 60)] <- NA
  d$z[sample.int(n, 60)] <- NA

  imp <- imp_c(d, method = c(y = "cp", z = "bp", x1 = ""), m = 2, maxit = 1,
               censor = list(y = 8), bounds = list(z = c(0, 5)))
  expect_true(all(unlist(imp$imp$y) <= 8))
  expect_true(all(unlist(imp$imp$z) <= 5))
  expect_true(all(unlist(imp$imp$z) >= 0))

  ## and the B58 rule still holds for the method that DOES read the argument:
  ## a list that omits a censored variable is an error, not a silent fallback
  expect_error(imp_c(d, method = c(y = "cp", z = "bp", x1 = ""), m = 1,
                     maxit = 1, censor = list(z = 8), bounds = list(z = c(0, 5))),
               "does not cover variable 'y'")
})


## --- 6: the internal-fit contract (B77) -------------------------------------

test_that("B83: the censored fit keeps the internal fit contract", {
  set.seed(836); n <- 200
  X <- cbind(1, rnorm(n))
  y <- pmin(rpois(n, exp(1.3 + 0.4 * X[, 2])), 5)
  f <- ci(".countimp_cs_fit")(X, y, rep(5, n), "negbin")
  expect_true(all(c("beta", "cov", "scale", "theta", "ll", "conv", "nobs",
                    "npar") %in% names(f)))
  expect_null(attr(f, "class"))
  expect_identical(names(f$theta), NULL)
  expect_equal(f$nobs, n)
  expect_equal(f$npar, 3L)
})


## --- B93: latent = TRUE without censored cells in `where` -------------------
##
## Reported from the simulation side, where it cost weeks. Study 10 called
## latent = TRUE but passed no `where`, so only the missing cells were drawn on
## the latent scale while the censored observations stayed at their limit --
## one column, two scales. Coverage against the latent truth was 6.6 % instead
## of 91.7 %, and nothing in the run said why.
##
## The package behaved as documented: the AUTOMATIC choice is recorded. But the
## choice that can go wrong is the one passed by argument, and that one was
## recorded nowhere. Note the asymmetry that made it worse: the mirror-image
## combination (where = censored plus latent = FALSE) has always been refused
## outright.

test_that("B93: latent = TRUE without censored cells in `where` is reported", {
  set.seed(931); n <- 400L
  x <- stats::rnorm(n)
  y_lat <- stats::rpois(n, exp(1.2 + 0.4 * x))
  limit <- 4
  d <- data.frame(y = pmin(y_lat, limit), x = x)
  d$y[sample.int(n, 100L)] <- NA
  stopifnot(sum(!is.na(d$y) & d$y >= limit) > 0)   # censoring must bind

  countimp_diagnostics(enable = TRUE)
  countimp_diagnostics(reset = TRUE)
  ## The warning fires once per session, so the flag has to be cleared for the
  ## test to see it regardless of what ran before.
  ##
  ## Through a local variable, not as `ci(".countimp_state")$flag <- FALSE`:
  ## ci() returns the value, and R cannot assign into the result of a call.
  ## That line aborted the block before its first expectation -- and the silent
  ## reporter counts an aborted block as zero checks, not as a failure, so it
  ## looked green here while R CMD check reported an error. `.countimp_state`
  ## is an environment, so assigning through a local name reaches the original.
  state <- ci(".countimp_state")
  state$mixed_scale_warned <- FALSE

  expect_warning(
    countimp(d, method = c(y = "cp", x = ""), m = 1, maxit = 1, seed = 1,
             printFlag = FALSE, censor = limit, latent = TRUE),
    "two scales")

  dg <- countimp_diagnostics()
  expect_true(any(grepl("censor_scale", dg$method)))
  expect_match(paste(dg$problems, collapse = " "), "no censored cell is marked")
})

test_that("B93: the combination that is not contradictory stays quiet", {
  ## latent = TRUE WITH censored cells in `where` merely confirms the automatic
  ## choice. Warning there would train users to ignore the warning.
  set.seed(932); n <- 400L
  x <- stats::rnorm(n)
  y_lat <- stats::rpois(n, exp(1.2 + 0.4 * x))
  limit <- 4
  d <- data.frame(y = pmin(y_lat, limit), x = x)
  d$y[sample.int(n, 100L)] <- NA
  censored <- !is.na(d$y) & d$y >= limit
  wh <- matrix(FALSE, n, 2L, dimnames = list(NULL, c("y", "x")))
  wh[, "y"] <- is.na(d$y) | censored

  state <- ci(".countimp_state")
  state$mixed_scale_warned <- FALSE
  expect_no_warning(
    countimp(d, method = c(y = "cp", x = ""), m = 1, maxit = 1, seed = 1,
             printFlag = FALSE, censor = limit, latent = TRUE, where = wh))
})

test_that("B93: the report does not change the imputations", {
  ## A log entry and a warning must not touch the draw. Same seed, same values.
  set.seed(933); n <- 300L
  x <- stats::rnorm(n)
  d <- data.frame(y = pmin(stats::rpois(n, exp(1.2 + 0.4 * x)), 4), x = x)
  d$y[sample.int(n, 80L)] <- NA
  draw_ <- function() suppressWarnings(unlist(
    countimp(d, method = c(y = "cp", x = ""), m = 2, maxit = 2, seed = 77,
             printFlag = FALSE, censor = 4, latent = TRUE)$imp$y))
  expect_identical(draw_(), draw_())
})
