## B84: two-level Conway-Maxwell-Poisson (2l.cmp and its three variants)
##
## What these tests defend, in the order in which each was easy to get wrong:
##
##  1. THE SHAPE ORIENTATION. glmmTMB reports 1/nu through sigma(). If a future
##     release flips that, imputations would be drawn at inverted dispersion and
##     still look like plausible counts -- the failure mode the register calls
##     the dangerous one. countimp measures the orientation instead of believing
##     it, and the mutation probe below simulates exactly that flip with a spy.
##  2. THE DISPERSION FLOOR. .countimp_draw_2l() floors the drawn dispersion at
##     0.25 to guard the NB tail (B09a). For compois the reported dispersion is
##     1/nu, so that floor caps nu at 4 -- and an ordinary underdispersed fit
##     sits below it (sigma = 0.172 for the data below, i.e. nu = 5.8 capped to
##     4; analyse/k25_2lcmp_boden.csv). The floor is therefore switched off for
##     this family, and the test shows both halves: that it WOULD have bound,
##     and that it does not.
##  3. The draw comes from countimp's sampler, not from glmmTMB::simulate(),
##     which ignores newdata and returns values for the fitting cases. A return
##     length equal to sum(wy) rather than nrow(fit data) is what shows it.
##  4. The imputations actually carry the underdispersion. Poisson cannot go
##     below a ratio of one and NB can only go above it, so both are the wrong
##     answer here: measured at nu ~ 10, the dispersion of the imputed cells was
##     0.092 for 2l.cmp against 1.106 for 2l.poisson and 1.069 for 2l.nb, on
##     data whose own dispersion is 0.098 (analyse/k25_2lcmp_dispersion.csv).

skip_if_not_installed("glmmTMB")

## Small on purpose: a compois fit costs 1-3 s, and R CMD check runs this.
.b84_daten <- function(ngrp = 12L, nj = 15L, seed = 84L, nmis = 60L,
                       groesse = 12L, p0 = 1.6, b = 0.25, tau = 0.15) {
  set.seed(seed)
  n <- ngrp * nj
  id <- rep(seq_len(ngrp), each = nj)
  u <- stats::rnorm(ngrp, 0, tau)
  x1 <- stats::rnorm(n)
  y <- stats::rbinom(n, groesse, stats::plogis(p0 + b * x1 + u[id]))
  d <- data.frame(y = y, x1 = x1, id = id)
  d$wahr <- y
  d$y[sample.int(n, nmis)] <- NA
  d
}

.b84_pred <- function() {
  p <- matrix(0L, 3L, 3L,
              dimnames = list(c("y", "x1", "id"), c("y", "x1", "id")))
  p["y", ] <- c(0L, 1L, -2L)
  p
}

## Conditional Pearson dispersion: Var(Y|x) / E(Y|x), the quantity the count
## families differ in. The marginal ratio is not it -- see the note in
## R/fitdiag.R, where a marginally overdispersed sample was conditionally
## underdispersed throughout.
.b84_disp <- function(dat) {
  g <- suppressWarnings(stats::glm(y ~ x1 + factor(id),
                                   family = stats::poisson(), data = dat))
  sum(stats::residuals(g, "pearson")^2) / g$df.residual
}

## Pearson dispersion OF THE IMPUTED CELLS, measured against a fit on the
## OBSERVED ones. This is the sharp version of the question "do the imputations
## carry the underdispersion": a Poisson draw lands at 1 whatever the fitted
## means are, a COM-Poisson draw at roughly 1/nu.
##
## The blunt version -- comparing the dispersion of the completed data against
## 2l.poisson's -- was measured and does NOT separate the two: with the mean
## structure taken from the compois fit, a Poisson draw substituted for the
## COM-Poisson one passed it (mutation probe M3, 22 August 2026). The imputed
## cells are only a quarter of the data, so their contribution is diluted by
## the observed ones.
.b84_disp_imp <- function(d, imp, m) {
  o <- d[!is.na(d$y), ]
  g <- suppressWarnings(stats::glm(y ~ x1 + factor(id),
                                   family = stats::poisson(), data = o))
  nd <- data.frame(x1 = d$x1[is.na(d$y)],
                   id = factor(d$id[is.na(d$y)], levels = levels(factor(o$id))))
  muh <- suppressWarnings(stats::predict(g, newdata = nd, type = "response"))
  ok <- is.finite(muh) & muh > 0
  mean(vapply(seq_len(m), function(i)
    mean((imp$imp$y[[i]][ok] - muh[ok])^2 / muh[ok]), numeric(1)))
}

.b84_fit <- function(d) {
  o <- d[!is.na(d$y), ]
  glmmTMB::glmmTMB(y ~ x1 + (1 | id),
                   data = data.frame(y = o$y, x1 = o$x1, id = factor(o$id)),
                   family = glmmTMB::compois())
}


## --- 1: the orientation is measured, not assumed ---------------------------

test_that("B84: the shape is read as 1/sigma, and the reading is reported", {
  d <- .b84_daten()
  f <- .b84_fit(d)
  o <- d[!is.na(d$y), ]
  mu <- as.numeric(stats::predict(f, type = "response"))
  sh <- ci(".countimp_cmp_nu_glmm")(f, o$y, mu)

  s <- as.numeric(stats::sigma(f))
  expect_true(sh$rezi)
  expect_equal(sh$nu, 1 / s, tolerance = 1e-8)
  ## and it is the reading the data support: the alternative fits far worse
  ll <- ci(".countimp_cmp_ll_mu")
  expect_gt(ll(o$y, mu, 1 / s), ll(o$y, mu, s))
  ## the measured shape is near the reciprocal of the conditional dispersion,
  ## which is what Var ~ mu/nu says it should be
  expect_lt(abs(sh$nu - 1 / .b84_disp(d[!is.na(d$y), ])), 2)
})

test_that("B84: a flipped sigma convention is detected rather than believed", {
  ## Mutation probe. The spy makes .countimp_theta() return nu where glmmTMB
  ## returns 1/nu -- exactly what a future release could do. The rule must then
  ## pick the OTHER reading, and say so in the diagnostics log.
  d <- .b84_daten()
  f <- .b84_fit(d)
  o <- d[!is.na(d$y), ]
  mu <- as.numeric(stats::predict(f, type = "response"))
  wahr <- ci(".countimp_cmp_nu_glmm")(f, o$y, mu)$nu

  countimp_diagnostics(enable = TRUE, reset = TRUE)
  ci_spy(".countimp_theta", function(fit, summ = NULL, pflicht = TRUE)
    1 / as.numeric(stats::sigma(fit)))
  sh <- ci(".countimp_cmp_nu_glmm")(f, o$y, mu)

  expect_false(sh$rezi)                       # the flip was noticed
  expect_equal(sh$nu, wahr, tolerance = 1e-6) # and the shape is still right
  log <- countimp_diagnostics()
  expect_true(any(grepl("compois_sigma_orientation", log$method)))
  countimp_diagnostics(reset = TRUE)
})


## --- 2: the dispersion floor is off for this family ------------------------

test_that("B84: the NB dispersion floor would bind here, and is switched off", {
  d <- .b84_daten()
  f <- .b84_fit(d)
  ## An ordinary underdispersed fit reports sigma below the NB floor: the floor
  ## is not a remote corner case for this family, it is the normal case.
  expect_lt(as.numeric(stats::sigma(f)), 0.25)

  nd <- data.frame(x1 = d$x1[is.na(d$y)], id = factor(d$id[is.na(d$y)]))
  lev <- levels(factor(d$id[!is.na(d$y)]))
  set.seed(1)
  frei <- replicate(5, ci(".countimp_draw_2l")(f, nd, grp = "id",
                                               obs_levels = lev,
                                               theta.min = 0)$theta)
  set.seed(1)
  boden <- replicate(5, ci(".countimp_draw_2l")(f, nd, grp = "id",
                                                obs_levels = lev,
                                                theta.min = 0.25)$theta)
  expect_true(all(boden >= 0.25))            # the floor does what it says
  expect_lt(stats::median(frei), 0.25)       # and 2l.cmp is not subject to it
})


## --- 3: the draw is countimp's, and it is a draw for the right rows --------

test_that("B84: imputations are returned for the imputed cases, not the fitted", {
  d <- .b84_daten(nmis = 60L)
  imp <- suppressWarnings(countimp(d[, c("y", "x1", "id")],
                                   method = c(y = "2l.cmp", x1 = "", id = ""),
                                   predictorMatrix = .b84_pred(), m = 1,
                                   maxit = 1, printFlag = FALSE))
  v <- imp$imp$y[[1L]]
  expect_equal(length(v), sum(is.na(d$y)))   # not nrow of the fitting data
  expect_false(anyNA(v))
  expect_true(all(v >= 0))
  expect_true(all(v == round(v)))
})

test_that("B84: all four exported variants run and stay in range", {
  d <- .b84_daten(nmis = 45L)[, c("y", "x1", "id")]
  ## .noint needs a random slope: without one the model has an empty random
  ## part, which countimp refuses on purpose (passing (0 | id) to glmmTMB
  ## crashes the session). So the noint variants get a slope in `type`.
  p_noint <- .b84_pred()
  p_noint["y", ] <- c(0L, 2L, -2L)
  for (meth in c("2l.cmp", "2l.cmp.boot")) {
    imp <- suppressWarnings(countimp(d, method = c(y = meth, x1 = "", id = ""),
                                     predictorMatrix = .b84_pred(), m = 1,
                                     maxit = 1, printFlag = FALSE))
    expect_false(anyNA(unlist(imp$imp$y)), info = meth)
  }
  for (meth in c("2l.cmp.noint", "2l.cmp.noint.boot")) {
    imp <- suppressWarnings(countimp(d, method = c(y = meth, x1 = "", id = ""),
                                     predictorMatrix = p_noint, m = 1,
                                     maxit = 1, printFlag = FALSE))
    expect_false(anyNA(unlist(imp$imp$y)), info = meth)
  }
})


## --- 4: the imputations carry the underdispersion --------------------------

test_that("B84: 2l.cmp keeps the conditional dispersion, Poisson cannot", {
  ## Strongly underdispersed: 12 trials at p ~ 0.9 gives Var/E ~ 0.1, i.e.
  ## nu ~ 10. Chosen so that the three regimes are far apart and the thresholds
  ## below are not deciding on noise:
  ##
  ##   COM-Poisson draw     0.092     (measured, m = 3)
  ##   nu floored at 4      ~0.25     (what theta.min = 0.25 would produce)
  ##   Poisson draw         1.106     (measured)
  ## Reproduced by analyse/k25_2lcmp_belege.R, which uses this same design.
  d <- .b84_daten(ngrp = 15L, nj = 16L, nmis = 90L, groesse = 12L,
                  p0 = 2.2, b = 0.3, tau = 0.2, seed = 841L)
  wahr <- .b84_disp(data.frame(y = d$wahr, x1 = d$x1, id = d$id))
  expect_lt(wahr, 0.15)                       # the data really are underdispersed

  dd <- d[, c("y", "x1", "id")]
  hol <- function(meth) {
    set.seed(4)
    imp <- suppressWarnings(countimp(dd, method = c(y = meth, x1 = "", id = ""),
                                     predictorMatrix = .b84_pred(), m = 3,
                                     maxit = 1, printFlag = FALSE))
    .b84_disp_imp(d, imp, 3)
  }
  d_cmp <- hol("2l.cmp")
  d_poi <- hol("2l.poisson")

  ## below the floor's own value: this is what fails if theta.min is ever put
  ## back for this family, and it is why the threshold is 0.20 and not 0.8
  expect_lt(d_cmp, 0.20)
  expect_gt(d_cmp, 0.03)                      # and not collapsed to a point mass
  expect_gt(d_poi, 0.80)                      # Poisson cannot go below one
  expect_lt(abs(d_cmp - wahr), abs(d_poi - wahr))
})
