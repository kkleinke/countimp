## B30: the two-level zero-inflated methods fitted a hurdle and drew as if the
## count part were unrestricted.
##
## Up to 2.6.0, mice.impute.2l.zip / .zinb ran a binomial GLMM on 1{y = 0},
## then an UNTRUNCATED count GLMM on the positive subset, and drew a count for
## every case the binomial part assigned to the count process. Two errors
## compound: the count model ignores the truncation of its own data, and a
## count-process case is still imputed as zero with probability
## (theta/(theta+mu))^theta. Both push the same way -- too many zeros,
## conditional mean too small. Measured coverage of the nominal 95 % interval
## for the marginal covariate effect: 0.885 under a true zero-inflation process
## and 0.870 under a true hurdle process (R = 200; see
## skripte/b30_zi_comparison.R). Since 3.0.0 these methods fit a genuine
## zero-inflation model, which restores coverage to 0.940 / 0.915 for zinb /
## zip on the package code (skripte/b30_paket_validierung.R).
##
## These tests guard the three properties that made the old code wrong, so a
## future refactor cannot quietly reintroduce any of them.

zi_data <- function(n = 500L, ngrp = 20L, seed = 1L, p.zero = 0.35) {
  set.seed(seed)
  nj <- n %/% ngrp
  g <- rep(seq_len(ngrp), each = nj)
  u <- stats::rnorm(ngrp, 0, sqrt(0.3))[g]
  x1 <- stats::rnorm(length(g)); x2 <- stats::rnorm(length(g))
  zi <- stats::rbinom(length(g), 1L, p.zero) == 1L
  y <- ifelse(zi, 0L, stats::rnbinom(length(g), size = 2,
                                     mu = exp(1 + 0.6 * x1 + 0.3 * x2 + u)))
  list(y = y, x = data.frame(g = g, x1 = x1, x2 = x2),
       type = stats::setNames(c(-2L, 1L, 1L), c("g", "x1", "x2")))
}

test_that("the ZI methods fit one model with a zero-inflation part, not a hurdle", {
  skip_if_not_installed("glmmTMB")
  ## The engine must pass a ziformula. A hurdle fit (two separate models, the
  ## count one on the positive subset) is what B30 was about.
  src <- deparse(ci(".countimp_2l_zi"))
  expect_true(any(grepl("ziformula", src, fixed = TRUE)))
  ## and it must NOT subset the data to the positives, which the old code did
  expect_false(any(grepl("Y > 0", src, fixed = TRUE)))
})

test_that("imputations are drawn, not deterministic, and are valid counts", {
  skip_if_not_installed("glmmTMB")
  d <- zi_data()
  ry <- rep(TRUE, length(d$y)); set.seed(5); ry[sample(length(d$y), 150L)] <- FALSE

  set.seed(1); a <- suppressWarnings(
    mice.impute.2l.zinb(y = d$y, ry = ry, x = d$x, type = d$type))
  set.seed(2); b <- suppressWarnings(
    mice.impute.2l.zinb(y = d$y, ry = ry, x = d$x, type = d$type))

  expect_length(a, sum(!ry))
  expect_true(all(is.finite(a)))
  expect_true(all(a >= 0))
  expect_equal(a, floor(a))
  expect_false(identical(a, b))          # two draws must differ
})

test_that("the zero proportion of the imputations tracks the observed one", {
  skip_if_not_installed("glmmTMB")
  skip_on_cran()
  ## The old code produced too many zeros: 0.531 imputed against 0.488
  ## observed in the reference run. This is the property that broke coverage.
  ## Averaged over draws the imputed proportion must sit near the observed one.
  d <- zi_data(n = 600L, ngrp = 24L, seed = 3L)
  ry <- rep(TRUE, length(d$y)); set.seed(7); ry[sample(length(d$y), 200L)] <- FALSE
  p.obs <- mean(d$y[ry] == 0)

  set.seed(11)
  p.imp <- mean(replicate(10L, {
    im <- suppressWarnings(mice.impute.2l.zip(y = d$y, ry = ry, x = d$x,
                                             type = d$type))
    mean(im == 0)
  }))
  ## generous band: this guards a systematic excess of zeros, not Monte Carlo noise
  expect_lt(abs(p.imp - p.obs), 0.10)
})

test_that("intercept.c and intercept.z act on separate model parts", {
  skip_if_not_installed("glmmTMB")
  ## Before 2.6.0 intercept.c governed both parts, so .noint.zero silently
  ## removed the COUNT intercept as well. Check the generated formulas.
  type <- stats::setNames(c(-2L, 2L, 4L), c("g", "x1", "x2"))
  dec <- ci(".countimp_decode_type")(type, c("g", "x1", "x2"))
  f <- function(part, ic) deparse(
    ci(".countimp_2l_formula")(dec, part, "Y", intercept = ic)[[3L]])

  expect_match(f("count", TRUE),  "(1 + x1 + x2 | g)", fixed = TRUE)
  expect_match(f("count", FALSE), "(0 + x1 + x2 | g)", fixed = TRUE)
  expect_match(f("zero",  TRUE),  "(1 + x1 | g)", fixed = TRUE)
  expect_match(f("zero",  FALSE), "(0 + x1 | g)", fixed = TRUE)
})

test_that("type codes 3-6 assign predictors to the intended part only", {
  dec <- ci(".countimp_decode_type")(
    c(-2L, 1L, 3L, 5L, 4L, 6L), c("g", "shared", "cnt", "zro", "cntR", "zroR"))
  expect_setequal(dec$c.fixed, c("shared", "cnt", "cntR"))
  expect_setequal(dec$z.fixed, c("shared", "zro", "zroR"))
  expect_setequal(dec$c.slope, "cntR")
  expect_setequal(dec$z.slope, "zroR")
  expect_identical(dec$group, "g")
})

test_that("a missing class variable is reported; several are grouping levels", {
  expect_error(ci(".countimp_decode_type")(c(1L, 1L), c("a", "b")),
               "no class variable", fixed = TRUE)
  ## Two -2 entries were an error up to 3.0.0 ("only one class allowed!") and
  ## are the three-level case since: (1 | school) + (1 | class) compiles to two
  ## -2 codes, in column order, and group[1] is the level that carries the
  ## random slopes. Changed deliberately -- this assertion is what caught it.
  dec <- ci(".countimp_decode_type")(c(-2L, -2L), c("a", "b"))
  expect_identical(dec$group, c("a", "b"))
  expect_error(ci(".countimp_decode_type")(c(-2L, 1L), c("a", "b", "c")),
               "length(type)", fixed = TRUE)
})

test_that("all 16 exported ZI variants exist and share one engine", {
  nms <- as.vector(outer(
    as.vector(outer(c("mice.impute.2l.zinb", "mice.impute.2l.zip"),
                    c("", ".noint.both", ".noint.count", ".noint.zero"), paste0)),
    c("", ".boot"), paste0))
  expect_length(nms, 16L)
  for (nm in nms) {
    f <- get(nm, envir = asNamespace("countimp"))
    expect_true(is.function(f), info = nm)
    ## each is a thin wrapper: the body must call the shared engine
    expect_true(any(grepl(".countimp_2l_zi", deparse(body(f)), fixed = TRUE)),
                info = nm)
  }
})
