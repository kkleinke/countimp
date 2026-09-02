## B95 -- a draw far above the observed values is reported.
##
## Found while building the conference example. On the delinquency data
## (crim4w) `nb` imputes values up to 955 where the largest observed one is 16
## -- and the diagnostics said nothing, although they compute the ratio.
##
## The cause was not the family but the link: count variables used as
## predictors in a log-link model extrapolate exponentially. At the observed
## maxima of BCRIM and DCRIM, eta reaches 7.2, so mu = 1380 -- against a
## largest observed value of 16. Drawing beta widens this further.
##
## Nothing is rejected: large values are the tail of the distribution, and
## discarding them until a small one arrives would pull the imputations
## downwards. It is reported so that somebody sees it.
##
## The threshold of 10 is calibrated. Across families and data situations,
## each with a model that fits, the ratio lies between 1.0 and 3.7:
##
##   Poisson data, poisson  1.3      NB data (theta=2), nb     1.3
##   Poisson data, nb       2.6      NB data (theta=0.5), nb   1.0
##   ZI data, zip           1.7      ZI data, nb               3.7

test_that("B95: a draw far above the data is flagged", {
  check <- ci(".countimp_check_draws")
  ## The ratio at issue: largest imputation to largest observed value.
  p <- check(c(1, 2, 955), c(0, 1, 2, 16))
  expect_true("draw_far_above_observed" %in% as.character(p))
  expect_equal(attr(p, "draw_ratio"), 955 / 16)
})

test_that("B95: ordinary draws pass without comment", {
  ## Otherwise the report would be worthless: a warning that always fires is
  ## read past. The ratio 3.7 is the highest value measured in a healthy run.
  check <- ci(".countimp_check_draws")
  for (r in c(1.0, 1.3, 2.6, 3.7, 9.9)) {
    p <- check(c(1, 2, r * 16), c(0, 1, 2, 16))
    expect_false("draw_far_above_observed" %in% as.character(p),
                 info = paste("ratio", r, "was flagged"))
  }
})

test_that("B95: the report is a flag, not a hard failure", {
  ## "flagged" means, in this package: usable but unusual -- recorded,
  ## reported and KEPT. A hard failure would repeat the draw, which would be
  ## wrong here: one would keep drawing until a small value turned up.
  ## action = "silent": record and return, do not warn. ("none" is not one of
  ## the choices -- this line used to say so, and the block then aborted
  ## before its first expectation. testthat's silent reporter counts an
  ## aborted block as zero checks rather than as a failure, so it looked
  ## green; R CMD check reports it as an error. A block that runs no
  ## expectation is a broken test, not a passing one.)
  z <- countimp_check(imp = c(1, 2, 955), y_obs = c(0, 1, 2, 16),
                      method = "nb", action = "silent")
  expect_identical(z$status, "flagged")
  expect_match(z$problems, "draw_far_above_observed")
})

test_that("B95: on the data where it surfaced", {
  skip_if_not_installed("MASS")
  ## The single-level zero-inflated and hurdle methods are fitted through
  ## pscl, which is in Suggests. Without this the block fails on the missing
  ## dependency instead of on the property under test -- as it did on the
  ## no-mice CI job, which installs hard dependencies only.
  skip_if_not_installed("pscl")
  data(crim4w)
  countimp_diagnostics(enable = TRUE)
  countimp_diagnostics(reset = TRUE)
  invisible(suppressWarnings(countimp(crim4w,
    formulas = list(CCRIM = CCRIM ~ FEMALE + GY + ACRIM + BCRIM + DCRIM),
    family = list(CCRIM = nb()), m = 5, maxit = 5, seed = 2026,
    printFlag = FALSE)))
  dg <- countimp_diagnostics()
  flagged <- dg[dg$method == "negbin" & dg$status == "flagged", ]
  expect_gt(nrow(flagged), 0)
  expect_match(paste(flagged$problems, collapse = " "), "draw_far_above_observed")

  ## And with a family that fits these data, it stays quiet.
  countimp_diagnostics(reset = TRUE)
  invisible(suppressWarnings(countimp(crim4w,
    formulas = list(CCRIM = CCRIM ~ FEMALE + GY + ACRIM + BCRIM + DCRIM),
    family = list(CCRIM = zi_poisson()), m = 5, maxit = 5, seed = 2026,
    printFlag = FALSE)))
  dg2 <- countimp_diagnostics()
  zi <- dg2[dg2$method == "zeroinfl.poisson", ]
  expect_true(all(zi$status == "ok"))
})
