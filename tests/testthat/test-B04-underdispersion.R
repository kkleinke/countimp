## B04: quasi-Poisson imputation with phi < 1.
##
## The old code drew from a negative binomial with theta = mu/(phi - 1), which
## is negative when phi < 1 and yields NaN for every value. Underdispersed
## counts are common (bounded scales, binomial-like processes), so the failure
## was not exotic: it produced 270 NaN out of 270 imputations.

test_that("underdispersed counts impute to finite values", {
  d <- sim_underdispersed()
  im <- imputed_values(quiet_impute(d, method = c("quasipoisson", ""),
                                    m = 3, maxit = 3, seed = 5))
  expect_true(all(is.finite(im)))
  expect_true(all(im >= 0))
  expect_true(all(im == round(im)))
})

test_that("the underdispersion fallback warns exactly once per session", {
  ## Bind the environment to a local name first. `pkg:::env$x <- value` is not
  ## a legal complex assignment -- R rewrites it through `:::<-`, which fails
  ## with a packageNotFoundError that testthat reports as a skip rather than a
  ## failure, so the test would silently never run.
  st <- ci(".countimp_state")
  st$underdisp_warned <- FALSE
  d <- sim_underdispersed()
  w <- character(0)
  withCallingHandlers(
    countimp(d, method = c("quasipoisson", ""), m = 3, maxit = 3, seed = 5,
             printFlag = FALSE),
    warning = function(cond) {
      w <<- c(w, conditionMessage(cond)); invokeRestart("muffleWarning")
    })
  expect_equal(sum(grepl("underdispers", w)), 1L)
})

test_that(".countimp_rqpois covers all dispersion regimes", {
  set.seed(2)
  for (phi in c(0.43, 1.0, 1.0000001, 2.5)) {
    v <- suppressWarnings(ci(".countimp_rqpois")(rep(4, 2000), phi, quiet = TRUE))
    expect_true(all(is.finite(v)), info = paste("phi =", phi))
    expect_true(all(v >= 0), info = paste("phi =", phi))
  }
  ## overdispersion must still be overdispersed
  set.seed(2)
  v <- ci(".countimp_rqpois")(rep(4, 20000), 2.5, quiet = TRUE)
  expect_gt(var(v), 4)
})

test_that("non-finite input is an error, not a silent NA", {
  expect_error(ci(".countimp_rqpois")(c(4, Inf), 2))
  expect_error(ci(".countimp_rqpois")(c(4, 5), NA_real_))
})

test_that("bootstrap variant is wired identically", {
  d <- sim_underdispersed()
  im <- imputed_values(quiet_impute(d, method = c("quasipoisson.boot", ""),
                                    m = 2, maxit = 2, seed = 9))
  expect_true(all(is.finite(im)))
})


## The abort message must not blame the predictors -------------------------
##
## Reported from the simulation side on 28 August 2026, from study 5. `bnb`
## aborts on underdispersed data -- correctly, because a negative binomial has
## no finite theta there: it contains the Poisson only as the limit
## theta -> Inf, so for variance < mean the likelihood keeps rising. Measured
## with MASS::glm.nb on the same data (K = 8, var/mean = 0.657, N = 500):
## theta = 141167, SE = 885097, iteration limit reached.
##
## What was wrong is the message. It named collinear predictors as the usual
## cause and advised dropping one -- on two INDEPENDENT standard normal
## predictors. And a user with sum scores over K items has underdispersed data
## always, so the advice would send them to delete predictors that are fine.
##
## The fit still fails; nothing was made to work that cannot work. The
## diagnosis is added where it applies, and the misleading half of the inner
## message is dropped with it.

test_that("B04: an underdispersed NB failure names the dispersion, not the predictors", {
  ## Driven through .countimp_draw_retry() directly: the message is what is
  ## under test, and a real fit would make the test depend on which of several
  ## failures happens to come first.
  scheitert <- function() stop(structure(
    class = c("countimp_error", "error", "condition"),
    list(message = paste("the coefficient covariance matrix is not positive",
                         "definite, so no coefficients can be drawn.",
                         "Collinear predictors are the usual cause; drop one",
                         "of them or use ridge."), call = NULL)))
  set.seed(4)
  y.unter <- rbinom(400L, 8L, 0.5)            # variance/mean well below 1
  vm <- stats::var(y.unter) / mean(y.unter)
  expect_lt(vm, 1)

  err <- tryCatch(ci(".countimp_draw_retry")(scheitert, y_obs = y.unter,
                    method = "negbin"), error = conditionMessage)
  expect_match(err, "UNDERDISPERSED", fixed = TRUE)
  expect_match(err, sprintf("variance/mean = %.2f", vm), fixed = TRUE)
  ## names the way out, and names it as a method
  expect_match(err, '"cmp"', fixed = TRUE)
  expect_match(err, '"bp"', fixed = TRUE)
  ## and the advice that would send the user to delete good predictors is gone
  expect_no_match(err, "Collinear predictors", fixed = TRUE)
  expect_no_match(err, "drop collinear", fixed = TRUE)
})

test_that("B04: an OVERdispersed failure keeps the general advice", {
  ## The guard against over-correcting: underdispersion and a genuinely
  ## collinear design can occur together, and where dispersion is not the
  ## cause the message must stay as it was.
  scheitert <- function() stop(structure(
    class = c("countimp_error", "error", "condition"),
    list(message = "the coefficient covariance matrix is not positive definite.",
         call = NULL)))
  set.seed(4)
  y.ueber <- rnbinom(400L, mu = 3, size = 0.8)
  expect_gt(stats::var(y.ueber) / mean(y.ueber), 1)
  err <- tryCatch(ci(".countimp_draw_retry")(scheitert, y_obs = y.ueber,
                    method = "negbin"), error = conditionMessage)
  expect_no_match(err, "UNDERDISPERSED", fixed = TRUE)
  expect_match(err, "rescale or drop collinear predictors", fixed = TRUE)

  ## The boundary is variance/mean = 1, and it is not arbitrary: that is where
  ## the negative binomial degenerates to the Poisson. Data just ABOVE it are
  ## overdispersed, however slightly, and a finite theta exists -- so the
  ## diagnosis must not fire. Without this case the threshold could drift
  ## upwards unnoticed (a mutation probe raising it from 1 to 2 went unseen).
  set.seed(4)
  y.knapp <- rnbinom(400L, mu = 4, size = 30)
  vmk <- stats::var(y.knapp) / mean(y.knapp)
  expect_gt(vmk, 1); expect_lt(vmk, 1.5)
  errk <- tryCatch(ci(".countimp_draw_retry")(scheitert, y_obs = y.knapp,
                     method = "negbin"), error = conditionMessage)
  expect_no_match(errk, "UNDERDISPERSED", fixed = TRUE)
})

test_that("B04: a Poisson failure is not blamed on underdispersion", {
  ## A Poisson model has no theta to run away, so underdispersion is not its
  ## problem -- saying so would be a different wrong diagnosis.
  scheitert <- function() stop(structure(
    class = c("countimp_error", "error", "condition"),
    list(message = "the coefficient covariance matrix is not positive definite.",
         call = NULL)))
  set.seed(4)
  y.unter <- rbinom(400L, 8L, 0.5)
  err <- tryCatch(ci(".countimp_draw_retry")(scheitert, y_obs = y.unter,
                    method = "poisson"), error = conditionMessage)
  expect_no_match(err, "UNDERDISPERSED", fixed = TRUE)
  ## the censored NB does carry a theta, so it does get the diagnosis
  err2 <- tryCatch(ci(".countimp_draw_retry")(scheitert, y_obs = y.unter,
                     method = "cnb"), error = conditionMessage)
  expect_match(err2, "UNDERDISPERSED", fixed = TRUE)
})
