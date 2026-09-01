## B94 -- theta is not drawn where the data do not determine it.
##
## Found while building rcountimp, in the compiled core first, and present here
## for the same reason: .countimp_draw_theta() drew from
## exp(rnorm(1, log(theta), se)) whatever se was.
##
## Fitting a negative binomial to data that carry no overdispersion -- the
## cautious choice, and therefore a frequent one -- sends theta towards
## infinity along a flat likelihood, and se(log theta) with it. The draw then
## spans orders of magnitude, and one near theta.min means Var = mu + 4 mu^2.
## Because every FCS cycle refits on the previous cycle's imputations, an
## unlucky draw feeds itself. Measured before the guard, on Poisson data with
## n = 800 and m = 10: one imputation reached sd 5.21 and a largest value of 44,
## against sd 1.92 and 11 in the data.
##
## .countimp_stabilize() covers this for the multilevel paths by fixing theta;
## the single-level path over MASS::glm.nb had nothing.

skip_if_not_installed("MASS")

test_that("B94: a theta with a huge standard error is not drawn", {
  ## The fit is built directly, so the test states the condition it is about
  ## rather than hoping some data set produces it.
  ## theta is passed explicitly: on a negbin object stats::sigma() returns the
  ## residual scale, not theta -- the function documents that, and the test
  ## must not fall into it.
  fit <- structure(list(theta = 8160, SE.theta = 8160 * 14.66),
                   class = c("negbin", "glm"))
  draw_ <- ci(".countimp_draw_theta")
  gezogen <- vapply(1:200, function(i) draw_(fit, theta = fit$theta, quiet = TRUE), numeric(1))
  expect_true(all(gezogen == 8160),
              info = "theta was drawn although se(log theta) is 14.66")
})

test_that("B94: a theta with an ordinary standard error is still drawn", {
  ## The guard must not swallow the identified case -- that is what the draw is
  ## for. Measured across MASS::glm.nb replications, se(log theta) stays below
  ## 0.3 for theta = 0.5 and 2, so 0.15 is an ordinary value.
  fit <- structure(list(theta = 2, SE.theta = 2 * 0.15),
                   class = c("negbin", "glm"))
  draw_ <- ci(".countimp_draw_theta")
  set.seed(94)
  gezogen <- vapply(1:200, function(i) draw_(fit, theta = fit$theta, quiet = TRUE), numeric(1))
  expect_gt(length(unique(gezogen)), 100)
  ## and on the log scale it scatters as it should
  expect_lt(abs(sd(log(gezogen)) - 0.15), 0.03)
  expect_lt(abs(mean(log(gezogen)) - log(2)), 0.05)
})

test_that("B94: the threshold sits between the two regimes", {
  ## Not a round number chosen for looks: no identified fit reaches it, and at
  ## se = 5 the 95 % range of the draw already spans four orders of magnitude.
  draw_ <- ci(".countimp_draw_theta")
  fest <- function(se) {
    fit <- structure(list(theta = 100, SE.theta = 100 * se),
                     class = c("negbin", "glm"))
    all(vapply(1:50, function(i) draw_(fit, theta = fit$theta, quiet = TRUE), numeric(1)) == 100)
  }
  set.seed(941)
  expect_false(fest(0.3))    # theta = 2 or 0.5, any n
  expect_false(fest(1.8))    # theta = 20, n = 200
  expect_true(fest(7.9))     # Poisson data, median
  expect_true(fest(14.66))   # the measured case
})

test_that("B94: imputations no longer explode on data without overdispersion", {
  set.seed(940); n <- 800L
  d <- data.frame(x = stats::rnorm(n))
  d$y <- stats::rpois(n, exp(0.7 + 0.5 * d$x))
  d$y[sample.int(n, 240L)] <- NA
  missing <- is.na(d$y)
  obs_sd <- stats::sd(d$y, na.rm = TRUE)

  imp <- suppressWarnings(countimp(d, method = c(y = "nb", x = ""), m = 10,
                                   maxit = 5, seed = 4711, printFlag = FALSE))
  sds <- vapply(seq_len(10), function(k)
    stats::sd(countimp_complete(imp, k)$y[missing]), numeric(1))

  ## Every imputation, not the average: the defect showed up in single passes
  ## while the others looked fine.
  expect_true(all(sds < 2 * obs_sd),
              info = paste("worst imputed sd", round(max(sds), 2),
                           "against", round(obs_sd, 2), "observed"))
})
