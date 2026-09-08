## B97: the formula route must split the model parts for the SINGLE-level
## two-part families too.
##
## The defect this guards against was silent. `zero = ~ z1` with zi_poisson()
## produced Y ~ 1 | x1 + x2: the count part lost every predictor, they moved
## into the zero part, and z1 was dropped -- no error, no warning, and the
## method name stayed correct at "zip". Only the fitted model showed it.
##
## Hence these tests look at the model that gets built, never at the method
## name. test-B53-spec.R checks the same call and cannot see the defect,
## because m1() reads only s$method[[y]].

zeug_b97 <- function(n = 200L, saat = 20260908L) {
  set.seed(saat)
  x1 <- stats::rnorm(n); x2 <- stats::rnorm(n); z1 <- stats::rnorm(n)
  y <- stats::rpois(n, exp(1 + .3 * x1 + .3 * x2))
  y[stats::runif(n) < stats::plogis(2 * z1)] <- 0L
  d <- data.frame(y, x1, x2, z1)
  d$y[stats::runif(n) < .2] <- NA
  d
}

test_that("B97: `zero` separates the parts for every single-level two-part family", {
  skip_if_not_installed("pscl")
  d <- zeug_b97()
  for (nm in c("zi_poisson", "zi_nb", "hurdle_poisson", "hurdle_nb")) {
    fam <- get(nm, envir = asNamespace("countimp"))()
    im <- suppressWarnings(countimp(d, formulas = list(y ~ x1 + x2),
            zero = ~ z1, family = fam, m = 1L, maxit = 1L, printFlag = FALSE))
    ty <- im$predictorMatrix["y", ]
    ty <- ty[names(ty) != "y"]
    ## The count predictors carry mice's code 2, the zero predictor code 3.
    expect_identical(unname(ty[c("x1", "x2")]), c(2, 2), info = nm)
    expect_identical(unname(ty["z1"]), 3, info = nm)
    ## And the model the method builds from them -- the size that matters.
    f <- ci(".countimp_zi_formula")(names(ty), ty)
    expect_identical(deparse(f), deparse(Y ~ x1 + x2 | z1), info = nm)
  }
})

test_that("B97: the two-level route keeps its own coding", {
  skip_if_not_installed("pscl")
  skip_if_not_installed("glmmTMB")
  d <- zeug_b97(); d$id <- rep(seq_len(20L), each = 10L)
  im <- suppressWarnings(countimp(d,
          formulas = list(y ~ x1 + x2 + (1 | id)),
          family = zi_poisson(), m = 1L, maxit = 1L, printFlag = FALSE))
  ty <- im$predictorMatrix["y", ]
  ## unchanged: 1 = both parts, -2 = grouping. The translation must not reach
  ## the two-level methods, which read 3/4/5/6.
  expect_identical(unname(ty[c("x1", "x2")]), c(1, 1))
  expect_identical(unname(ty["id"]), -2)
})

test_that("B97: without `zero` both parts keep every predictor", {
  skip_if_not_installed("pscl")
  d <- zeug_b97()
  im <- suppressWarnings(countimp(d, formulas = list(y ~ x1 + x2),
          family = zi_poisson(), m = 1L, maxit = 1L, printFlag = FALSE))
  ty <- im$predictorMatrix["y", ]
  ty <- ty[names(ty) != "y"]
  f <- ci(".countimp_zi_formula")(names(ty), ty)
  expect_identical(deparse(f), deparse(Y ~ x1 + x2 | x1 + x2))
})
