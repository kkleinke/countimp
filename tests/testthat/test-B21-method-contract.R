## B21 -- the univariate-method contract (ci(".countimp_call_method"))
##
## countimp's own FCS engine calls imputation methods that may come from
## countimp, from the user's workspace, or from mice. The calling convention is
## therefore an interface, and these tests pin it down. Regression guarded:
## the sampler inherited from mice 2.46 passed `x` as a data.frame, which broke
## every method that requires a matrix (mice.impute.lasso.norm, via glmnet).

contract <- ci(".countimp_call_method")

test_that("predictors reach the method as a numeric matrix", {
  set.seed(1); n <- 80
  y <- rnorm(n); ry <- rep(c(TRUE, FALSE), each = n / 2)
  xdf <- data.frame(a = rnorm(n), b = rnorm(n))
  seen <- NULL
  f <- function(y, ry, x, wy = NULL, type = NULL, ...) {
    seen <<- list(cls = class(x)[1], nms = colnames(x))
    rep(0, sum(if (is.null(wy)) !ry else wy))
  }
  contract(f, y = y, ry = ry, x = xdf, type = c(a = 1, b = 1))
  expect_identical(seen$cls, "matrix")
  expect_identical(seen$nms, c("a", "b"))
  ## a matrix is passed through untouched
  contract(f, y = y, ry = ry, x = as.matrix(xdf), type = c(a = 1, b = 1))
  expect_identical(seen$cls, "matrix")
})

test_that("the first three arguments are matched by name, not by position", {
  set.seed(2); n <- 40
  y <- rnorm(n); ry <- rep(c(TRUE, FALSE), each = n / 2)
  x <- matrix(rnorm(2 * n), n, 2, dimnames = list(NULL, c("a", "b")))
  ## arguments deliberately declared in a different order
  f <- function(x, ry, y, wy = NULL, ...) {
    stopifnot(is.matrix(x), is.logical(ry), is.numeric(y))
    rep(99, sum(if (is.null(wy)) !ry else wy))
  }
  expect_equal(unique(contract(f, y = y, ry = ry, x = x)), 99)
})

test_that("`type` is withheld from methods that cannot take it", {
  set.seed(3); n <- 40
  y <- rnorm(n); ry <- rep(c(TRUE, FALSE), each = n / 2)
  x <- matrix(rnorm(2 * n), n, 2, dimnames = list(NULL, c("a", "b")))
  ## no `type`, no `...`: passing type would be an unused-argument error
  f_notype <- function(y, ry, x, wy = NULL) rep(1, sum(if (is.null(wy)) !ry else wy))
  expect_silent(contract(f_notype, y = y, ry = ry, x = x, type = c(a = 1, b = 1)))
  ## declares `type`: it must arrive, aligned with the columns of x
  got <- NULL
  f_type <- function(y, ry, x, wy = NULL, type = NULL) {
    got <<- type; rep(1, sum(if (is.null(wy)) !ry else wy))
  }
  contract(f_type, y = y, ry = ry, x = x, type = c(a = 1, b = -2))
  expect_equal(got, c(a = 1, b = -2))
})

test_that("misaligned `type` is an error, not a silent mis-modelling", {
  set.seed(4); n <- 30
  y <- rnorm(n); ry <- rep(c(TRUE, FALSE), each = n / 2)
  x <- matrix(rnorm(2 * n), n, 2, dimnames = list(NULL, c("a", "b")))
  f <- function(y, ry, x, wy = NULL, type = NULL) rep(1, sum(!ry))
  expect_error(contract(f, y = y, ry = ry, x = x, type = 1),
               "1 entries but `x` has 2 columns")
})

test_that("an unexpanded factor is refused rather than modelled as codes", {
  set.seed(5); n <- 30
  y <- rnorm(n); ry <- rep(c(TRUE, FALSE), each = n / 2)
  xdf <- data.frame(a = rnorm(n), g = factor(sample(letters[1:3], n, TRUE)))
  f <- function(y, ry, x, wy = NULL, type = NULL) rep(1, sum(!ry))
  expect_error(contract(f, y = y, ry = ry, x = xdf), "unexpanded")
})

test_that("stray dots are not forwarded to methods without ...", {
  set.seed(6); n <- 30
  y <- rnorm(n); ry <- rep(c(TRUE, FALSE), each = n / 2)
  x <- matrix(rnorm(2 * n), n, 2, dimnames = list(NULL, c("a", "b")))
  f <- function(y, ry, x, wy = NULL, donors = 5) rep(donors, sum(!ry))
  ## `ridge` is not a formal of f and must be dropped; `donors` must arrive
  expect_equal(unique(contract(f, y = y, ry = ry, x = x, donors = 3, ridge = 1e-5)), 3)
})

test_that("mice methods run through countimp's own engine", {
  ## pscl is in Suggests, not in Imports: the one-level zero-inflation and
  ## hurdle methods need it, the rest of the package does not. Without this
  ## skip the test fails in an environment without Suggests -- and it fails on
  ## the missing dependency, not on the property under test.
  skip_if_not_installed("pscl")
  skip_if_not_installed("mice")
  skip_if_not_installed("glmnet")
  skip_on_cran()
  set.seed(11); n <- 150
  v  <- rnorm(n)
  y  <- rpois(n, exp(1 + 0.4 * v)); y[rbinom(n, 1, 0.3) == 1] <- 0
  b  <- factor(rbinom(n, 1, plogis(0.3 * v)))
  d  <- data.frame(y = y, v = v, b = b)
  d$y[sample.int(n, 30)] <- NA
  d$v[sample.int(n, 30)] <- NA
  pm <- matrix(0, 3, 3, dimnames = list(names(d), names(d)))
  pm["y", ] <- c(0, 1, 1); pm["v", ] <- c(1, 0, 1)
  ## "zinb" forces countimp's own sampler; lasso.norm comes from mice and
  ## requires a matrix -- this is the combination that used to fail.
  imp <- suppressWarnings(countimp(d, method = c("zinb", "lasso.norm", ""),
                                   predictorMatrix = pm, m = 2, maxit = 1,
                                   printFlag = FALSE, seed = 3))
  expect_s3_class(imp, "mids")
  expect_true(all(!is.na(unlist(imp$imp$v))))
  expect_true(all(!is.na(unlist(imp$imp$y))))
})
