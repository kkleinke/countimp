## B84: the zero-inflated and hurdle methods must refuse zero-free data with a
## message that names the cause -- the mirror image of B73 for ztp/ztnb.
##
## Before this, the check lived inside the fit and was raised from within the
## redraw closure. .countimp_draw_retry() repeated it three times and reported
## "the imputation model ... failed on 3 successive draws (... Extreme predictor
## values, a separated model or too few observed cases are the usual causes ...)"
## -- the wrong diagnosis for a user who simply picked the wrong method. Worse,
## the pscl message that DID name the cause ("invalid dependent variable,
## minimum count is not zero") was destroyed on the way out: the extraction used
## sub("^Error[^:]*:", "") on the try-error text, which stops at the first colon
## and therefore kept the failing CALL and dropped the reason.

test_that("zero-free data is refused by name, not as a convergence failure", {
  ## pscl is in Suggests, not in Imports: the one-level zero-inflation and
  ## hurdle methods need it, the rest of the package does not. Without this
  ## skip the test fails in an environment without Suggests -- and it fails on
  ## the missing dependency, not on the property under test.
  skip_if_not_installed("pscl")
  set.seed(4711)
  n <- 400L
  x1 <- stats::rnorm(n); x2 <- stats::rnorm(n)
  ## Zero-truncated Poisson: no zeros by construction.
  y <- ci(".countimp_rktp")(n, k = 0, mu = exp(1 + 0.4 * x1))
  testthat::expect_true(min(y) >= 1L)
  d <- data.frame(y = y, x1 = x1, x2 = x2)
  d$y[sample.int(n, 100L)] <- NA

  for (meth in c("zip", "zinb", "hp", "hnb")) {
    e <- tryCatch(suppressWarnings(
      countimp(d, method = c(y = meth, x1 = "", x2 = ""),
               m = 1L, maxit = 1L, printFlag = FALSE)),
      error = function(e) e)
    testthat::expect_s3_class(e, "error")
    msg <- conditionMessage(e)
    ## The cause, and the way out.
    testthat::expect_match(msg, "none of the [0-9]+ observed value\\(s\\) is zero")
    testthat::expect_match(msg, "ztp", fixed = TRUE)
    ## NOT the retry wrapper and NOT the misleading predictor advice.
    testthat::expect_false(grepl("successive draws", msg, fixed = TRUE))
    testthat::expect_false(grepl("Extreme predictor values", msg, fixed = TRUE))
  }
})

test_that("data with zeros still runs", {
  ## pscl is in Suggests, not in Imports: the one-level zero-inflation and
  ## hurdle methods need it, the rest of the package does not. Without this
  ## skip the test fails in an environment without Suggests -- and it fails on
  ## the missing dependency, not on the property under test.
  skip_if_not_installed("pscl")
  set.seed(5)
  n <- 400L
  x1 <- stats::rnorm(n); z1 <- stats::rnorm(n)
  y <- ifelse(stats::runif(n) < stats::plogis(2 * z1), 0L,
              stats::rpois(n, exp(1 + 0.3 * x1)))
  testthat::expect_gt(sum(y == 0L), 0L)
  d <- data.frame(y = y, x1 = x1, z1 = z1)
  d$y[sample.int(n, 80L)] <- NA
  for (meth in c("zip", "hp")) {
    im <- suppressWarnings(countimp(d, method = c(y = meth, x1 = "", z1 = ""),
                                    m = 1L, maxit = 1L, printFlag = FALSE))
    testthat::expect_false(anyNA(unlist(im$imp$y)))
  }
})

test_that("the pscl message survives extraction intact", {
  ## pscl is in Suggests, not in Imports: the one-level zero-inflation and
  ## hurdle methods need it, the rest of the package does not. Without this
  ## skip the test fails in an environment without Suggests -- and it fails on
  ## the missing dependency, not on the property under test.
  skip_if_not_installed("pscl")
  ## Direct probe of the extraction path: a genuine pscl failure must arrive
  ## with its own wording, not with the call fragment.
  zf <- ci(".countimp_zi_fit")
  ## Zeros present, so the pre-check passes and pscl itself is reached. Negative
  ## counts make pscl refuse with a message of its own -- chosen over a
  ## too-small data set, which pscl still fits (with a singular Hessian) rather
  ## than rejecting.
  set.seed(1)
  dat <- data.frame(Y = c(0L, -1L, 2L, 3L, 0L, 4L),
                    X = stats::rnorm(6))
  e <- tryCatch(suppressWarnings(zf(Y ~ X | X, dat, "zeroinfl", "poisson")),
                error = function(e) e)
  testthat::expect_s3_class(e, "error")
  msg <- conditionMessage(e)
  testthat::expect_match(msg, "The original message was")
  ## The reason itself must arrive, not just the wrapper.
  testthat::expect_match(msg, "negative counts")
  ## The old sub() left a leading ":" plus the call; neither may appear.
  testthat::expect_false(grepl("original message was: :", msg, fixed = TRUE))
  testthat::expect_false(grepl("zeroinfl(form", msg, fixed = TRUE))
})

test_that("the check function itself is exact about the boundary", {
  ck <- ci(".countimp_zi_check_zeros")
  testthat::expect_error(ck(c(1L, 2L, 3L), "zeroinfl"), "is zero")
  testthat::expect_true(ck(c(0L, 1L, 2L), "zeroinfl"))
  ## A single zero is enough -- the model is then estimable in principle.
  testthat::expect_true(ck(c(0L, 5L, 5L), "hurdle"))
  ## NA must not count as a zero.
  testthat::expect_error(ck(c(NA, 1L, 2L), "hurdle"), "is zero")
  ## An all-NA vector has nothing to judge; that is the caller's problem.
  testthat::expect_true(ck(c(NA_integer_, NA_integer_), "hurdle"))
})
