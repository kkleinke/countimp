## B35 -- pscl reads theta BY NAME, and a nameless theta fails silently.
##
## pscl::predict.hurdle() reads object$theta["count"]. A hurdle fit carries
## theta as c(count = .); a zeroinfl fit carries it unnamed. So:
##
##   g <- f; g$theta <- 2.7        # hurdle fit
##   predict(g, type = "zero")     # NaN for every row -- and NO error
##
## No warning, no condition, just unusable numbers that would flow straight
## into the imputations. countimp does not currently replace theta on a pscl
## fit -- it only reads it via ci(".countimp_zi_theta")() -- so this is a guard
## against a future change, and a guard on the reader's own name handling.
##
## Found while writing the k03 measurement script, which did exactly this and
## lost all 400 hurdle replications to it.

testthat::test_that("B35: the theta reader handles both pscl shapes", {
  ## pscl is in Suggests, not in Imports: the one-level zero-inflation and
  ## hurdle methods need it, the rest of the package does not. Without this
  ## skip the test fails in an environment without Suggests -- and it fails on
  ## the missing dependency, not on the property under test.
  skip_if_not_installed("pscl")
  ## zeroinfl: unnamed scalar. hurdle: c(count = .). The reader must return the
  ## count-part theta in both cases, preferring the name over the position.
  mk <- function(theta, dist) {
    f <- list(theta = theta, dist = dist)
    class(f) <- "zeroinfl"
    f
  }
  testthat::expect_equal(ci(".countimp_zi_theta")(mk(2.5, "negbin")), 2.5)
  testthat::expect_equal(ci(".countimp_zi_theta")(mk(c(count = 2.5), "negbin")), 2.5)
  ## a named vector where count is NOT first: the name must win over position
  testthat::expect_equal(
    ci(".countimp_zi_theta")(mk(c(zero = 9.9, count = 2.5), list(count = "negbin"))), 2.5)
  ## Poisson count part -> Inf, whatever theta says
  testthat::expect_equal(ci(".countimp_zi_theta")(mk(c(count = 2.5), "poisson")), Inf)
})

testthat::test_that("B35: an invalid theta stops instead of imputing", {
  mk <- function(theta) {
    f <- list(theta = theta, dist = "negbin"); class(f) <- "zeroinfl"; f
  }
  ## NA is what a name-based read of a nameless theta produces
  testthat::expect_error(ci(".countimp_zi_theta")(mk(c(count = NA_real_))), "did not converge")
  testthat::expect_error(ci(".countimp_zi_theta")(mk(-1)), "not a valid dispersion")
  testthat::expect_error(ci(".countimp_zi_theta")(mk(0)), "not a valid dispersion")
  testthat::expect_error(ci(".countimp_zi_theta")(mk(numeric(0))), "carries no theta")
})

testthat::test_that("B35: the package never assigns to a pscl fit's theta", {
  ## pscl is in Suggests, not in Imports: the one-level zero-inflation and
  ## hurdle methods need it, the rest of the package does not. Without this
  ## skip the test fails in an environment without Suggests -- and it fails on
  ## the missing dependency, not on the property under test.
  skip_if_not_installed("pscl")
  ## Reading theta is safe; replacing it is the trap. If a future change starts
  ## replacing it, this test fails and points at the name requirement above.
  cand <- c("../../R", "../../../R", "R", "paket/R")
  hit <- cand[vapply(cand, function(p)
    dir.exists(p) && length(list.files(p, pattern = "[.]R$")) > 0L, logical(1))]
  testthat::skip_if(length(hit) == 0L, "package sources not locatable")
  files <- list.files(hit[1], pattern = "[.]R$", full.names = TRUE)
  offenders <- character(0)
  for (f in files) {
    src <- sub("#.*$", "", readLines(f, warn = FALSE))
    src <- gsub('"[^"]*"', '""', src)
    ## fit$theta <- ... / fit[["theta"]] <- ...  on a pscl object
    hits <- grep("(fit|f)[[:alnum:]._]*\\$theta[[:space:]]*<-|\\[\\[\"theta\"\\]\\][[:space:]]*<-",
                 src)
    for (i in hits)
      offenders <- c(offenders, sprintf("%s:%d %s", basename(f), i, trimws(src[i])))
  }
  testthat::expect_equal(offenders, character(0),
                         info = paste(offenders, collapse = "\n"))
})

testthat::test_that("B35: pscl still reads theta by name (live check)", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("pscl")
  ## If pscl ever stops reading theta by name this test fails, which is the
  ## signal that the reader's name preference can be simplified.
  set.seed(20260821)
  n <- 300
  x1 <- stats::rnorm(n); x2 <- stats::rnorm(n)
  y <- ifelse(stats::rbinom(n, 1L, stats::plogis(-0.4 + 0.35 * x2)) == 1L, 0L,
              stats::rnbinom(n, size = 1.8, mu = exp(1 + 0.45 * x1)))
  d <- data.frame(Y = y, x1 = x1, x2 = x2)
  f <- suppressWarnings(pscl::hurdle(Y ~ x1 + x2 | x1 + x2, data = d, dist = "negbin"))
  testthat::expect_true("count" %in% names(f$theta))
  ## the reader gets it right on the real object
  testthat::expect_equal(ci(".countimp_zi_theta")(f), as.numeric(f$theta[["count"]]))
  ## and pscl reports an SE for it, so a proper theta draw IS available
  testthat::expect_true(is.numeric(f$SE.logtheta) && is.finite(f$SE.logtheta))
})
