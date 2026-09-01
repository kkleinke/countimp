## B83: the categorical base methods must return the same TYPE they were given.
##
## Before this was fixed, polr/polyreg/logreg always returned a factor. With a
## NUMERIC column the imputation loop inserts it via as.numeric(), which yields
## the LEVEL INDEX rather than the value: y = 3 + Poisson(3), observed 3..13,
## came back as 1..11 -- every value shifted by the distance of the lower bound
## to 1, silently. With y starting at 0 the shift was 1 and easy to miss.
##
## The test asserts the property (imputed values lie inside the observed set),
## not the return class of the internal helper -- a later refactor may move
## where the conversion happens.

test_that("polr and polyreg keep numeric counts on their own scale", {
  set.seed(9)
  ## Offsets chosen so a level-index bug cannot pass by coincidence: with an
  ## offset of 0 indices and values differ by only 1.
  for (versatz in c(0L, 3L, 10L)) {
    n <- 300L
    x1 <- stats::rnorm(n); x2 <- stats::rnorm(n)
    y  <- versatz + stats::rpois(n, 3)
    d  <- data.frame(y = y, x1 = x1, x2 = x2)
    d$y[sample.int(n, 80L)] <- NA
    obs_ <- sort(unique(stats::na.omit(d$y)))
    for (meth in c("polr", "polyreg")) {
      im <- suppressWarnings(countimp(d, method = c(y = meth, x1 = "", x2 = ""),
                                      m = 3L, printFlag = FALSE))
      z <- unlist(im$imp$y)
      testthat::expect_true(is.numeric(z))
      ## Every imputed value must be one of the observed categories.
      testthat::expect_true(all(z %in% obs_))
    }
  }
})

test_that("logreg keeps a numeric 0/1 variable numeric", {
  set.seed(2)
  n <- 400L
  x <- stats::rnorm(n)
  d <- data.frame(y = stats::rbinom(n, 1L, stats::plogis(x)), x = x)
  d$y[sample.int(n, 100L)] <- NA
  im <- suppressWarnings(countimp(d, method = c(y = "logreg", x = ""),
                                  m = 3L, printFlag = FALSE))
  z <- unlist(im$imp$y)
  testthat::expect_true(is.numeric(z))
  testthat::expect_setequal(sort(unique(z)), c(0, 1))
})

test_that("factors stay factors with their levels intact", {
  set.seed(3)
  n <- 300L
  d <- data.frame(y = factor(sample(c("a", "b", "c"), n, TRUE)),
                  x = stats::rnorm(n))
  d$y[sample.int(n, 80L)] <- NA
  for (meth in c("polyreg", "polr")) {
    im <- suppressWarnings(countimp(d, method = c(y = meth, x = ""),
                                    m = 2L, printFlag = FALSE))
    v <- countimp_complete(im, 1L)$y
    testthat::expect_true(is.factor(v))
    testthat::expect_setequal(levels(v), c("a", "b", "c"))
    testthat::expect_false(anyNA(v))
  }
})

test_that("the shortcut for a single observed category keeps the type", {
  ## polyreg returns early when one category covers every observed case. That
  ## path had the same defect and is not reached by the tests above.
  cd <- ci(".countimp_draw_category")
  post <- matrix(c(0.2, 0.8, 0.5, 0.5), nrow = 2L, byrow = TRUE)
  set.seed(1)
  znum <- cd(post, c("7", "8"), c("7", "8"), want_numeric = TRUE)
  testthat::expect_true(is.numeric(znum))
  testthat::expect_true(all(znum %in% c(7, 8)))
  set.seed(1)
  zfac <- cd(post, c("7", "8"), c("7", "8"), want_numeric = FALSE)
  testthat::expect_true(is.factor(zfac))
  ## Same draw, same values -- only the type differs.
  testthat::expect_equal(as.numeric(as.character(zfac)), znum)
})

test_that("non-numeric levels fall back to a factor instead of losing data", {
  cd <- ci(".countimp_draw_category")
  post <- matrix(c(0.3, 0.7), nrow = 1L)
  z <- suppressWarnings(cd(post, c("low", "high"), c("low", "high"),
                           want_numeric = TRUE))
  testthat::expect_true(is.factor(z))
})
