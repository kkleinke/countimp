## B58 -- the documented return class must be the measured one.
##
## Until this test, @return of countimp() promised "class mids ... when mice is
## available, otherwise the package's own equivalent". B55 had replaced that by
## an unconditional c("countimp_mids", "mids"), so the sentence described a
## branch that no longer existed, on the package's main help page. The same
## sentence sent users to countimp_check() for diagnostics -- a function that no
## imputation method calls (B56) and whose fit/mu arguments cannot be obtained
## from a finished imputation object.
##
## These checks compare the CLAIM in the Rd files against the MEASURED class, so
## a future change to either side has to touch this file.

testthat::test_that("countimp() returns countimp_mids inheriting from mids", {
  d <- .b58_daten()
  imp <- suppressWarnings(countimp(d, method = c("", "nb"), m = 2L, maxit = 2L,
                                   printFlag = FALSE, seed = 1L))
  testthat::expect_identical(class(imp), c("countimp_mids", "mids"))
  testthat::expect_true(inherits(imp, "mids"))
})

testthat::test_that("with() returns countimp_mira inheriting from mira", {
  d <- .b58_daten()
  imp <- suppressWarnings(countimp(d, method = c("", "nb"), m = 2L, maxit = 2L,
                                   printFlag = FALSE, seed = 1L))
  a <- with(imp, stats::lm(y ~ x))
  testthat::expect_identical(class(a), c("countimp_mira", "mira"))
  testthat::expect_true(inherits(a, "mira"))
})

testthat::test_that("the return class does not depend on whether mice is loaded", {
  ## The retired claim was conditional on mice. Measure the assignment itself:
  ## it must be unconditional in the source, not wrapped in a requireNamespace
  ## branch.
  quelle <- .countimp_quelle("R/zimice.R")
  testthat::skip_if(is.null(quelle), "package sources not locatable")
  zeile <- grep('oldClass(midsobj)', quelle, fixed = TRUE, value = TRUE)
  testthat::expect_length(zeile, 1L)
  testthat::expect_match(zeile, 'c\\("countimp_mids", *"mids"\\)')
  ## No mice check anywhere in the function that assigns the class.
  testthat::expect_false(any(grepl('requireNamespace\\("mice"', quelle)))
})

testthat::test_that("the documented return class matches the measured one", {
  rd <- .countimp_quelle("man/countimp.Rd")
  testthat::skip_if(is.null(rd), "Rd file not available")
  wert <- paste(rd, collapse = " ")
  testthat::expect_match(wert, 'countimp_mids')
  ## The retired sentence promised a plain mids object conditional on mice.
  testthat::expect_false(grepl("otherwise the package's own equivalent", wert, fixed = TRUE))
})

testthat::test_that("countimp() does not point users at countimp_check()", {
  ## countimp_check() needs fit and mu, which a user cannot obtain from a
  ## finished mids object; countimp_diagnostics() is the reachable entry point.
  quelle <- .countimp_quelle("R/countimp.R")
  testthat::skip_if(is.null(quelle), "package sources not locatable")
  testthat::expect_false(any(grepl("link{countimp_check}", quelle, fixed = TRUE)))
  testthat::expect_true(any(grepl("link{countimp_diagnostics}", quelle, fixed = TRUE)))
})
