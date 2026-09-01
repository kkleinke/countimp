## B55 -- countimp's S3 methods must not displace mice's methods for mice's
## own objects.
##
## Until 3.0.0 countimp registered with.mids, print.mira and summary.mira
## directly. With mice attached, loading countimp printed "Registered S3 methods
## overwritten by 'countimp'" and silently changed the behaviour of summary()
## on objects created by mice itself: mice::summary.mira takes a `type` argument
## and returns a tidy data frame, countimp's returned a plain list and ignored
## `type`. Any script that loaded both packages got the countimp behaviour,
## whichever object it was working on.
##
## The fix subclasses: countimp() returns c("countimp_mids", "mids") and with()
## returns c("countimp_mira", "mira"), with the methods registered on the
## countimp classes. mice keeps its methods for its own objects, and mice
## functions that test inherits(x, "mids") still accept countimp objects.

testthat::test_that("B55: countimp objects carry the subclass, mice's class is inherited", {
  d <- .b55_daten(seed = 11)
  imp <- suppressWarnings(countimp(d, method = c("", "nb"), m = 2, maxit = 1,
                                   printFlag = FALSE, seed = 4))
  testthat::expect_s3_class(imp, "countimp_mids")
  testthat::expect_s3_class(imp, "mids")
  testthat::expect_identical(class(imp), c("countimp_mids", "mids"))

  fit <- with(imp, stats::lm(y ~ x))
  testthat::expect_s3_class(fit, "countimp_mira")
  testthat::expect_s3_class(fit, "mira")
  testthat::expect_identical(class(fit), c("countimp_mira", "mira"))
  ## the spurious "matrix" class is gone: inherits() and is.matrix() agree
  testthat::expect_false(inherits(fit, "matrix"))
  testthat::expect_false(is.matrix(fit))
})

testthat::test_that("B55: the methods are registered on the countimp classes only", {
  home <- ci_home()
  for (nm in c("with.countimp_mids", "print.countimp_mira", "summary.countimp_mira"))
    testthat::expect_true(exists(nm, envir = home, inherits = FALSE),
                          info = paste(nm, "must exist"))
  ## the mice-class methods must NOT be defined here -- that was the defect
  for (nm in c("with.mids", "print.mira", "summary.mira"))
    testthat::expect_false(exists(nm, envir = home, inherits = FALSE),
                           info = paste(nm, "must not be defined by countimp"))

  ## and they must not be registered for the mice classes either. The namespace
  ## registration table is a character matrix: generic, class, method name.
  tab <- .getNamespaceInfo(asNamespace("countimp"), "S3methods")
  testthat::expect_true(is.matrix(tab))
  reg <- paste(tab[, 1L], tab[, 2L], sep = ".")
  testthat::expect_false(any(c("with.mids", "print.mira", "summary.mira") %in% reg))
  ## every method countimp registers is for a countimp class
  testthat::expect_true(all(grepl("^countimp_", tab[, 2L])))
})

testthat::test_that("B55: countimp's methods dispatch on countimp objects", {
  d <- .b55_daten(seed = 12)
  imp <- suppressWarnings(countimp(d, method = c("", "nb"), m = 2, maxit = 1,
                                   printFlag = FALSE, seed = 6))
  fit <- with(imp, stats::lm(y ~ x))
  ## countimp's summary returns one summary per analysis, as a plain list
  s <- summary(fit)
  testthat::expect_type(s, "list")
  testthat::expect_length(s, 2L)
  testthat::expect_s3_class(s[[1L]], "summary.lm")
  ## print returns its argument invisibly and mentions the number of analyses
  ausgabe <- utils::capture.output(r <- print(fit))
  testthat::expect_identical(r, fit)
  testthat::expect_true(any(grepl("2 imputed data sets", ausgabe, fixed = TRUE)))
  ## pooling works on it
  testthat::expect_type(miinference(fit), "list")
})

testthat::test_that("B55: loading countimp leaves mice's own methods untouched", {
  testthat::skip_if_not_installed("mice")
  ## The decisive property: summary() on an object created BY MICE must be
  ## mice's tidy data frame, not countimp's list.
  d <- .b55_daten(seed = 13)
  imp_mice <- mice::mice(d, m = 2, maxit = 2, printFlag = FALSE, seed = 3)
  testthat::expect_false(inherits(imp_mice, "countimp_mids"))

  fit_mice <- with(imp_mice, stats::lm(y ~ x))
  testthat::expect_identical(class(fit_mice), "mira")
  s <- summary(fit_mice)
  testthat::expect_s3_class(s, "data.frame")   # mice returns a tidy frame
  testthat::expect_true("estimate" %in% names(s))
  testthat::expect_false(is.null(nrow(s)))

  ## and countimp's objects still reach mice's own machinery
  imp_ci <- suppressWarnings(countimp(d, method = c("", "nb"), m = 2, maxit = 1,
                                      printFlag = FALSE, seed = 8))
  testthat::expect_true(mice::is.mids(imp_ci))
  testthat::expect_s3_class(mice::complete(imp_ci, 1L), "data.frame")
  testthat::expect_s3_class(mice::pool(with(imp_ci, stats::lm(y ~ x))), "mipo")
})
