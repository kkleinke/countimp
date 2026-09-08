## B98: countimp() must not need a random stream that does not exist yet.
##
## The mids object records .Random.seed. In a fresh session that object is
## absent until something draws a random number -- and where there is nothing
## to impute, nothing draws. countimp() then died while assembling its result,
## after the work was done. Inherited from mice, which still fails the same way
## in 3.19.
##
## Run in a separate R process: the defect only shows before the stream exists,
## and testthat itself has long since created it.

test_that("B98: a complete data set does not need an existing .Random.seed", {
  skip_on_cran()                      # starts a second R process
  bib <- .libPaths()
  skript <- tempfile(fileext = ".R")
  writeLines(c(
    sprintf(".libPaths(%s)", paste0("c(", paste(sprintf('"%s"', bib), collapse = ", "), ")")),
    "suppressMessages(library(countimp))",
    "if (exists('.Random.seed', envir = globalenv())) stop('stream already there')",
    "d <- data.frame(y = c(1L, 2L, 3L, 4L, 5L), x = c(1, 2, 3, 4, 5))",
    "invisible(countimp(d, m = 1L, maxit = 1L, printFlag = FALSE))",
    "cat('OK\\n')"), skript)
  aus <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
           c("--vanilla", shQuote(skript)), stdout = TRUE, stderr = TRUE))
  expect_true(any(grepl("^OK$", aus)),
              info = paste(utils::tail(aus, 4), collapse = " | "))
})

test_that("B98: the seed field is filled and usable", {
  d <- data.frame(y = c(1L, NA, 3L, NA, 5L, 2L), x = c(1, 2, 3, 4, 5, 6))
  im <- suppressWarnings(countimp(d, m = 1L, maxit = 1L, printFlag = FALSE))
  expect_false(is.null(im$lastSeedValue))
  expect_type(im$lastSeedValue, "integer")
})
