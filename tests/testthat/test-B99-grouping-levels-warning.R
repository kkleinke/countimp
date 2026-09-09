## B99: a grouping term with very few levels is almost always a mistyped
## two-part model.
##
## `y ~ x | z` is how pscl writes count | zero. In an R formula `|` is the
## grouping operator, so countimp reads it as a random effect over z and picks
## a two-level method -- with a dichotomous z, over two clusters. Nothing said
## so; glmmTMB reported a non-positive-definite Hessian, which names the
## arithmetic and not the mistake.
##
## The warning explains, it does not change the choice: the method stays what
## the formula asked for.

test_that("B99: a two-level grouping term warns and points at `zero`", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("pscl")
  utils::data("crim4w", package = "countimp")
  d <- crim4w[, c("BCRIM", "FEMALE", "RE", "GY")]
  expect_equal(length(unique(stats::na.omit(d$GY))), 2L)   # the premise
  w <- character(0)
  im <- withCallingHandlers(
    suppressMessages(countimp(d, formulas = list(BCRIM ~ FEMALE + RE | GY),
      family = hurdle_poisson(), m = 1L, maxit = 1L, printFlag = FALSE)),
    warning = function(x) { w <<- c(w, conditionMessage(x))
                            invokeRestart("muffleWarning") })
  eigene <- grep("grouping term", w, value = TRUE)
  expect_length(eigene, 1L)
  ## it must name the way out, not just complain about the level count
  expect_match(eigene, "`zero`", fixed = TRUE)
  expect_match(eigene, "grouping operator", fixed = TRUE)
  expect_match(eigene, "GY", fixed = TRUE)
  ## and it explains rather than overrides
  expect_identical(im$method[["BCRIM"]], "2l.hp")
})

test_that("B99: a real multilevel call stays silent", {
  skip_if_not_installed("glmmTMB")
  utils::data("crim4l", package = "countimp")
  d <- crim4l[, c("DELINQ", "FEMALE", "TIME", "ID")]
  expect_gt(length(unique(stats::na.omit(d$ID))), 5L)      # the premise
  w <- character(0)
  withCallingHandlers(
    suppressMessages(countimp(d, formulas = list(DELINQ ~ FEMALE + TIME + (1 | ID)),
      family = hurdle_nb(), m = 1L, maxit = 1L, printFlag = FALSE)),
    warning = function(x) { w <<- c(w, conditionMessage(x))
                            invokeRestart("muffleWarning") })
  expect_length(grep("grouping term", w, value = TRUE), 0L)
})
