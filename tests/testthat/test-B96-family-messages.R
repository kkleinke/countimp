## B96 -- what happens when the wrong family is handed in
##
## Two occasions, both found while closing B12.
##
## 1. `family = poisson()` is the likeliest slip of all: it is the spelling
##    from glm(), it exists, and it returns an object of class "family". The
##    message answering it read "must be a family object ... not an object of
##    class family" -- a contradiction no reader can resolve. (I walked into it
##    myself, with `family = list(CCRIM = poisson())` in a check script.)
##
## 2. The message listed the available constructors and named six of thirteen
##    -- the scale families from B89 were missing. The same hand-maintained
##    list sat in `analyse/k27_b12_stand.R`, where it cost more: it made twelve
##    method names appear "not reachable through formulas" when they had been
##    reachable for weeks, and the project record carried the wrong figure as
##    an open item for two days.
##
## The first block therefore holds the list AGAINST THE NAMESPACE. A
## hand-maintained list is defensible when a test catches up with it.

test_that("B96: the family list knows every exported constructor", {
  required <- list(censored_poisson = list(censor = 10),
                  censored_nb      = list(censor = 10),
                  bounded_poisson  = list(bounds = c(0, 5)),
                  bounded_nb       = list(bounds = c(0, 5)))
  ## Every export returning a countimp_family_spec is a constructor.
  found <- Filter(function(n) {
    f <- try(do.call(n, if (!is.null(required[[n]])) required[[n]] else list()),
             silent = TRUE)
    !inherits(f, "try-error") && inherits(f, "countimp_family_spec")
  }, getNamespaceExports("countimp"))
  expect_setequal(found, ci(".countimp_familien")())
})

test_that("B96: a stats family object is named as such", {
  m <- ci(".countimp_familie_falsch")(poisson())
  expect_match(m, "stats/glm", fixed = TRUE)
  expect_match(m, "poisson_count()", fixed = TRUE)
  ## The contradiction of the old message must not return.
  expect_false(grepl("must be a family object", m, fixed = TRUE))
})

test_that("B96: quasipoisson is refused with a reason, not merely refused", {
  ## The refusal is a design decision: quasi-likelihood has no posterior for
  ## `draw = "bayes"` to draw from. Whoever needs the method gets it through
  ## the classic interface -- and the message has to name that route, or the
  ## refusal is a dead end.
  m <- ci(".countimp_familie_falsch")(quasipoisson())
  expect_match(m, "no quasipoisson family on purpose", fixed = TRUE)
  expect_match(m, "posterior", fixed = TRUE)
  expect_match(m, 'method = c(<variable> = "quasipoisson")', fixed = TRUE)
})

test_that("B96: inside a list it says WHICH element is wrong", {
  m <- ci(".countimp_familie_falsch")(list(a = zi_poisson(), b = poisson()))
  expect_match(m, "`family$b`", fixed = TRUE)
  expect_match(m, "poisson_count()", fixed = TRUE)
  ## and does not blame the good element
  expect_false(grepl("`family$a`", m, fixed = TRUE))
})

test_that("B96: a character string points to `method`", {
  m <- ci(".countimp_familie_falsch")("zip")
  expect_match(m, 'method = c(<variable> = "zip")', fixed = TRUE)
})

test_that("B96: the messages also come out of countimp() itself", {
  d <- data.frame(y = c(0, 1, 2, NA, 3, 1), x = c(1, 2, 3, 4, 5, 6))
  expect_error(countimp(d, formulas = list(y ~ x), family = quasipoisson()),
               "no quasipoisson family on purpose", fixed = TRUE)
  expect_error(countimp(d, formulas = list(y ~ x), family = list(y = poisson())),
               "poisson_count()", fixed = TRUE)
})

test_that("B96: no function is defined twice in the package", {
  ## The general lesson from that finding. While adding the messages I defined
  ## `.countimp_short_message` -- a name diagnostics.R already uses (it shortens message
  ## texts there). R kept the existing definition, mine had no effect, and the
  ## message came out truncated. With the opposite load order my version would
  ## have won and broken message shortening throughout the package -- a defect
  ## no test in this suite would have looked for.
  sources <- Sys.glob(file.path("..", "..", "R", "*.R"))
  skip_if(!length(sources), "source files not reachable (installed package)")
  names_seen <- unlist(lapply(sources, function(f) {
    z <- readLines(f, warn = FALSE)
    t <- grep("^[a-zA-Z._][a-zA-Z0-9._]* *<- *function", z, value = TRUE)
    sub(" *<-.*$", "", t)
  }))
  twice <- names(which(table(names_seen) > 1L))
  expect_identical(twice, character(0),
                   info = paste("defined more than once:",
                                paste(twice, collapse = ", ")))
})
