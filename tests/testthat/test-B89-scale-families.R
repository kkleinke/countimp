## B89: the scale families in the formula interface
##
## `zerotrunc_poisson()`, `bounded_poisson(bounds)`, `censored_poisson(censor)`
## and their NB counterparts put the scale INTO the family object, so that it
## stands next to the variable it describes instead of in a parallel list.
##
## What these tests defend, in the order in which each was easy to get wrong:
##
##  1. THE LEVEL TRAP. .countimp_method_name() composes "2l." with the stem, so
##     a one-level-only family in a two-level formula produced "2l.cp" -- a
##     name that does not exist. Measured before the fix: the translation
##     succeeded and the run failed later with "Imputation method '2l.cp' was
##     not found", i.e. in the user's run rather than at the point of the
##     mistake. test-B53 checks the same property from the other side (every
##     generatable name exists).
##  2. The scale actually reaching the method. It travels family -> spec$args
##     -> countimp() -> per-variable argument, four steps, and a break anywhere
##     in the chain does not fail: the method falls back to its own default or
##     errors about a missing argument, both far from the cause.
##  3. Double specification. The same scale can now come from two places. It is
##     an error, deliberately, and not a precedence rule -- the two can differ,
##     and imputations under a bound that was not meant still look like
##     plausible counts.
##  4. The mandatory argument. bounded_poisson() without a bound would leave
##     the scale in two places again, which is the coupling this interface
##     exists to remove.

imp_s <- function(...) suppressWarnings(suppressMessages(countimp(..., printFlag = FALSE)))

.b89_daten <- function(n = 300L, seed = 89L) {
  set.seed(seed)
  x1 <- stats::rnorm(n)
  d <- data.frame(
    y = pmin(stats::rpois(n, exp(1.4 + 0.5 * x1)), 8L),   # top-coded at 8
    z = pmin(stats::rpois(n, exp(1.0 - 0.3 * x1)), 5L),   # scale 0..5
    w = 1L + stats::rpois(n, exp(0.7 + 0.3 * x1)),        # cannot be zero
    x1 = x1)
  for (v in c("y", "z", "w")) d[[v]][sample.int(n, 50L)] <- NA
  d
}


## --- 1: the family objects ---------------------------------------------------

test_that("B89: the six scale families are well-formed and one-level", {
  fams <- list(zerotrunc_poisson(), zerotrunc_nb(),
               bounded_poisson(c(0, 8)), bounded_nb(c(0, 8)),
               censored_poisson(10), censored_nb(10))
  stems <- c("ztp", "ztnb", "bp", "bnb", "cp", "cnb")
  for (i in seq_along(fams)) {
    f <- fams[[i]]
    expect_s3_class(f, "countimp_family_spec")
    expect_identical(f$stem, stems[i])
    expect_false(f$twopart)
    ## none of them has a two-level counterpart
    expect_identical(f$at_levels, "1l", info = stems[i])
  }
  ## the scale is carried in the family, the truncated ones carry nothing
  expect_identical(bounded_poisson(c(0, 8))$args, list(bounds = c(0, 8)))
  expect_identical(censored_nb(10)$args, list(censor = 10))
  expect_length(zerotrunc_poisson()$args, 0L)
})

test_that("B89: the scale argument is mandatory", {
  expect_error(bounded_poisson(), "needs its `bounds` argument")
  expect_error(bounded_nb(), "needs its `bounds` argument")
  expect_error(censored_poisson(), "needs its `censor` argument")
  expect_error(censored_nb(), "needs its `censor` argument")
  ## and the message points at the way to supply it separately
  expect_error(bounded_poisson(), "classic interface")
})


## --- 2: the level trap -------------------------------------------------------

test_that("B89: a one-level family in a two-level formula is refused", {
  d <- .b89_daten()
  d$id <- rep(seq_len(30), each = 10)
  for (fam in list(censored_poisson(10), bounded_poisson(c(0, 8)),
                   zerotrunc_poisson())) {
    expect_error(
      suppressMessages(countimp_spec(d, y ~ x1 + (1 | id), family = fam)),
      "no two-level counterpart")
  }
  ## the message names the method that does not exist, so the user can look it
  ## up rather than guess
  msg <- tryCatch(suppressMessages(
           countimp_spec(d, y ~ x1 + (1 | id), family = censored_poisson(10))),
         error = function(e) conditionMessage(e))
  expect_match(msg, "2l.cp", fixed = TRUE)

  ## a family that HAS both levels is unaffected
  expect_silent(suppressMessages(
    countimp_spec(d, y ~ x1 + (1 | id), family = poisson_count())))
})


## --- 3: the scale reaches the method ----------------------------------------

test_that("B89: the spec carries the scale arguments per variable", {
  d <- .b89_daten()
  sp <- suppressMessages(countimp_spec(d, list(y ~ x1, z ~ x1, w ~ x1),
          family = list(y = censored_poisson(censor = 8),
                        z = bounded_poisson(bounds = c(0, 5)),
                        w = zerotrunc_poisson())))
  expect_identical(unname(sp$method[c("y", "z", "w")]), c("cp", "bp", "ztp"))
  expect_identical(sp$args$censor, list(y = 8))
  expect_identical(sp$args$bounds, list(z = c(0, 5)))
})

test_that("B89: three scale families in one call impute on their own scales", {
  d <- .b89_daten()
  imp <- imp_s(d, formulas = list(y ~ x1, z ~ x1, w ~ x1),
               family = list(y = censored_poisson(censor = 8),
                             z = bounded_poisson(bounds = c(0, 5)),
                             w = zerotrunc_poisson()), m = 2)
  vy <- unlist(imp$imp$y); vz <- unlist(imp$imp$z); vw <- unlist(imp$imp$w)
  expect_true(all(vy >= 0 & vy <= 8))
  expect_true(all(vz >= 0 & vz <= 5))
  expect_true(all(vw >= 1))                 # zero-truncated: never zero
  expect_true(any(vy == 8))                 # the censoring limit is reached
})

test_that("B89: the bootstrap variants come from draw = 'boot'", {
  d <- .b89_daten()
  sp <- suppressMessages(countimp_spec(d, list(y ~ x1),
          family = list(y = censored_poisson(censor = 8)), draw = "boot"))
  expect_identical(unname(sp$method[["y"]]), "cp.boot")
  imp <- imp_s(d[, c("y", "x1")], formulas = list(y ~ x1),
               family = censored_poisson(censor = 8), draw = "boot", m = 1)
  expect_true(all(unlist(imp$imp$y) <= 8))
})


## --- 4: double specification -------------------------------------------------

test_that("B89: giving the same scale twice is an error, not a precedence rule", {
  d <- .b89_daten()
  expect_error(
    imp_s(d[, c("y", "x1")], formulas = list(y ~ x1),
          family = censored_poisson(censor = 8), censor = list(y = 6), m = 1),
    "given twice")
  ## also when the two agree: one rule, no exception
  expect_error(
    imp_s(d[, c("y", "x1")], formulas = list(y ~ x1),
          family = censored_poisson(censor = 8), censor = list(y = 8), m = 1),
    "given twice")
  ## and for bounds as well
  expect_error(
    imp_s(d[, c("z", "x1")], formulas = list(z ~ x1),
          family = bounded_poisson(bounds = c(0, 5)), bounds = c(0, 5), m = 1),
    "given twice")
})

test_that("B89: the classic interface is untouched", {
  ## the scale families change nothing about method + censor/bounds.
  ## `w` is dropped here on purpose: leaving an incomplete variable with
  ## method "" in the data makes it an incomplete PREDICTOR, and the engine
  ## then returns NA for every row where it is missing -- correct behaviour,
  ## but it would hide what this test is about.
  d <- .b89_daten()[, c("y", "z", "x1")]
  imp <- imp_s(d, method = c(y = "cp", z = "bp", x1 = ""), m = 1,
               maxit = 1, censor = list(y = 8), bounds = list(z = c(0, 5)))
  expect_false(anyNA(unlist(imp$imp$y)))
  expect_true(all(unlist(imp$imp$y) <= 8))
  expect_true(all(unlist(imp$imp$z) <= 5))
})
