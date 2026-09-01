## test-B53-spec.R -- the formula interface (countimp_spec + the six families).
##
## Why this file exists: k10_exportabdeckung.R wrapped every exported function
## in a counter and ran the whole suite. Result -- countimp_spec() and all six
## family constructors were called ZERO times by 806 checks. The headline
## feature of this version was untested, and that is exactly where B53 sat:
## a one-sided `zero = ~ ...` formula (the notation used by the documentation
## and every example) was silently dropped, the zero part inherited the count
## part's random structure, and the method name came out without its .noint
## suffix. Wrong model, no error, no warning.
##
## The checks are cheap on purpose -- countimp_spec() only compiles arguments,
## it fits nothing, so this whole file runs in well under a second and can sit
## in the CRAN-time budget.

zeug <- function(n = 200L, seed = 1L) {
  set.seed(seed)
  d <- data.frame(y  = rpois(n, 2), y2 = rpois(n, 1),
                  x1 = rnorm(n), x2 = rnorm(n), z1 = rnorm(n),
                  g  = factor(rep(seq_len(n %/% 10L), each = 10L)))
  d$y[sample(n, n %/% 5L)]  <- NA
  d$y2[sample(n, n %/% 5L)] <- NA
  d
}
m1 <- function(s, y = "y") as.character(s$method[[y]])

## countimp_spec(), quietly. zeug() has TWO incomplete variables while most
## checks here give a formula for one of them, so the (correct, wanted) message
## "no formula given for incomplete variable(s) y2" fires on nearly every call
## and put 17 lines of noise into the R CMD check log. Swallowing it here keeps
## the log readable; the one test that is ABOUT the message calls
## countimp_spec() directly so expect_message() still sees it.
cs_still <- function(...) suppressMessages(ci_home()$countimp_spec(...))

## Every method name a family constructor can generate must exist as an export.
## Without this, adding a constructor for a family that has no two-level (or no
## .noint) counterpart compiles silently and fails later with "Imputation method
## 'X' was not found" -- at imputation time, in the user's run, not here.
testthat::test_that("B53: every generatable method name exists as an export", {
  fams <- c("hurdle_nb", "hurdle_poisson", "zi_nb", "zi_poisson",
            "nb", "poisson_count", "compois")
  nm <- ci(".countimp_method_name")
  namen <- character(0)
  for (f in fams) {
    fam <- get(f, envir = ci_home())()
    for (twol in c(FALSE, TRUE))
      for (ic in c(TRUE, FALSE))
        for (iz in c(TRUE, FALSE))
          for (dr in c("bayes", "boot"))
            namen <- c(namen, nm(fam, list(twolevel = twol, rint = ic),
                                 list(rint = iz), draw = dr)$method)
  }
  namen <- unique(namen)
  expect_gt(length(namen), 30L)
  fehlt <- namen[!vapply(paste0("mice.impute.", namen), function(x)
    exists(x, envir = ci_home(), inherits = FALSE), logical(1))]
  expect_identical(fehlt, character(0))
})

testthat::test_that("B53: the family constructors return well-formed specs", {
  fams <- c("hurdle_nb", "hurdle_poisson", "zi_nb", "zi_poisson",
            "nb", "poisson_count", "compois")
  for (f in fams) {
    o <- get(f, envir = ci_home())()
    testthat::expect_s3_class(o, "countimp_family_spec")
    testthat::expect_true(all(c("stem", "twopart", "label") %in% names(o)),
                          info = f)
    testthat::expect_true(is.character(o$stem) && length(o$stem) == 1L, info = f)
    testthat::expect_true(is.logical(o$twopart) && length(o$twopart) == 1L,
                          info = f)
  }
  ## twopart must match the family, not be constant
  tp <- vapply(fams, function(f) get(f, envir = ci_home())()$twopart, logical(1))
  testthat::expect_true(all(tp[c("hurdle_nb", "hurdle_poisson",
                                 "zi_nb", "zi_poisson")]))
  testthat::expect_false(any(tp[c("nb", "poisson_count", "compois")]))
})

testthat::test_that("B53: single-level formulas map onto the bare stems", {
  ## pscl is in Suggests, not in Imports: the one-level zero-inflation and
  ## hurdle methods need it, the rest of the package does not. Without this
  ## skip the test fails in an environment without Suggests -- and it fails on
  ## the missing dependency, not on the property under test.
  skip_if_not_installed("pscl")
  d <- zeug()
  cs <- cs_still
  testthat::expect_identical(m1(cs(d, y ~ x1 + x2, family = nb())), "nb")
  testthat::expect_identical(m1(cs(d, y ~ x1, family = poisson_count())),
                             "poisson")
  testthat::expect_identical(m1(cs(d, y ~ x1, zero = ~ z1,
                                   family = zi_poisson())), "zip")
  testthat::expect_identical(m1(cs(d, y ~ x1, zero = ~ z1,
                                   family = hurdle_nb())), "hnb")
  ## draw = "boot" appends the suffix, and only that
  testthat::expect_identical(m1(cs(d, y ~ x1, family = nb(), draw = "boot")),
                             "nb.boot")
})

testthat::test_that("B53: a random-slope-only count part yields .noint", {
  d <- zeug()
  cs <- cs_still
  testthat::expect_identical(m1(cs(d, y ~ x1 + (1 | g), family = nb())),
                             "2l.nb")
  testthat::expect_identical(m1(cs(d, y ~ x1 + (1 + x1 | g), family = nb())),
                             "2l.nb")
  testthat::expect_identical(m1(cs(d, y ~ x1 + (0 + x1 | g), family = nb())),
                             "2l.nb.noint")
})

## --- B53 proper -------------------------------------------------------------
testthat::test_that("B53: a one-sided `zero` formula reaches the zero part", {
  d <- zeug()
  cs <- cs_still
  ## count keeps its intercept, zero drops it -> .noint.zero
  testthat::expect_identical(
    m1(cs(d, y ~ x1 + (1 | g), zero = ~ z1 + (0 + z1 | g),
          family = hurdle_nb())),
    "2l.hnb.noint.zero")
  ## the mirror image -> .noint.count (this is the one that used to come out
  ## as .noint.both, because the count flag was applied to both parts)
  testthat::expect_identical(
    m1(cs(d, y ~ x1 + (0 + x1 | g), zero = ~ z1 + (1 | g),
          family = hurdle_nb())),
    "2l.hnb.noint.count")
  ## both drop it
  testthat::expect_identical(
    m1(cs(d, y ~ x1 + (0 + x1 | g), zero = ~ z1 + (0 + z1 | g),
          family = hurdle_nb())),
    "2l.hnb.noint.both")
  ## neither
  testthat::expect_identical(
    m1(cs(d, y ~ x1 + (1 | g), zero = ~ z1 + (1 | g), family = hurdle_nb())),
    "2l.hnb")
})

testthat::test_that("B53: all three zero notations agree", {
  d <- zeug()
  cs <- cs_still
  a <- cs(d, y ~ x1 + (1 | g), zero = ~ z1 + (0 + z1 | g),
          family = zi_nb())                                   # one-sided
  b <- cs(d, y ~ x1 + (1 | g), zero = y ~ z1 + (0 + z1 | g),
          family = zi_nb())                                   # two-sided
  cc <- cs(d, y ~ x1 + (1 | g), zero = list(y = ~ z1 + (0 + z1 | g)),
           family = zi_nb())                                  # named list
  testthat::expect_identical(m1(a), "2l.zinb.noint.zero")
  testthat::expect_identical(m1(a), m1(b))
  testthat::expect_identical(m1(a), m1(cc))
  ## and the compiled arguments, not just the name
  testthat::expect_identical(a$type[["y"]], b$type[["y"]])
  testthat::expect_identical(a$type[["y"]], cc$type[["y"]])
  testthat::expect_identical(a$predictorMatrix, b$predictorMatrix)
})

testthat::test_that("B53: the zero part actually enters the type vector", {
  d <- zeug()
  cs <- cs_still
  ## z1 appears ONLY in the zero part; if the zero formula were dropped, z1
  ## would carry type 0 and the fitted model would not see it at all.
  s <- cs(d, y ~ x1 + (1 | g), zero = ~ z1 + (1 | g), family = hurdle_nb())
  ty <- s$type[["y"]]
  testthat::expect_true("z1" %in% names(ty))
  testthat::expect_false(ty[["z1"]] == 0L,
    info = paste("z1 type =", ty[["z1"]]))
  ## and the predictor matrix must list it
  testthat::expect_equal(s$predictorMatrix["y", "z1"], 1L,
                         ignore_attr = TRUE)
  ## the group variable is flagged -2, whichever part named it
  testthat::expect_equal(s$predictorMatrix["y", "g"], -2L,
                         ignore_attr = TRUE)
})

testthat::test_that("B53: an unattributable one-sided zero formula is refused", {
  d <- zeug()
  cs <- cs_still
  ## two targets, one nameless zero formula: guessing would be wrong, so stop
  testthat::expect_error(
    cs(d, list(y ~ x1 + (1 | g), y2 ~ x2 + (1 | g)),
       zero = ~ z1 + (1 | g), family = hurdle_nb()),
    "Cannot tell which target")
  ## with a name it goes through, and only to that target
  s <- cs(d, list(y ~ x1 + (1 | g), y2 ~ x2 + (1 | g)),
          zero = list(y2 = ~ z1 + (0 + z1 | g)), family = hurdle_nb())
  testthat::expect_identical(m1(s, "y2"), "2l.hnb.noint.zero")
  testthat::expect_identical(m1(s, "y"),  "2l.hnb")
  ## a zero formula for a target with no count formula is an error
  testthat::expect_error(
    cs(d, y ~ x1, zero = list(y2 = ~ z1), family = hurdle_nb()),
    "no count-part formula")
  ## two zero formulas for one target likewise
  testthat::expect_error(
    cs(d, y ~ x1, zero = list(y = ~ z1, y = ~ x2), family = hurdle_nb()),
    "More than one")
})

testthat::test_that("B53: a zero formula for a one-part family warns", {
  d <- zeug()
  cs <- cs_still
  testthat::expect_warning(
    s <- cs(d, y ~ x1 + (1 | g), zero = ~ z1, family = nb()),
    "no zero part")
  testthat::expect_identical(m1(s), "2l.nb")
})

testthat::test_that("B53: incomplete variables without a formula are reported", {
  d <- zeug()
  cs <- ci_home()$countimp_spec
  testthat::expect_message(
    s <- cs(d, y ~ x1 + (1 | g), family = nb()),
    "no formula given for incomplete variable")
  ## y2 is incomplete and unmodelled -> a default method, not an empty entry
  testthat::expect_true(nzchar(m1(s, "y2")))
  testthat::expect_identical(m1(s, "y2"), "pmm")
})


## A grouping variable may be a factor on the formula route -----------------
##
## Reported from the simulation side on 30 August 2026. The formula syntax
## `(1 | grp)` is lme4's and glmmTMB's, where a factor is the normal case, so a
## user who wrote their analysis model there and reuses the formula for the
## imputation arrives with a factor -- and got
##
##   Class variable (column 5) cannot be factor. Convert to numeric by as.integer()
##
## a message that names a column number instead of a variable, and points at a
## predictorMatrix coding the formula route exists to hide.
##
## The check itself stays: on the classic route the -2 is the user's own, and a
## factor grouping can be expanded into indicator columns further down -- 50
## groups would silently become 49 predictors. What changed is that the formula
## route converts before the check runs, because there the variable is marked
## as a grouping SYNTACTICALLY, so the conversion is unambiguous.

testthat::test_that("B53: a factor grouping is accepted on the formula route", {
  testthat::skip_if_not_installed("glmmTMB")
  set.seed(5)
  n <- 400L
  d <- data.frame(y = rpois(n, 3), x = rnorm(n), z = rnorm(n),
                  grp = rep(1:20, each = 20))
  d$y[sample.int(n, 80L)] <- NA
  run <- function(dd) {
    set.seed(11)
    suppressWarnings(countimp(dd, formulas = list(y ~ x + z + (1 | grp)),
      family = nb(), m = 1L, maxit = 1L, printFlag = FALSE))
  }
  df <- d; df$grp <- factor(df$grp)
  ## not merely "it runs": the factor must give the SAME imputations as the
  ## integer coding, or the conversion changed the model
  testthat::expect_identical(run(df)$imp$y, run(d)$imp$y)
})

testthat::test_that("B53: the conversion does not depend on the labels", {
  testthat::skip_if_not_installed("glmmTMB")
  ## An identifier carries neither order nor metric, so non-consecutive codes
  ## and character labels must land on the same imputations as 1..20.
  set.seed(5)
  n <- 400L
  d <- data.frame(y = rpois(n, 3), x = rnorm(n), z = rnorm(n),
                  grp = rep(1:20, each = 20))
  d$y[sample.int(n, 80L)] <- NA
  run <- function(dd) {
    set.seed(11)
    suppressWarnings(countimp(dd, formulas = list(y ~ x + z + (1 | grp)),
      family = nb(), m = 1L, maxit = 1L, printFlag = FALSE))$imp$y
  }
  ref <- run(d)
  dl <- d; dl$grp <- c(3, 7, 12, 18, 25, 31, 40, 44, 50, 57,
                       63, 70, 77, 81, 88, 92, 95, 101, 107, 110)[d$grp]
  ds <- d; ds$grp <- letters[d$grp]
  testthat::expect_identical(run(dl), ref)
  testthat::expect_identical(run(ds), ref)
})

testthat::test_that("B53: the classic route still refuses a factor grouping", {
  ## The guard against over-correcting. Here the -2 is the user's own, so the
  ## check must stay -- and now name the variable and the way out.
  set.seed(5)
  n <- 200L
  d <- data.frame(y = rpois(n, 3), x = rnorm(n), z = rnorm(n),
                  grp = factor(rep(1:10, each = 20)))
  d$y[sample.int(n, 40L)] <- NA
  p <- matrix(0L, 4L, 4L, dimnames = list(names(d), names(d)))
  p["y", ] <- c(0L, 2L, 1L, -2L)
  err <- tryCatch(suppressWarnings(countimp(d, method = c("2l.nb2", "", "", ""),
           predictorMatrix = p, m = 1L, maxit = 1L, printFlag = FALSE)),
         error = conditionMessage)
  testthat::expect_type(err, "character")
  ## names the variable, not just a column number, and gives a runnable fix
  testthat::expect_match(err, "`grp`", fixed = TRUE)
  testthat::expect_match(err, "as.integer(data$grp)", fixed = TRUE)
  testthat::expect_match(err, "10 levels", fixed = TRUE)
})

testthat::test_that("B53: three levels, both groupings as factors", {
  testthat::skip_if_not_installed("glmmTMB")
  set.seed(5)
  n <- 400L
  d <- data.frame(y = rpois(n, 3), x = rnorm(n), z = rnorm(n),
                  classroom = rep(1:20, each = 20), school = rep(1:5, each = 80))
  d$y[sample.int(n, 80L)] <- NA
  run <- function(dd) {
    set.seed(11)
    suppressWarnings(countimp(dd,
      formulas = list(y ~ x + z + (1 | school) + (1 | classroom)),
      family = nb(), m = 1L, maxit = 1L, printFlag = FALSE))$imp$y
  }
  df <- d; df$classroom <- factor(df$classroom); df$school <- factor(df$school)
  testthat::expect_identical(run(df), run(d))
})
