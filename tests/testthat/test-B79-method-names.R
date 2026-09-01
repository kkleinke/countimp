## B79/B80: a named `method` vector
##
## Both bugs were inherited from the mice 2.46 engine this code was derived
## from, and both were silent -- the imputations came out, under the wrong
## model or not at all, and the returned object looked correct.
##
##  B79  `method` was matched by POSITION even when named. With data in the
##       order (a, b), method = c(b = "poisson", a = "norm") fitted "poisson"
##       to `a`. The returned $method carried the data's names, so the object
##       reported the assignment the user asked for, not the one that ran.
##  B80  a named vector that mentioned only some variables left the others
##       unimputed, because check.method() filled in defaults only when EVERY
##       entry was "". The named variable then came out incomplete as well,
##       since its own predictors still had missing values.

imp_q <- function(...) suppressWarnings(countimp(..., printFlag = FALSE))

dat3 <- function(seed = 1, n = 300) {
  set.seed(seed)
  d <- data.frame(a = rnorm(n), b = stats::rpois(n, 3), cc = rnorm(n))
  d$a[sample(n, 50)] <- NA
  d$b[sample(n, 60)] <- NA
  d$cc[sample(n, 40)] <- NA
  d
}

## --- B79: matching ----------------------------------------------------------

test_that("B79: a named method vector is matched by name", {
  d <- dat3()
  ## names in an order that does not match the data
  r <- imp_q(d, method = c(cc = "norm", b = "poisson", a = "norm"), m = 1)
  expect_equal(r$method[["b"]], "poisson")
  expect_equal(r$method[["a"]], "norm")
  expect_equal(r$method[["cc"]], "norm")
})

test_that("B79: the model that ran is the model that was asked for", {
  ## The decisive check: not what $method reports, but what the imputations
  ## look like. Under the positional bug "poisson" was fitted to a normal
  ## variable, which is detectable in the draws (integers vs continuous).
  d <- data.frame(a = NA_real_, b = NA_real_)
  set.seed(2); n <- 300
  d <- data.frame(a = rnorm(n), b = stats::rpois(n, 3))
  d$a[sample(n, 50)] <- NA
  d$b[sample(n, 60)] <- NA
  r <- imp_q(d, method = c(b = "poisson", a = "norm"), m = 1)
  zb <- unlist(r$imp$b); za <- unlist(r$imp$a)
  expect_true(all(zb == round(zb)))    # poisson went to b
  expect_false(all(za == round(za)))   # norm went to a
})

test_that("B79: an unnamed method vector keeps positional semantics", {
  d <- dat3()
  r <- imp_q(d, method = c("norm", "poisson", "norm"), m = 1)
  expect_equal(unname(r$method), c("norm", "poisson", "norm"))
  ## and a single unnamed method is still recycled over all columns
  r2 <- imp_q(d, method = "norm", m = 1)
  expect_true(all(r2$method == "norm"))
})

test_that("B79: an ambiguous or wrong method specification is rejected", {
  d <- dat3()
  expect_error(imp_q(d, method = c(zz = "poisson"), m = 1),
               "not in the data")
  expect_error(imp_q(d, method = c(b = "poisson", "norm"), m = 1),
               "partly named")
  expect_error(imp_q(d, method = c(b = "poisson", b = "nb"), m = 1),
               "more than once")
})

## --- B80: unmentioned variables --------------------------------------------

test_that("B80: variables absent from a named method vector get the default", {
  d <- dat3()
  r <- imp_q(d, method = c(b = "poisson"), m = 1)
  expect_equal(r$method[["b"]], "poisson")
  expect_true(nzchar(r$method[["a"]]))
  expect_true(nzchar(r$method[["cc"]]))
  ## and the data really are complete afterwards -- including the named
  ## variable, which was left with holes while its predictors stayed missing
  cc <- countimp_complete(r, 1)
  expect_equal(sum(is.na(cc)), 0L)
})

test_that("B80: an explicit empty string still means 'do not impute'", {
  d <- dat3()
  r <- imp_q(d, method = c(b = "poisson", a = "", cc = ""), m = 1)
  expect_equal(r$method[["a"]], "")
  expect_equal(r$method[["cc"]], "")
  cc <- countimp_complete(r, 1)
  expect_equal(sum(is.na(cc$a)), 50L)
  expect_equal(sum(is.na(cc$cc)), 40L)
})

test_that("B80: a complete column gets no method", {
  d <- dat3()
  d$full <- rnorm(nrow(d))
  r <- imp_q(d, method = c(b = "poisson"), m = 1)
  expect_equal(r$method[["full"]], "")
})

test_that("B80: the default respects the variable's type", {
  set.seed(3); n <- 300
  d <- data.frame(num = rnorm(n),
                  bin = factor(sample(c("x", "y"), n, TRUE)),
                  cnt = stats::rpois(n, 2))
  d$num[sample(n, 40)] <- NA
  d$bin[sample(n, 40)] <- NA
  d$cnt[sample(n, 40)] <- NA
  r <- imp_q(d, method = c(cnt = "poisson"), m = 1)
  ## defaultMethod is c(pmm, logreg, polyreg, polr): a two-level factor must
  ## not be handed the numeric default
  expect_equal(r$method[["num"]], "pmm")
  expect_equal(r$method[["bin"]], "logreg")
  expect_equal(r$method[["cnt"]], "poisson")
})

test_that("B79/B80: the formula interface is unaffected", {
  ## The named-method branch must not intercept the formula interface, which
  ## builds `method` itself.
  d <- dat3()
  r <- imp_q(d, formulas = list(b ~ a + cc), family = poisson_count(), m = 1)
  expect_equal(r$method[["b"]], "poisson")
  expect_equal(sum(is.na(countimp_complete(r, 1))), 0L)
})

## --- B81: a formula call that does not compile must not fall through --------

test_that("B81: a formula interface call with a bad family is refused", {
  ## Found while writing the test above: the dispatch condition required a
  ## countimp_family_spec, and anything else fell through the block in silence.
  ## `formulas` and `family` were then dropped as unused arguments and every
  ## variable was imputed by the default method -- countimp() reported success
  ## and returned pmm imputations for a call that asked for poisson.
  d <- dat3()
  expect_error(imp_q(d, formulas = list(b ~ a + cc), family = "poisson", m = 1),
               "character string")
  ## the message names the classic-interface equivalent
  expect_error(imp_q(d, formulas = list(b ~ a + cc), family = "poisson", m = 1),
               "method = c\\(<variable> = \\\"poisson\\\"\\)")
  ## B96 replaced the message here. It used to say "not an object of class
  ## family" -- for `family = poisson()`, which IS an object of class family,
  ## a contradiction. It now names the origin and the countimp equivalent.
  expect_error(imp_q(d, formulas = list(b ~ a + cc), family = stats::poisson(),
                     m = 1), "family object of stats/glm", fixed = TRUE)
  expect_error(imp_q(d, formulas = list(b ~ a + cc), family = stats::poisson(),
                     m = 1), "poisson_count()", fixed = TRUE)
  expect_error(imp_q(d, formulas = list(b ~ a + cc), m = 1),
               "needs a matching `family`")
  expect_error(imp_q(d, family = poisson_count(), m = 1),
               "without `formulas`")
})

test_that("B81: a list of families is still accepted", {
  d <- dat3()
  r <- imp_q(d, formulas = list(b = b ~ a + cc), family = list(b = poisson_count()),
             m = 1)
  expect_equal(r$method[["b"]], "poisson")
})
