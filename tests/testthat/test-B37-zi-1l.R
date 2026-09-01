## B37 -- the eight one-level ZI/hurdle methods are wrappers around one engine.
##
## Guards the consolidation of R/mice.impute.zip.R (369 -> 161 lines) onto
## ci(".countimp_1l_zi")() in R/impute1lzi.R. The properties tested are the ones the
## copied bodies could each get wrong independently.

zi.data <- function(seed = 20260821, n = 300, p.obs = 0.75, k = 2L) {
  set.seed(seed)
  X <- as.data.frame(matrix(stats::rnorm(n * k), n, k))
  colnames(X) <- paste0("x", seq_len(k))
  eta <- 1.0 + 0.45 * X[[1]]
  pz <- stats::plogis(-0.4 + if (k >= 2L) 0.35 * X[[2]] else 0)
  y <- ifelse(stats::rbinom(n, 1L, pz) == 1L, 0L,
              stats::rnbinom(n, size = 1.8, mu = exp(eta)))
  list(y = y, ry = stats::rbinom(n, 1L, p.obs) == 1L, x = X,
       type = stats::setNames(rep(1L, k), colnames(X)))
}

meths <- c("zip", "zip.boot", "zinb", "zinb.boot",
           "hp", "hp.boot", "hnb", "hnb.boot")

testthat::test_that("B37: all eight methods exist and are one-line wrappers", {
  ## pscl is in Suggests, not in Imports: the one-level zero-inflation and
  ## hurdle methods need it, the rest of the package does not. Without this
  ## skip the test fails in an environment without Suggests -- and it fails on
  ## the missing dependency, not on the property under test.
  skip_if_not_installed("pscl")
  for (m in meths) {
    f <- get(paste0("mice.impute.", m))
    testthat::expect_true(is.function(f), info = m)
    testthat::expect_named(formals(f), c("y", "ry", "x", "type", "wy"), info = m)
    ## the body must delegate: one call to the engine, no pscl call of its own
    b <- paste(deparse(body(f)), collapse = " ")
    testthat::expect_true(grepl(".countimp_1l_zi", b, fixed = TRUE), info = m)
    testthat::expect_false(grepl("pscl::", b, fixed = TRUE), info = m)
  }
})

testthat::test_that("B37: the three switches reach the engine correctly", {
  ## Each wrapper must pass its own engine/dist/draw combination. Read them
  ## back from the call rather than trusting the file layout.
  want <- list(
    zip       = c("zeroinfl", "poisson", "bayes"),
    zip.boot  = c("zeroinfl", "poisson", "boot"),
    zinb      = c("zeroinfl", "negbin",  "bayes"),
    zinb.boot = c("zeroinfl", "negbin",  "boot"),
    hp        = c("hurdle",   "poisson", "bayes"),
    hp.boot   = c("hurdle",   "poisson", "boot"),
    hnb       = c("hurdle",   "negbin",  "bayes"),
    hnb.boot  = c("hurdle",   "negbin",  "boot"))
  for (m in names(want)) {
    b <- paste(deparse(body(get(paste0("mice.impute.", m)))), collapse = " ")
    for (v in want[[m]])
      testthat::expect_true(grepl(paste0('"', v, '"'), b, fixed = TRUE),
                            info = paste(m, "->", v))
  }
})

testthat::test_that("B37: every method imputes non-negative integers", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("pscl")
  d <- zi.data()
  for (m in meths) {
    set.seed(4)
    r <- suppressWarnings(get(paste0("mice.impute.", m))(d$y, d$ry, d$x, d$type))
    testthat::expect_equal(length(r), sum(!d$ry), info = m)
    testthat::expect_true(all(is.finite(r)), info = m)
    testthat::expect_true(all(r >= 0), info = m)
    testthat::expect_true(all(abs(r - round(r)) < 1e-8), info = m)
  }
})

testthat::test_that("B37: a single predictor keeps its real name", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("pscl")
  ## The shipped code did X <- x[ry, ] without drop = FALSE, so with one
  ## predictor the frame collapsed to a vector and data.frame() renamed the
  ## column to "X". It only worked because newdata was renamed the same way,
  ## so the mismatch cancelled. The engine keeps the frame instead; this test
  ## fails if the drop = FALSE is ever lost.
  d <- zi.data(k = 1L)
  for (m in meths) {
    set.seed(4)
    r <- suppressWarnings(get(paste0("mice.impute.", m))(d$y, d$ry, d$x, d$type))
    testthat::expect_equal(length(r), sum(!d$ry), info = m)
    testthat::expect_true(all(is.finite(r) & r >= 0), info = m)
  }
  ## and the engine's formula must name the predictor, not "X"
  f <- ci(".countimp_zi_formula")("x1", stats::setNames(1L, "x1"))
  testthat::expect_equal(deparse(f), "Y ~ x1 | x1")
})

testthat::test_that("B37: the two-part formula honours the type codes", {
  nam <- c("a", "b", "c", "d")
  ## 1 = both, 2 = count only, 3 = zero only, 0 = unused
  f <- ci(".countimp_zi_formula")(nam, stats::setNames(c(1L, 2L, 3L, 0L), nam))
  testthat::expect_equal(deparse(f), "Y ~ a + b | a + c")
  ## an empty side becomes intercept-only rather than an error
  f0 <- ci(".countimp_zi_formula")(nam, stats::setNames(c(0L, 0L, 3L, 0L), nam))
  testthat::expect_equal(deparse(f0), "Y ~ 1 | c")
  f1 <- ci(".countimp_zi_formula")(nam, stats::setNames(c(2L, 0L, 0L, 0L), nam))
  testthat::expect_equal(deparse(f1), "Y ~ a | 1")
  ## column ORDER follows the data, not the order of the codes
  f2 <- ci(".countimp_zi_formula")(nam, stats::setNames(c(3L, 1L, 1L, 0L), nam))
  testthat::expect_equal(deparse(f2), "Y ~ b + c | a + b + c")
})

testthat::test_that("B37: the beta draw is taken before the slots are written", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("pscl")
  ## The shipped code re-read length(fit$coefficients$count) AFTER overwriting
  ## that slot, so the split of b.star depended on its own side effect. Harmless
  ## while the lengths matched; this checks the split is right when the two
  ## parts have DIFFERENT lengths.
  d <- zi.data(k = 3L)
  d$type <- stats::setNames(c(1L, 2L, 3L), colnames(d$x))   # count 2, zero 2
  dat <- data.frame(Y = d$y[d$ry], d$x[d$ry, , drop = FALSE])
  form <- ci(".countimp_zi_formula")(colnames(d$x), d$type)
  fit <- suppressWarnings(ci(".countimp_zi_fit")(form, dat, "zeroinfl", "negbin"))
  nc0 <- length(fit$coefficients$count); nz0 <- length(fit$coefficients$zero)
  set.seed(1)
  g <- ci(".countimp_zi_draw_beta")(fit)
  testthat::expect_equal(length(g$coefficients$count), nc0)
  testthat::expect_equal(length(g$coefficients$zero), nz0)
  ## names survive the draw, and the values actually moved
  testthat::expect_named(g$coefficients$count, names(fit$coefficients$count))
  testthat::expect_false(isTRUE(all.equal(unname(g$coefficients$count),
                                          unname(fit$coefficients$count))))
})

testthat::test_that("B37: a rank-deficient model fails with an explanation", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("pscl")
  ## A duplicated predictor makes vcov singular. The engine must say what to do
  ## rather than let chol() error deep inside.
  d <- zi.data(k = 2L)
  d$x$x2 <- d$x$x1                       # exact collinearity
  dat <- data.frame(Y = d$y[d$ry], d$x[d$ry, , drop = FALSE])
  form <- ci(".countimp_zi_formula")(colnames(d$x), d$type)
  ## pscl aborts here with "non-finite value supplied by optim", which tells a
  ## user nothing. Whichever of the two layers catches it, the message must name
  ## countimp and say what to change.
  err <- tryCatch({
    fit <- suppressWarnings(ci(".countimp_zi_fit")(form, dat, "zeroinfl", "poisson"))
    ci(".countimp_zi_draw_beta")(fit)
    NA_character_
  }, error = function(e) conditionMessage(e))
  testthat::expect_false(is.na(err))
  testthat::expect_match(err, "^countimp: ")
  testthat::expect_match(err,
    "positive definite|coefficients but a|failed to fit the imputation model")
})

testthat::test_that("B37: pscl is called in exactly one place", {
  ## pscl is in Suggests, not in Imports: the one-level zero-inflation and
  ## hurdle methods need it, the rest of the package does not. Without this
  ## skip the test fails in an environment without Suggests -- and it fails on
  ## the missing dependency, not on the property under test.
  skip_if_not_installed("pscl")
  ## The point of the consolidation: one call site per fitting function, so a
  ## pscl change is a one-line fix. Fails if a new direct call appears.
  cand <- c("../../R", "../../../R", "R", "paket/R")
  hit <- cand[vapply(cand, function(p)
    dir.exists(p) && length(list.files(p, pattern = "[.]R$")) > 0L, logical(1))]
  testthat::skip_if(length(hit) == 0L, "package sources not locatable")
  files <- list.files(hit[1], pattern = "[.]R$", full.names = TRUE)
  sites <- character(0)
  for (f in files) {
    src <- sub("#.*$", "", readLines(f, warn = FALSE))
    src <- gsub('"[^"]*"', '""', src)
    for (i in grep("pscl::(zeroinfl|hurdle)\\(", src))
      sites <- c(sites, sprintf("%s:%d", basename(f), i))
  }
  ## both calls live in ci(".countimp_zi_fit")()
  testthat::expect_true(all(grepl("^impute1lzi[.]R:", sites)),
                        info = paste(sites, collapse = ", "))
  testthat::expect_true(length(sites) <= 2L, info = paste(sites, collapse = ", "))
})
