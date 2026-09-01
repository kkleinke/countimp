## B56/B57 -- the retry loop, and the three action values that must differ.
##
## Before this change countimp_check() was exported, documented and called by
## nothing; .countimp_max_redraw <- 3L had no reader; and action = "redraw"
## returned identical() results to action = "silent" (B57). A hard draw failure
## aborted the whole run -- at m = 5, maxit = 20 that discards 100 draws because
## of one.
##
## The contract these tests pin down:
##   1. a HARD failure is retried, up to .countimp_max_redraw times
##   2. a FLAGGED draw is never retried -- retrying large-but-finite draws
##      would bias the imputations towards the centre
##   3. healthy runs are bit-for-bit unchanged
##   4. the three action values behave differently

testthat::test_that("B56: a hard failure is retried and a later good draw wins", {
  n <- 0L
  draw_ <- function() {
    n <<- n + 1L
    if (n <= 2L) list(imp = c(1, NA, 3), fit = NULL, mu = c(1, 2, 3))
    else         list(imp = c(1, 2, 3),  fit = NULL, mu = c(1, 2, 3))
  }
  r <- ci(".countimp_draw_retry")(draw_, y_obs = c(0, 1, 2, 3, 4), method = "test")
  testthat::expect_identical(n, 3L)
  testthat::expect_identical(r, c(1, 2, 3))
})

testthat::test_that("B56: the loop stops after .countimp_max_redraw attempts", {
  m <- 0L
  draw_ <- function() { m <<- m + 1L; list(imp = c(1, NA), fit = NULL, mu = c(1, 2)) }
  testthat::expect_error(
    ci(".countimp_draw_retry")(draw_, y_obs = c(0, 1, 2), method = "test2"),
    "failed on 3 successive draws")
  testthat::expect_identical(m, 3L)
  ## the constant is the reader, not a hard-coded 3
  testthat::expect_identical(ci(".countimp_max_redraw"), 3L)
})

testthat::test_that("B56: max_redraw is honoured", {
  for (k in c(1L, 2L, 5L)) {
    m <- 0L
    draw_ <- function() { m <<- m + 1L; list(imp = c(1, NA), fit = NULL, mu = c(1, 2)) }
    testthat::expect_error(
      ci(".countimp_draw_retry")(draw_, y_obs = c(0, 1, 2), method = "k",
                                 max_redraw = k))
    testthat::expect_identical(m, k, info = paste("max_redraw =", k))
  }
})

testthat::test_that("B56: a flagged but finite draw is NOT retried", {
  ## Large finite draws are the tail of the target distribution, not a defect.
  ## Retrying them would quietly shrink the imputations towards the centre.
  n <- 0L
  draw_ <- function() {
    n <<- n + 1L
    list(imp = c(1, 5000), fit = NULL, mu = c(1, 5000))   # mu_extrapolated
  }
  r <- suppressWarnings(
    ci(".countimp_draw_retry")(draw_, y_obs = c(0, 1, 2), method = "flagged"))
  testthat::expect_identical(n, 1L)          # drawn once, kept
  testthat::expect_identical(r, c(1, 5000))
})

testthat::test_that("B56: an error inside the engine is retried, not propagated", {
  n <- 0L
  draw_ <- function() {
    n <<- n + 1L
    if (n < 3L) stop("countimp: synthetic engine failure", call. = FALSE)
    list(imp = c(2, 2), fit = NULL, mu = c(2, 2))
  }
  r <- ci(".countimp_draw_retry")(draw_, y_obs = c(1, 2, 3), method = "boom")
  testthat::expect_identical(n, 3L)
  testthat::expect_identical(r, c(2, 2))
})

testthat::test_that("B56: a persistent engine error reports its own reason once", {
  draw_ <- function() stop("countimp: synthetic persistent failure", call. = FALSE)
  err <- tryCatch(ci(".countimp_draw_retry")(draw_, y_obs = c(1, 2), method = "boom2"),
                  error = function(e) conditionMessage(e))
  testthat::expect_match(err, "failed on 3 successive draws")
  testthat::expect_match(err, "synthetic persistent failure")
  ## the nested "countimp:" prefix is stripped -- exactly one in the message
  testthat::expect_identical(lengths(regmatches(err, gregexpr("countimp:", err, fixed = TRUE))), 1L)
})

testthat::test_that("B56: every retry is recorded, for both failure branches", {
  for (fall in c("nonfinite", "error")) {
    countimp_diagnostics(enable = TRUE, reset = TRUE)
    draw_ <- if (identical(fall, "nonfinite"))
      function() list(imp = c(1, NA), fit = NULL, mu = c(1, 2))
    else function() stop("countimp: synthetic", call. = FALSE)
    invisible(tryCatch(ci(".countimp_draw_retry")(draw_, y_obs = c(0, 1, 2),
                                                  method = "rec"),
                       error = function(e) NULL))
    d <- countimp_diagnostics()
    testthat::expect_gte(sum(grepl("draw_repeated", d$problems)), 2L)
  }
  countimp_diagnostics(enable = FALSE, reset = TRUE)
})

testthat::test_that("B57: the three action values differ", {
  y <- c(0, 1, 2, 3, 4)
  hart <- c(1, NA, 3)          # draw_nonfinite -> hard
  weich <- c(1, 5000)          # mu_extrapolated via mu -> flagged, not hard

  ## silent never warns; warn always warns; redraw warns for flagged only,
  ## because a hard failure is about to be repeated by the caller.
  testthat::expect_silent(countimp_check(hart, y, method = "t", action = "silent"))
  testthat::expect_warning(countimp_check(hart, y, method = "t", action = "warn"))
  testthat::expect_silent(countimp_check(hart, y, method = "t", action = "redraw"))

  testthat::expect_warning(countimp_check(weich, y, mu = c(1, 5000),
                                          method = "t", action = "redraw"))

  ## and the returned objects are no longer identical between silent and redraw
  a <- countimp_check(weich, y, mu = c(1, 5000), method = "t", action = "silent")
  b <- suppressWarnings(countimp_check(weich, y, mu = c(1, 5000), method = "t",
                                       action = "redraw"))
  testthat::expect_identical(a$status, b$status)   # the verdict is the same ...
  testthat::expect_false(isTRUE(a$redraw))         # ... and neither asks for a redraw
})

testthat::test_that("B59: check_fit accepts countimp's own fit list", {
  ## .countimp_1l_fit() returns list(beta, cov, scale, theta) -- no class, no
  ## $coefficients. coef()/vcov() return NULL on it, which made check_fit
  ## report se_nonfinite for a sound fit and skip convergence entirely.
  gut <- list(beta = c(1, 0.5), cov = diag(2), scale = 1, theta = 2)
  testthat::expect_length(ci(".countimp_check_fit")(gut, "nb"), 0L)

  faelle <- list(
    coef_nonfinite        = list(beta = c(1, Inf), cov = diag(2), scale = 1, theta = 2),
    coef_extreme          = list(beta = c(1, 80),  cov = diag(2), scale = 1, theta = 2),
    se_nonfinite          = list(beta = c(1, 0.5), cov = diag(2), scale = -1, theta = 2),
    theta_degenerate_high = list(beta = c(1, 0.5), cov = diag(2), scale = 1, theta = 1e9))
  for (nm in names(faelle))
    testthat::expect_true(nm %in% ci(".countimp_check_fit")(faelle[[nm]], "nb"),
                          info = nm)
})

testthat::test_that("B56: healthy runs are unchanged by the retry loop", {
  ## The regression guard for the whole change: with no failures, the loop must
  ## draw exactly once per call and return what the engine returned.
  d <- zaehl_daten()
  ry <- !is.na(d$y)
  for (m in c("poisson", "quasipoisson", "nb")) {
    f <- get(paste0("mice.impute.", m), envir = ci_home())
    set.seed(77); a <- suppressWarnings(f(d$y, ry, d[, c("x1", "x2")]))
    set.seed(77); b <- suppressWarnings(f(d$y, ry, d[, c("x1", "x2")]))
    testthat::expect_identical(a, b, info = m)     # same seed, same draws
    testthat::expect_length(a, sum(!ry))
    testthat::expect_true(all(is.finite(a)), info = m)
  }
})

testthat::test_that("B56: an empty draw is legitimate when there is nothing to impute", {
  ## A variable with no missing values yields length-0 draws. That is correct,
  ## not a failure: no_draws must stay FLAGGED, never HARD, or the retry loop
  ## would repeat it three times and then abort a perfectly sound run.
  ## Complete cases only -- ry = TRUE on data that still contains NA would claim
  ## the missing values are observed and hand NA to the fit, which is a broken
  ## test rather than a legitimate empty draw.
  d <- zaehl_daten()
  d <- d[!is.na(d$y), , drop = FALSE]
  z <- suppressWarnings(
    ci(".countimp_1l_count")(d$y, rep(TRUE, nrow(d)), d[, c("x1", "x2")],
                             dist = "poisson", bayes = TRUE))
  testthat::expect_length(z, 0L)

  ## and the classification itself
  r <- suppressWarnings(countimp_check(numeric(0), d$y[!is.na(d$y)], method = "e"))
  testthat::expect_identical(r$status, "flagged")
  testthat::expect_false(isTRUE(r$redraw))
})

testthat::test_that("B56: warnings of the accepted draw reach the user", {
  ## .countimp_quietly() muffles everything so retries stay quiet. Wrapping the
  ## engines therefore silenced countimp's own deliberate warnings -- B04's
  ## underdispersion notice among them. The accepted draw's warnings must pass
  ## through; a discarded attempt's must not.
  n <- 0L
  draw_ <- function() {
    n <<- n + 1L
    warning("countimp: synthetic notice from attempt ", n, call. = FALSE)
    if (n <= 1L) list(imp = c(1, NA), fit = NULL, mu = c(1, 2))
    else         list(imp = c(1, 2),  fit = NULL, mu = c(1, 2))
  }
  w <- character(0)
  r <- withCallingHandlers(
    ci(".countimp_draw_retry")(draw_, y_obs = c(0, 1, 2), method = "warn"),
    warning = function(cond) {
      w <<- c(w, conditionMessage(cond)); invokeRestart("muffleWarning")
    })
  testthat::expect_identical(r, c(1, 2))
  testthat::expect_identical(sum(grepl("synthetic notice", w)), 1L)
  ## the one that got through belongs to attempt 2, the accepted draw
  testthat::expect_true(any(grepl("from attempt 2", w)))
  testthat::expect_false(any(grepl("from attempt 1", w)))
})
