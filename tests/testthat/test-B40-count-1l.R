## B40 -- single-level count family on one core, and the quasi-Poisson
## coefficient covariance.
##
## Ten exported names (poisson/pois/quasipoisson/qpois/nb, each with .boot)
## were six distinct bodies plus four byte-identical copies. They now delegate
## to .countimp_1l_count() with three switches.
##
## The substantive bug: the Bayesian variants drew beta from
## summary.glm(fit)$cov.unscaled for all three families. For quasi-Poisson
## with estimated dispersion phi the sampling covariance is phi *
## cov.unscaled, so the draw understated parameter uncertainty by sqrt(phi).
## Measured consequence: coverage of nominal 95% intervals for the slopes was
## 0.836 / 0.838 (MCSE 0.017) before the fix -- see
## analyse/k05_ueberdeckung_{VORHER,NACHHER}.csv.

## zaehl_daten() now lives in helper-internals.R -- B56 needs it too, and a
## function defined in a test file is not visible from another test file.

alle_10 <- c("poisson", "poisson.boot", "pois", "pois.boot",
             "quasipoisson", "quasipoisson.boot", "qpois", "qpois.boot",
             "nb", "nb.boot")

testthat::test_that("B40: all ten names delegate to the shared core", {
  for (m in alle_10) {
    f <- get(paste0("mice.impute.", m), envir = ci_home())
    src <- paste(deparse(body(f)), collapse = " ")
    ist_alias <- grepl("mice\\.impute\\.(poisson|quasipoisson)", src)
    testthat::expect_true(grepl(".countimp_1l_count", src, fixed = TRUE) || ist_alias,
                          info = m)
  }
})

testthat::test_that("B40: the aliases are thin and agree with their targets", {
  paare <- list(c("pois", "poisson"), c("pois.boot", "poisson.boot"),
                c("qpois", "quasipoisson"), c("qpois.boot", "quasipoisson.boot"))
  d <- zaehl_daten()
  for (p in paare) {
    a <- get(paste0("mice.impute.", p[1]), envir = ci_home())
    b <- get(paste0("mice.impute.", p[2]), envir = ci_home())
    ## the alias must call the target, not carry a second copy of the body
    testthat::expect_lte(length(deparse(body(a))), 6L)
    ## same random stream in, same imputations out
    set.seed(99); ia <- a(d$y, !is.na(d$y), d[, c("x1", "x2")])
    set.seed(99); ib <- b(d$y, !is.na(d$y), d[, c("x1", "x2")])
    testthat::expect_equal(ia, ib, info = paste(p, collapse = " vs "))
  }
})

testthat::test_that("B40: each switch value reaches the core", {
  gesehen <- list()
  echt <- ci(".countimp_1l_count")
  spion <- function(y, ry, x, wy = NULL, dist, bayes, EV = FALSE, ...) {
    gesehen[[length(gesehen) + 1L]] <<- list(dist = dist, bayes = bayes)
    echt(y, ry, x, wy = wy, dist = dist, bayes = bayes, EV = EV, ...)
  }
  ## via ci_spy(), not assign(): under R CMD check the namespace is sealed and
  ## the plain assignment errors with "cannot change value of locked binding"
  ci_spy(".countimp_1l_count", spion)

  d <- zaehl_daten()
  ry <- !is.na(d$y)
  for (m in alle_10)
    get(paste0("mice.impute.", m), envir = ci_home())(d$y, ry, d[, c("x1", "x2")])

  paare <- unique(vapply(gesehen, function(g) paste(g$dist, g$bayes), character(1)))
  testthat::expect_setequal(paare,
    c("poisson TRUE", "poisson FALSE", "quasipoisson TRUE",
      "quasipoisson FALSE", "negbin TRUE", "negbin FALSE"))
})

testthat::test_that("B40: imputations are non-negative whole numbers", {
  d <- zaehl_daten()
  ry <- !is.na(d$y)
  for (m in alle_10) {
    im <- get(paste0("mice.impute.", m), envir = ci_home())(d$y, ry, d[, c("x1", "x2")])
    testthat::expect_equal(length(im), sum(!ry), info = m)
    testthat::expect_true(all(is.finite(im)), info = m)
    testthat::expect_true(all(im >= 0), info = m)
    testthat::expect_true(all(im == round(im)), info = m)
  }
})

testthat::test_that("B40: a single predictor keeps its dimension", {
  ## x[ry, ] without drop = FALSE collapses to a vector for one predictor.
  ## cbind(1, ...) used to hide this; the core is explicit about it.
  d <- zaehl_daten()
  ry <- !is.na(d$y)
  for (m in alle_10) {
    im <- get(paste0("mice.impute.", m), envir = ci_home())(d$y, ry, d[, "x1", drop = FALSE])
    testthat::expect_equal(length(im), sum(!ry), info = m)
    testthat::expect_true(all(is.finite(im)), info = m)
  }
})

testthat::test_that("B40: quasi-Poisson scales the coefficient covariance by phi", {
  ## The regression guard. Fit a deliberately overdispersed model, then check
  ## that the reported scale equals the dispersion -- not 1.
  set.seed(11); n <- 400
  x <- cbind(1, x1 = stats::rnorm(n), x2 = stats::rnorm(n))
  y <- MASS::rnegbin(n, mu = exp(1 + 0.6 * x[, 2] - 0.4 * x[, 3]), theta = 1.5)

  fq <- ci(".countimp_1l_fit")(x, y, "quasipoisson")
  testthat::expect_gt(fq$scale, 1.5)                 # clearly overdispersed
  testthat::expect_equal(fq$scale, fq$dispersion)

  ## Poisson and negbin fix phi = 1 -- scaling must stay neutral there.
  fp <- ci(".countimp_1l_fit")(x, y, "poisson")
  testthat::expect_equal(fp$scale, 1)
  fn <- ci(".countimp_1l_fit")(x, y, "negbin")
  testthat::expect_equal(fn$scale, 1)

  ## and the draw actually uses it: sd of the drawn beta must grow by sqrt(phi)
  set.seed(3)
  draw_ <- function(scale) {
    replicate(400, ci(".countimp_1l_draw_beta")(fq$beta, fq$cov, scale, TRUE)[2])
  }
  s1 <- stats::sd(draw_(1)); sp <- stats::sd(draw_(fq$scale))
  testthat::expect_equal(sp / s1, sqrt(fq$scale), tolerance = 0.12)
})

testthat::test_that("B40: a non-positive dispersion stops with an explanation", {
  b <- c(1, 0.5); cv <- diag(2) * 0.01
  testthat::expect_error(ci(".countimp_1l_draw_beta")(b, cv, 0, TRUE),
                         "dispersion estimate")
  testthat::expect_error(ci(".countimp_1l_draw_beta")(b, cv, NA_real_, TRUE),
                         "dispersion estimate")
  ## singular covariance
  sing <- matrix(1, 2, 2)
  testthat::expect_error(ci(".countimp_1l_draw_beta")(b, sing, 1, TRUE),
                         "positive definite")
  ## bootstrap route ignores all of it and returns the estimate
  testthat::expect_equal(ci(".countimp_1l_draw_beta")(b, sing, 0, FALSE), b)
})

testthat::test_that("B40: the glm summary wrapper has a working fallback", {
  set.seed(2); n <- 150
  x <- cbind(1, x1 = stats::rnorm(n))
  y <- stats::rpois(n, exp(1 + 0.5 * x[, 2]))
  fit <- stats::glm.fit(x, y, family = stats::poisson(link = "log"))

  ref <- ci(".countimp_glm_summary")(fit, 2L)
  testthat::expect_equal(ref$dispersion, 1)

  ## the fallback path must reproduce cov.unscaled from the QR factor
  cu <- chol2inv(fit$qr$qr[1:2, 1:2, drop = FALSE])
  testthat::expect_equal(unname(cu), unname(ref$cov.unscaled), tolerance = 1e-8)
})

testthat::test_that("B40: foreign calls are confined to the core", {
  pfad <- file.path(if (dir.exists("../../R")) "../../R" else "R")
  testthat::skip_if_not(dir.exists(pfad), "package sources not locatable")
  src <- sub("#.*$", "", readLines(file.path(pfad, "mice.impute.poisson.R"),
                                  warn = FALSE))
  testthat::expect_false(any(grepl("MASS::|summary\\.glm|glm\\.fit", src)))
  kern <- sub("#.*$", "", readLines(file.path(pfad, "impute1lcount.R"), warn = FALSE))
  testthat::expect_true(any(grepl("MASS::glm\\.nb", kern)))
})
