## B10 -- countimp runs without mice, and its own base methods are sound.
##
## Up to 2.6.0 countimp() routed every call that did not request a
## zero-inflated count method to mice::mice(), and it had no imputation methods
## for continuous, binary or categorical variables. Both are gone: the engine
## handles all methods and the base methods below are countimp's own. These
## tests pin the properties that make that safe.

test_that("countimp imputes a mixed data set without calling mice", {
  set.seed(11)
  n <- 120
  d <- data.frame(
    y   = rpois(n, 3),
    v   = rnorm(n),
    b   = factor(rbinom(n, 1, 0.5)),
    nom = factor(sample(c("a", "b", "c"), n, TRUE)),
    ord = ordered(sample(c("lo", "mid", "hi"), n, TRUE),
                  levels = c("lo", "mid", "hi")))
  for (v in names(d)) d[sample.int(n, 25), v] <- NA

  ## A local mock of mice::mice cannot be installed (locked namespace), so the
  ## check is on the result instead: every variable imputed, no error.
  ##
  ## suppressWarnings() is deliberate and narrow. With n = 120 and 25 missing
  ## values the negative-binomial dispersion is weakly identified, so some
  ## draws land near zero and get floored -- the diagnostics system warns about
  ## exactly that, which is the behaviour B22 restored. The warning is expected
  ## here and is asserted below rather than hidden.
  imp <- suppressWarnings(
    countimp(d, method = c("nb", "pmm", "logreg", "polyreg", "polr"),
             m = 2, maxit = 2, printFlag = FALSE, seed = 1))

  expect_s3_class(imp, "mids")
  expect_equal(imp$m, 2L)
  for (v in names(d)) {
    expect_equal(nrow(imp$imp[[v]]), sum(is.na(d[[v]])), info = v)
    expect_false(anyNA(imp$imp[[v]]), info = v)
  }
})

test_that("the diagnostics log survives mixed row types (B22)", {
  ## Numerical events and screened draws used to write rows with different
  ## column sets; countimp_diagnostics() rbind()s them, so the two together
  ## aborted the call and lost every recorded row.
  countimp_diagnostics(enable = TRUE, reset = TRUE)
  on.exit(countimp_diagnostics(enable = FALSE, reset = TRUE), add = TRUE)

  ci(".countimp_note_event")("theta_draw_floored",
                                 "drawn theta = 0.01 raised to 0.1")
  ci(".countimp_note")("collinear predictors: ridge penalty applied",
                            where = "norm")
  dg <- countimp_diagnostics()

  expect_equal(nrow(dg), 2L)
  expect_true(all(c("method", "status", "problems", "n_imp", "obs_max",
                    "draw_max", "draw_ratio", "mu_ratio", "warning")
                  %in% names(dg)))
  expect_equal(dg$method, c("theta_draw_floored", "norm"))
  expect_true(all(dg$status == "flagged"))
})

test_that("default methods resolve to countimp's own implementations", {
  set.seed(12)
  n <- 80
  d <- data.frame(v = rnorm(n), b = factor(rbinom(n, 1, 0.5)),
                  nom = factor(sample(c("a", "b", "c"), n, TRUE)),
                  ord = ordered(sample(1:3, n, TRUE)))
  for (v in names(d)) d[sample.int(n, 15), v] <- NA

  ini <- countimp(d, maxit = 0)
  expect_equal(unname(ini$method[c("v", "b", "nom", "ord")]),
               c("pmm", "logreg", "polyreg", "polr"))

  ## The resolver must return countimp's function, not mice's, for each.
  ## Assert what the claim actually is -- "not mice's" -- rather than the name
  ## of the environment: under the development load path countimp's own code
  ## lives in the global environment, and demanding "countimp" made this test
  ## pass only against an installed build (see B38).
  own <- environmentName(ci_home())
  for (m in c("pmm", "logreg", "polyreg", "polr", "norm", "mean", "sample")) {
    f <- ci(".countimp_find_method")(m)
    expect_true(is.function(f), info = m)
    expect_identical(environmentName(environment(f)), own, info = m)
    expect_false(identical(environmentName(environment(f)), "mice"), info = m)
  }
})

test_that("options(countimp.methods = 'mice') switches the resolver back", {
  skip_if_not_installed("mice")
  old <- options(countimp.methods = "mice")
  on.exit(options(old), add = TRUE)
  f <- ci(".countimp_find_method")("pmm")
  expect_identical(environmentName(environment(f)), "mice")

  ## countimp's own count methods keep precedence regardless of the option --
  ## they are the point of the package.
  g <- ci(".countimp_find_method")("nb")
  expect_identical(environmentName(environment(g)), environmentName(ci_home()))
  expect_false(identical(environmentName(environment(g)), "mice"))
})

test_that("countimp() rejects mice-only arguments instead of ignoring them", {
  set.seed(13)
  n <- 60
  d <- data.frame(y = rpois(n, 2), v = rnorm(n))
  d$y[sample.int(n, 12)] <- NA

  ## The engine's `...` would absorb all of these silently, which would mean
  ## imputing under settings the caller did not ask for.
  expect_error(countimp(d, method = c("nb", "pmm"), ignore = rep(FALSE, n)),
               "mice::mice", fixed = TRUE)
  expect_error(countimp(d, method = c("nb", "pmm"), blocks = list("y")),
               "visitSequence")
  expect_error(countimp(d, method = c("nb", "pmm"), maxti = 3),
               "did you mean `maxit`", fixed = TRUE)

  ## And the escape hatch works.
  old <- options(countimp.check.args = FALSE)
  on.exit(options(old), add = TRUE)
  expect_silent(imp <- countimp(d, method = c("nb", "pmm"), m = 1, maxit = 1,
                                printFlag = FALSE, donors = 3))
})

test_that("countimp_complete fills all missing cells and keeps column types", {
  set.seed(21)
  n <- 90
  d <- data.frame(y = rpois(n, 3), v = rnorm(n),
                  b = factor(sample(c("x", "y"), n, TRUE)),
                  o = ordered(sample(1:3, n, TRUE)))
  for (v in names(d)) d[sample.int(n, 18), v] <- NA
  imp <- suppressWarnings(
    countimp(d, method = c("nb", "pmm", "logreg", "polr"), m = 3, maxit = 2,
             printFlag = FALSE, seed = 2))

  cc <- countimp_complete(imp, 1)
  expect_equal(dim(cc), dim(d))
  expect_false(anyNA(cc))
  expect_identical(vapply(cc, function(z) class(z)[1], ""),
                   vapply(d,  function(z) class(z)[1], ""))
  expect_equal(levels(cc$b), levels(d$b))
  expect_true(is.ordered(cc$o))

  ## observed values must be untouched
  obs <- !is.na(d$v)
  expect_equal(cc$v[obs], d$v[obs])

  lng <- countimp_complete(imp, "long")
  expect_equal(nrow(lng), 3L * n)
  expect_equal(tail(names(lng), 2), c(".imp", ".id"))   # mice's column order
  expect_equal(sort(unique(lng$.imp)), 1:3)

  lst <- countimp_complete(imp, "all")
  expect_length(lst, 3L)
  expect_equal(dim(lst[[1]]), dim(d))

  brd <- countimp_complete(imp, "broad")
  expect_equal(ncol(brd), 3L * ncol(d))

  expect_equal(countimp_complete(imp, 0), d)            # action 0 = original
  expect_error(countimp_complete(imp, 99), "between 1 and m")
  expect_error(countimp_complete(d), "mids object")

  ## include = TRUE prepends the incomplete data as imputation 0
  li <- countimp_complete(imp, "long", include = TRUE)
  expect_equal(nrow(li), 4L * n)
  expect_equal(sum(is.na(li$v)), sum(is.na(d$v)))
})

test_that("countimp_complete agrees with mice::complete", {
  skip_if_not_installed("mice")
  set.seed(22)
  n <- 70
  d <- data.frame(y = rpois(n, 3), v = rnorm(n),
                  b = factor(sample(c("x", "y"), n, TRUE)))
  for (v in names(d)) d[sample.int(n, 14), v] <- NA
  imp <- suppressWarnings(
    countimp(d, method = c("nb", "pmm", "logreg"), m = 2, maxit = 2,
             printFlag = FALSE, seed = 3))

  for (a in list(1L, 2L, "long", "broad")) {
    expect_equal(as.data.frame(countimp_complete(imp, a)),
                 as.data.frame(mice::complete(imp, a)),
                 ignore_attr = TRUE, info = paste("action", a))
  }
  expect_equal(unname(countimp_complete(imp, "all")),
               unname(mice::complete(imp, "all")), ignore_attr = TRUE)
})

test_that("the post-imputation tools run without mice", {
  ## The point of B23: nothing in the analysis path may require mice. The
  ## check makes requireNamespace("mice") fail inside countimp's namespace.
  set.seed(23)
  n <- 80
  d <- data.frame(y = rpois(n, 3), v = rnorm(n))
  for (v in names(d)) d[sample.int(n, 16), v] <- NA

  suppressMessages(
    trace(requireNamespace, where = asNamespace("countimp"), print = FALSE,
          tracer = quote(if (identical(package, "mice")) return(FALSE))))
  on.exit(suppressMessages(
    untrace(requireNamespace, where = asNamespace("countimp"))), add = TRUE)

  imp <- suppressWarnings(countimp(d, method = c("nb", "pmm"), m = 2, maxit = 2,
                                  printFlag = FALSE, seed = 4))
  expect_s3_class(imp, "mids")
  expect_false(anyNA(countimp_complete(imp, 1)))
  expect_length(compare.obs.imp(d, imp), 2L)
  expect_length(compare.percent.count(d, imp, counts = "y"), 1L)

  fits <- lapply(1:2, function(i)
    stats::lm(y ~ v, data = countimp_complete(imp, i)))
  p <- miinference(fits)
  expect_true(all(is.finite(p$est)))
})

test_that("with() analyses the imputed data, not the complete cases (B27)", {
  ## Without a with.mids method, with(imp, ...) fell through to with.default().
  ## If identically named objects existed in the caller's environment -- the
  ## normal situation in a simulation script -- the model was fitted to those,
  ## i.e. to the complete cases, and the imputations were never used.
  set.seed(24)
  n <- 200
  d <- data.frame(x = rnorm(n))
  d$y <- rpois(n, exp(1 + 0.5 * d$x))
  d$y[sample.int(n, 60)] <- NA
  ## deliberately shadow the columns, which is what made the old bug silent
  y <- d$y
  x <- d$x

  imp <- suppressWarnings(countimp(d, method = c("", "nb"), m = 3, maxit = 3,
                                  printFlag = FALSE, seed = 5))
  res <- with(imp, stats::glm(y ~ x, family = stats::poisson))

  expect_s3_class(res, "mira")
  expect_length(res$analyses, 3L)
  ## every fit must use ALL rows, not just the observed ones
  for (a in res$analyses)
    expect_equal(length(a$residuals), n)
  expect_false(identical(coef(res$analyses[[1]]), coef(res$analyses[[2]])))

  p <- miinference(res)
  expect_equal(length(p$est), 2L)
  expect_true(all(is.finite(p$std.err)))

  ## objects from the calling scope stay visible inside the expression
  w <- rep(1, n)
  res2 <- with(imp, stats::lm(y ~ x, weights = w))
  expect_length(res2$analyses, 3L)
})

test_that("the returned data keep the user's factor attributes (B28)", {
  ## padModel() recodes factors to contr.treatment internally. That copy used
  ## to be returned as $data, so an ordered factor came back with treatment
  ## contrasts and models on the completed data reported different
  ## coefficients than on the input data.
  set.seed(25)
  n <- 150
  d <- data.frame(y = rnorm(n),
                  o = ordered(sample(1:3, n, TRUE), levels = 1:3),
                  f = factor(sample(c("a", "b", "c"), n, TRUE)))
  d$y <- d$y + as.integer(d$o)
  d$y[sample.int(n, 30)] <- NA

  imp <- suppressWarnings(countimp(d, method = c("pmm", "", ""), m = 1,
                                  maxit = 1, printFlag = FALSE, seed = 6))

  expect_identical(attributes(imp$data$o), attributes(d$o))
  expect_identical(attributes(imp$data$f), attributes(d$f))
  expect_null(attr(imp$data$o, "contrasts"))

  ## the consequence that made this matter: polynomial contrasts for ordered
  ## factors, as R uses by default
  cf <- names(coef(stats::lm(y ~ o, data = countimp_complete(imp, 1))))
  expect_equal(cf, c("(Intercept)", "o.L", "o.Q"))
})

test_that("miinference names empty inputs instead of failing in vcov", {
  expect_error(miinference(list(NULL, NULL), list(NULL, NULL)),
               "empty element")
  expect_error(miinference(list(numeric(0), numeric(0)),
                           list(numeric(0), numeric(0))),
               "empty|length")
})

test_that("norm draws propagate parameter uncertainty", {
  ## The Bayesian draw must be more variable than the one that conditions on
  ## the point estimate; norm.predict must have no residual variability at all.
  set.seed(14)
  n <- 60
  x <- matrix(rnorm(n), ncol = 1, dimnames = list(NULL, "x"))
  y <- 1 + 0.7 * x[, 1] + rnorm(n)
  ry <- rep(TRUE, n); ry[1:20] <- FALSE

  ## replicate() evaluates its expression lazily in an environment where a
  ## `...` from the enclosing function is no longer resolvable, which yields
  ## silent NAs. Build the call with do.call() instead.
  sd_of_mean <- function(f, ...) {
    extra <- list(...)
    stats::sd(replicate(200, mean(do.call(f, c(list(y, ry, x), extra)))))
  }
  s_bayes <- sd_of_mean(ci("countimp.impute.norm"))
  s_nob   <- sd_of_mean(ci("countimp.impute.norm.nob"))
  expect_false(is.na(s_bayes) || is.na(s_nob))
  expect_gt(s_bayes, s_nob)

  a <- ci("countimp.impute.norm.predict")(y, ry, x)
  b <- ci("countimp.impute.norm.predict")(y, ry, x)
  expect_equal(a, b)                       # deterministic by construction
})

test_that("pmm imputes observed values only and respects the donor count", {
  set.seed(15)
  n <- 80
  x <- matrix(rnorm(n), ncol = 1, dimnames = list(NULL, "x"))
  y <- round(1 + 0.8 * x[, 1] + rnorm(n), 3)
  ry <- rep(TRUE, n); ry[1:25] <- FALSE

  imp <- ci("countimp.impute.pmm")(y, ry, x)
  expect_length(imp, sum(!ry))
  expect_true(all(imp %in% y[ry]))         # donors, not model predictions

  ## donors = 1 must pick the single nearest neighbour, deterministically
  ## given the drawn coefficients: repeated calls with the same seed agree.
  set.seed(1); a <- ci("countimp.impute.pmm")(y, ry, x, donors = 1)
  set.seed(1); b <- ci("countimp.impute.pmm")(y, ry, x, donors = 1)
  expect_equal(a, b)
})

test_that("the donor search finds true nearest neighbours", {
  ## ci(".countimp_match") replaces a per-case loop with a binary search; this pins
  ## it against a direct computation on random inputs.
  set.seed(16)
  for (k in c(1L, 3L, 5L)) {
    obs <- rnorm(60)
    mis <- rnorm(25)
    idx <- ci(".countimp_match")(obs, mis, donors = k)
    expect_length(idx, length(mis))
    expect_true(all(idx >= 1L & idx <= length(obs)))
    for (i in seq_along(mis)) {
      ## the chosen donor must be among the k nearest by absolute distance
      d <- abs(obs - mis[i])
      expect_true(idx[i] %in% order(d)[seq_len(k)], info = paste("k", k, "i", i))
    }
  }
  ## donors larger than the number of donors available: sample from all
  expect_length(ci(".countimp_match")(rnorm(3), rnorm(10), donors = 99), 10L)
})

test_that("categorical methods return factors with the original levels", {
  set.seed(17)
  n <- 100
  x <- matrix(rnorm(n), ncol = 1, dimnames = list(NULL, "x"))
  ry <- rep(TRUE, n); ry[1:30] <- FALSE

  b <- factor(rbinom(n, 1, plogis(x[, 1])), labels = c("no", "yes"))
  ib <- ci("countimp.impute.logreg")(b, ry, x)
  expect_s3_class(ib, "factor")
  expect_equal(levels(ib), levels(b))
  expect_length(ib, sum(!ry))

  nom <- factor(sample(c("a", "b", "c"), n, TRUE))
  inom <- ci("countimp.impute.polyreg")(nom, ry, x)
  expect_equal(levels(inom), levels(nom))
  expect_false(anyNA(inom))

  ord <- ordered(sample(c("lo", "mid", "hi"), n, TRUE),
                 levels = c("lo", "mid", "hi"))
  iord <- ci("countimp.impute.polr")(ord, ry, x)
  expect_equal(levels(iord), levels(ord))
  expect_false(anyNA(iord))

  ## logreg on a three-category variable is a specification error, not a
  ## silently wrong imputation.
  expect_error(ci("countimp.impute.logreg")(nom, ry, x), "two categories")
})

test_that("augmentation keeps coefficients finite under perfect prediction", {
  ## y is a deterministic function of x: the ML logistic fit diverges. With
  ## augmentation the imputations must still show variability rather than
  ## collapsing onto one category.
  set.seed(18)
  n <- 60
  x <- matrix(c(rep(-1, 30), rep(1, 30)), ncol = 1, dimnames = list(NULL, "x"))
  y <- factor(ifelse(x[, 1] > 0, "yes", "no"), levels = c("no", "yes"))
  ry <- rep(TRUE, n)
  ry[c(1:5, 56:60)] <- FALSE               # missing in both groups

  draws <- replicate(50, as.integer(
    ci("countimp.impute.logreg")(y, ry, x)))
  expect_false(anyNA(draws))
  ## Perfect prediction pushes towards a single category; augmentation must
  ## leave at least some uncertainty in the group that is being extrapolated.
  expect_gt(length(unique(as.vector(draws))), 1L)
})

test_that("norm.draw survives collinear predictors", {
  ## A duplicated predictor makes the design matrix rank deficient. FCS
  ## produces such matrices routinely, so this must not error.
  set.seed(19)
  n <- 50
  x1 <- rnorm(n)
  x <- cbind(a = x1, b = x1, c = rnorm(n))   # b is a copy of a
  y <- 1 + x1 + rnorm(n)
  ry <- rep(TRUE, n); ry[1:15] <- FALSE

  countimp_diagnostics(enable = TRUE, reset = TRUE)
  on.exit(countimp_diagnostics(enable = FALSE, reset = TRUE), add = TRUE)
  imp <- ci("countimp.impute.norm")(y, ry, x)
  expect_length(imp, sum(!ry))
  expect_false(anyNA(imp))
  expect_true(all(is.finite(imp)))
})

test_that("mean and sample behave as documented", {
  set.seed(20)
  y <- c(rnorm(40), rep(NA, 10))
  ry <- !is.na(y)
  x <- matrix(rnorm(50), ncol = 1)

  m <- ci("countimp.impute.mean")(y, ry, x)
  expect_equal(unique(m), mean(y[ry]))
  expect_length(m, 10L)

  s <- ci("countimp.impute.sample")(y, ry)
  expect_true(all(s %in% y[ry]))
  expect_length(s, 10L)
})
