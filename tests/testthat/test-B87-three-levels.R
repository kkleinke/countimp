## B87: three (and more) grouping levels through the formula interface
##
## What these tests defend, in the order in which each was easy to get wrong:
##
##  1. The per-level variance. .countimp_tau() read the FIRST random-effect
##     term whatever level was asked for. With school SD 0.42 and class SD 0.20
##     that draws an unseen class effect from the school distribution -- twice
##     too wide, and invisible in the imputations, which stay plausible counts.
##  2. Partial novelty. A unit can sit in a KNOWN school and an UNKNOWN class.
##     Its school effect is estimated and must be kept; only the class effect is
##     drawn. The old code replaced the whole prediction by the population
##     average, which throws the school effect away.
##  3. The single-level path must be untouched. The rewrite generalises it (an
##     unseen effect is zero under allow.new.levels, so adding a draw to the
##     link-scale prediction IS the old "population average plus draw"), and
##     B20/B34 keep measuring that -- here we only check the two calling
##     conventions for `obs_levels` agree.
##  4. What the type codes cannot carry: a slope on an inner level, or an inner
##     level without an intercept. Both are refused by name rather than silently
##     dropped -- a dropped random term does not fail, it just shrinks the
##     between-cluster variance of the imputations.
##  5. The cluster bootstrap is NOT extended yet. It must say so instead of
##     resampling the outer level alone, which would treat nested classes as
##     independent.

skip_if_not_installed("glmmTMB")

## Pupils in classes in schools. Sizes kept small: three-level fits run in
## R CMD check.
.b87_daten <- function(nsch = 10L, nclass = 4L, nsub = 10L, seed = 87L,
                       tau_s = 0.40, tau_k = 0.25, nmis = 90L) {
  set.seed(seed)
  school <- rep(seq_len(nsch), each = nclass * nsub)
  classroom <- rep(seq_len(nsch * nclass), each = nsub)
  n <- length(school)
  u_s <- stats::rnorm(nsch, 0, tau_s)
  u_k <- stats::rnorm(nsch * nclass, 0, tau_k)
  x1 <- stats::rnorm(n)
  y <- stats::rpois(n, exp(0.9 + 0.5 * x1 + u_s[school] + u_k[classroom]))
  d <- data.frame(y = y, x1 = x1, school = school, classroom = classroom)
  d$wahr <- y
  d$y[sample.int(n, nmis)] <- NA
  d
}

.b87_fit <- function(d) {
  o <- d[!is.na(d$y), ]
  glmmTMB::glmmTMB(y ~ x1 + (1 | school) + (1 | classroom),
                   data = data.frame(y = o$y, x1 = o$x1,
                                     school = factor(o$school),
                                     classroom = factor(o$classroom)),
                   family = stats::poisson())
}


## --- 1: the compiler ---------------------------------------------------------

test_that("B87: two random-effect terms compile to two class codes", {
  d <- .b87_daten()[, c("y", "x1", "school", "classroom")]
  sp <- suppressMessages(countimp_spec(d, y ~ x1 + (1 | school) + (1 | classroom),
                                       family = poisson_count()))
  expect_identical(unname(sp$type$y[c("school", "classroom")]), c(-2L, -2L))
  expect_identical(unname(sp$type$y[["x1"]]), 1L)
  ## the method name is unchanged: the levels live in the type codes, so no new
  ## export is needed for a three-level model
  expect_identical(unname(sp$method[["y"]]), "2l.poisson")

  dec <- ci(".countimp_decode_type")(sp$type$y, names(sp$type$y))
  expect_identical(dec$group, c("school", "classroom"))
  ## and the formula the engine builds carries both. Compared without spaces:
  ## deparse() renders the bar as "1 | school", the builder writes "1|school",
  ## and the test is about the terms, not about either one's spacing.
  form <- gsub(" ", "", deparse1(ci(".countimp_2l_formula")(dec, "count",
                                                            response = "Y")))
  expect_match(form, "(1|school)", fixed = TRUE)
  expect_match(form, "(1|classroom)", fixed = TRUE)
})

test_that("B87: random slopes sit on the first level only", {
  d <- .b87_daten()[, c("y", "x1", "school", "classroom")]
  sp <- suppressMessages(countimp_spec(d,
          y ~ x1 + (1 + x1 | school) + (1 | classroom), family = poisson_count()))
  dec <- ci(".countimp_decode_type")(sp$type$y, names(sp$type$y))
  form <- gsub(" ", "", deparse1(ci(".countimp_2l_formula")(dec, "count",
                                                            response = "Y")))
  expect_match(form, "(1+x1|school)", fixed = TRUE)
  expect_match(form, "(1|classroom)", fixed = TRUE)

  ## what the type codes cannot carry is refused, not dropped
  expect_error(suppressMessages(countimp_spec(d,
                 y ~ x1 + (1 | school) + (1 + x1 | classroom),
                 family = poisson_count())),
               "Only the first random-effect term")
  expect_error(suppressMessages(countimp_spec(d,
                 y ~ x1 + (1 | school) + (0 + x1 | classroom),
                 family = poisson_count())),
               "Only the first random-effect term")
  ## the same factor twice is a writing error, not a second level
  expect_error(suppressMessages(countimp_spec(d,
                 y ~ x1 + (1 | school) + (0 + x1 | school),
                 family = poisson_count())),
               "more than one random-effect term")
})


## --- 2: the per-level variance ----------------------------------------------

test_that("B87: tau is read per grouping level, not always the first", {
  d <- .b87_daten()
  f <- .b87_fit(d)
  t_s <- ci(".countimp_tau")(f, "school")
  t_k <- ci(".countimp_tau")(f, "classroom")
  ## the two levels differ enough that confusing them matters
  expect_gt(t_s / t_k, 1.4)
  ## no name asked for: the first, as every single-level caller expects
  expect_equal(ci(".countimp_tau")(f), t_s, tolerance = 1e-12)
  ## a name the model does not have is a programming error, said by name
  expect_error(ci(".countimp_tau")(f, "gibtsnicht"), "no random-effect term")
})


## --- 3: unseen levels, per level --------------------------------------------

test_that("B87: a new class in a known school keeps the school effect", {
  d <- .b87_daten()
  f <- .b87_fit(d)
  o <- d[!is.na(d$y), ]
  nd <- data.frame(x1 = d$x1, school = factor(d$school), classroom = factor(d$classroom))

  ## hide one class, keep its school
  weg <- 3L
  lev <- list(school = unique(as.character(o$school)),
              classroom = setdiff(unique(as.character(o$classroom)), as.character(weg)))
  ziel <- which(d$classroom == weg)
  expect_gt(length(ziel), 5L)

  set.seed(11)
  eta <- as.numeric(ci(".countimp_predict_2l")(f, nd, type = "link",
                                               grp = c("school", "classroom"),
                                               obs_levels = lev))
  ## The school effect must still be in there. Reference: the prediction with
  ## BOTH levels known, minus the class effect that the fit holds for this
  ## class -- i.e. what remains is school + fixed part.
  roh <- as.numeric(stats::predict(f, newdata = nd, type = "link",
                                   allow.new.levels = TRUE))
  ## rows of other classes are untouched
  expect_equal(eta[-ziel], roh[-ziel], tolerance = 1e-10)
  ## the hidden class is shifted by ONE draw, the same for all its rows
  versatz <- eta[ziel] - roh[ziel]
  expect_equal(length(unique(round(versatz, 10))), 1L)
  expect_true(abs(versatz[1L]) > 0)

  ## and the shift is of the order of tau_klasse, not tau_schule: 200 draws,
  ## compared against both. Drawing from the wrong level is a factor ~2 error.
  set.seed(12)
  vs <- replicate(200, {
    e <- as.numeric(ci(".countimp_predict_2l")(f, nd[ziel[1L], , drop = FALSE],
           type = "link", grp = c("school", "classroom"), obs_levels = lev))
    e - roh[ziel[1L]]
  })
  t_k <- ci(".countimp_tau")(f, "classroom")
  t_s <- ci(".countimp_tau")(f, "school")
  expect_lt(abs(stats::sd(vs) - t_k), abs(stats::sd(vs) - t_s))
})

test_that("B87: both calling conventions for obs_levels agree", {
  d <- .b87_daten()
  f <- .b87_fit(d)
  o <- d[!is.na(d$y), ]
  nd <- data.frame(x1 = d$x1[1:20], school = factor(d$school[1:20]),
                   classroom = factor(d$classroom[1:20]))
  lev <- unique(as.character(o$school))
  set.seed(5)
  a <- as.numeric(ci(".countimp_predict_2l")(f, nd, type = "link",
                                             grp = "school", obs_levels = lev))
  set.seed(5)
  b <- as.numeric(ci(".countimp_predict_2l")(f, nd, type = "link",
                     grp = "school", obs_levels = list(school = lev)))
  expect_equal(a, b, tolerance = 1e-12)
  ## a bare vector with several levels cannot be attributed, so it is refused
  expect_error(ci(".countimp_predict_2l")(f, nd, type = "link",
                 grp = c("school", "classroom"), obs_levels = lev),
               "one entry per grouping level")
})


## --- 4: the imputation itself ------------------------------------------------

test_that("B87: three-level imputation runs and keeps both variance components", {
  d <- .b87_daten()
  dd <- d[, c("y", "x1", "school", "classroom")]
  imp <- suppressWarnings(countimp(dd,
           formulas = list(y ~ x1 + (1 | school) + (1 | classroom)),
           family = poisson_count(), m = 2, maxit = 1, printFlag = FALSE))
  v <- unlist(imp$imp$y)
  expect_equal(length(v), 2L * sum(is.na(dd$y)))
  expect_false(anyNA(v))
  expect_true(all(v >= 0))

  ## The imputations must not flatten the inner level. Measured on the
  ## completed data against the full data: the class SD survives.
  sd_kl <- function(dat) {
    f <- suppressWarnings(glmmTMB::glmmTMB(y ~ x1 + (1 | school) + (1 | classroom),
           data = transform(dat, school = factor(school), classroom = factor(classroom)),
           family = stats::poisson()))
    sqrt(glmmTMB::VarCorr(f)$cond$classroom[1L, 1L])
  }
  voll <- sd_kl(data.frame(y = d$wahr, x1 = d$x1, school = d$school,
                           classroom = d$classroom))
  fertig <- mean(vapply(1:2, function(i) sd_kl(countimp_complete(imp, i)),
                        numeric(1)))
  expect_gt(voll, 0.1)                       # the design really has that level
  expect_gt(fertig, 0.5 * voll)              # and it is not flattened away
})

test_that("B87: the two model parts must name the same grouping levels", {
  ## The type codes mark a variable as grouping (-2) without saying for which
  ## model part, so a zero formula with FEWER (1 | g) terms cannot be expressed.
  ## Measured before this check: `zero = ~ z1 + (1 | id)` beside
  ## `y ~ x1 + (1 | id) + (1 | kl)` produced the zero-part formula
  ## `Y ~ z1 + (1 | id) + (1 | kl)` -- a random effect the user did not write,
  ## invisible in the output. Same shape as B53.
  ##
  ## Only possible since three levels exist: with one grouping factor the two
  ## parts could not differ.
  set.seed(875); n <- 400L
  d <- data.frame(x1 = rnorm(n), z1 = rnorm(n),
                  id = rep(seq_len(20), each = 20L),
                  kl = rep(seq_len(40), each = 10L))
  d$y <- ifelse(runif(n) < 0.3, 0L, rpois(n, 3))
  d$y[sample.int(n, 80L)] <- NA

  expect_error(suppressMessages(countimp_spec(d, y ~ x1 + (1 | id) + (1 | kl),
                 zero = ~ z1 + (1 | id), family = zi_nb())),
               "different grouping levels")
  ## the other direction too
  expect_error(suppressMessages(countimp_spec(d, y ~ x1 + (1 | id),
                 zero = ~ z1 + (1 | id) + (1 | kl), family = zi_nb())),
               "different grouping levels")

  ## same levels: fine, and the zero part keeps its OWN predictors
  sp <- suppressMessages(countimp_spec(d, y ~ x1 + (1 | id) + (1 | kl),
          zero = ~ z1 + (1 | id) + (1 | kl), family = zi_nb()))
  expect_identical(unname(sp$method[["y"]]), "2l.zinb")
  expect_identical(unname(sp$type$y[["x1"]]), 3L)   # count part only
  expect_identical(unname(sp$type$y[["z1"]]), 5L)   # zero part only

  ## omitting `zero` reuses the count part, so the levels agree by construction
  expect_silent(suppressMessages(countimp_spec(d, y ~ x1 + (1 | id) + (1 | kl),
                                               family = zi_nb())))
})


## --- 5: the hierarchical cluster bootstrap ----------------------------------

test_that("B87: the bootstrap resamples whole schools, classes included", {
  set.seed(871)
  dat <- data.frame(Y = 1:40, school = rep(1:4, each = 10),
                    classroom = rep(1:8, each = 5))
  out <- ci(".countimp_boot_clusters")(dat, c("school", "classroom"))
  ## every school in the resample brings ALL of its classes -- that is what
  ## keeps the nesting intact
  for (sc in unique(out$school)) {
    expected <- sort(unique(dat$classroom[dat$school == sc]))
    actual  <- sort(unique(out$classroom[out$school == sc]))
    expect_identical(actual, expected, info = paste("school", sc))
  }
  ## and the resample is a resample: whole schools, drawn with replacement
  expect_equal(nrow(out) %% 10L, 0L)
  expect_true(all(out$classroom %in% dat$classroom))
})

test_that("B87: the outer level is found from the data, not from the order", {
  set.seed(872)
  dat <- data.frame(Y = 1:40, school = rep(1:4, each = 10),
                    classroom = rep(1:8, each = 5))
  expect_identical(ci(".countimp_outermost_level")(dat, c("school", "classroom")),
                   "school")
  ## the same answer when the caller lists them the other way round: the
  ## nesting is a property of the data, not of the formula's word order
  expect_identical(ci(".countimp_outermost_level")(dat, c("classroom", "school")),
                   "school")
  ## and the nesting test itself, now taking the two columns rather than the
  ## data frame and two names -- one function where there were two
  expect_true(ci(".countimp_is_nested")(dat$classroom, dat$school))
  expect_false(ci(".countimp_is_nested")(dat$school, dat$classroom))
})

test_that("B87: crossed factors are refused, not silently resampled", {
  ## pupils crossed with teachers: neither contains the other
  dat <- data.frame(Y = 1:12, a = rep(1:3, each = 4), b = rep(1:4, 3))
  expect_true(is.na(ci(".countimp_outermost_level")(dat, c("a", "b"))))
  expect_error(ci(".countimp_boot_clusters")(dat, c("a", "b")), "crossed")
  ## one level still works exactly as before
  out <- ci(".countimp_boot_clusters")(dat, "a")
  expect_equal(nrow(out), 12L)
})

test_that("B87: a three-level bootstrap imputation runs", {
  d <- .b87_daten(nmis = 60L)[, c("y", "x1", "school", "classroom")]
  imp <- suppressWarnings(countimp(d,
           formulas = list(y ~ x1 + (1 | school) + (1 | classroom)),
           family = poisson_count(), draw = "boot", m = 1, maxit = 1,
           printFlag = FALSE))
  v <- unlist(imp$imp$y)
  expect_equal(length(v), sum(is.na(d$y)))
  expect_false(anyNA(v))
  expect_true(all(v >= 0))
})

test_that("B87: 2l.pmm still refuses a second level, and says why", {
  d <- .b87_daten(nmis = 40L)[, c("y", "x1", "school", "classroom")]
  pred <- matrix(0L, 4L, 4L, dimnames = list(names(d), names(d)))
  pred["y", ] <- c(0L, 1L, -2L, -2L)
  expect_error(suppressWarnings(countimp(d, method = c(y = "2l.pmm", x1 = "",
                                                       school = "", classroom = ""),
                                         predictorMatrix = pred, m = 1,
                                         maxit = 1, printFlag = FALSE)),
               "one grouping level")
})


## A model richer than the data support says so -------------------------------
##
## Two columns coded -2 build a genuine three-level model, (1|g) + (1|G), and
## on data that carry that structure it runs. On data that do not, the
## covariance of the fixed effects comes out singular -- and until 1 September
## 2026 what the user saw was chol()'s own words:
##
##   the leading minor of order 1 is not positive
##
## which names the arithmetic and not the situation. The run stopped either
## way; only the message changes. It stops rather than drawing beta at its
## point estimate, because that would understate the parameter uncertainty
## without showing it (B01), and the core words the diagnosis the same way.
##
## Noted on 31 August as "a two-level method cannot carry two groupings" --
## which was wrong, and is recorded as such in ZUSTAND.md. The method carries
## them deliberately; it was the test data that carried nothing.

testthat::test_that("B87: a singular fixed-effect covariance is named, not shown raw", {
  testthat::skip_if_not_installed("glmmTMB")
  set.seed(1)
  n <- 200L
  d <- data.frame(y = stats::rpois(n, 3), x = stats::rnorm(n),
                  z = stats::rnorm(n), k = stats::rnorm(n),
                  g = rep(1:20, each = 10L), G = rep(1:5, each = 40L))
  d$y[sample.int(n, 40L)] <- NA          # plain counts: no cluster structure
  p <- matrix(0L, 6L, 6L, dimnames = list(names(d), names(d)))
  p["y", ] <- c(0L, 2L, 1L, 1L, -2L, -2L)
  err <- tryCatch(suppressWarnings(countimp(d,
           method = c("2l.nb2", "", "", "", "", ""), predictorMatrix = p,
           m = 1L, maxit = 1L, printFlag = FALSE)), error = conditionMessage)
  ## Reaching the message requires chol(vcov(fit)$cond) to fail, and that
  ## requires the fit NOT to converge -- the covariance is singular only
  ## because glmmTMB stops at a non-positive-definite Hessian. Whether it does
  ## is a property of the optimiser and the platform's BLAS, not of the data:
  ## win-builder R-devel (x86_64-w64-mingw32, 2026-09-02) fits these data
  ## cleanly and returns imputations, where macOS arm64 stops. Both are correct
  ## behaviour for countimp. So the branch is tested where it is reached, and
  ## the precondition is stated instead of assumed.
  if (!is.character(err))
    testthat::skip("the fit converged here, so the singular branch is not reached")
  testthat::expect_match(err, "covariance of the fixed effects is singular",
                         fixed = TRUE)
  ## names the situation the user is in -- two grouping levels, and which
  testthat::expect_match(err, "2 grouping levels", fixed = TRUE)
  testthat::expect_match(err, "g, G", fixed = TRUE)
  ## and says why it stops instead of carrying on
  testthat::expect_match(err, "B01", fixed = TRUE)
  ## the raw LAPACK wording is gone
  testthat::expect_no_match(err, "leading minor", fixed = TRUE)
})

testthat::test_that("B87: three levels still run where the data carry them", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("glmmTMB")
  ## The guard against over-correcting: the refusal must depend on the DATA,
  ## not on the number of -2 columns.
  set.seed(8)
  nS <- 15L; nK <- 4L; nP <- 12L
  school <- rep(seq_len(nS), each = nK * nP)
  classroom <- rep(seq_len(nS * nK), each = nP)
  n <- length(school)
  uS <- stats::rnorm(nS, 0, 0.5)[school]
  uK <- stats::rnorm(nS * nK, 0, 0.4)[classroom]
  x <- stats::rnorm(n)
  d <- data.frame(y = MASS::rnegbin(n, mu = exp(0.6 + 0.3 * x + uS + uK),
                                    theta = 3),
                  x = x, classroom = classroom, school = school)
  d$y[sample.int(n, 200L)] <- NA
  p <- matrix(0L, 4L, 4L, dimnames = list(names(d), names(d)))
  p["y", ] <- c(0L, 2L, -2L, -2L)
  testthat::expect_s3_class(suppressWarnings(countimp(d,
    method = c("2l.nb2", "", "", ""), predictorMatrix = p,
    m = 1L, maxit = 1L, printFlag = FALSE)), "mids")
})
