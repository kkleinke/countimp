## B50 -- mice.impute.2l.pmm: the one export no test touched
##
## The export audit mapped all 63 mice.impute.* names onto their engines:
## 16 route to .countimp_2l_hurdle, 16 to .countimp_2l_zi, 12 to
## .countimp_2l_count, 8 to .countimp_1l_zi, 6 to .countimp_1l_count, 4 are
## thin aliases -- and exactly one carries its own model code: 2l.pmm. It was
## also the only export not named by any test file, so every property below was
## unverified until now.
##
## Two defects were found and fixed here, both invisible to the family tests
## because 2l.pmm shares no code with the consolidated engines:
##
##   B51  glmmPQL(verbose = TRUE) is the default. Each call cat()s "Iteration i"
##        for ten PQL steps; inside countimp() with m = 5, maxit = 10 and three
##        incomplete variables that is ~300 lines of noise that bury real
##        warnings. Measured: 2 lines per call before the fix, 0 after.
##
## Also covered here for the first time: the unobserved-cluster path shared with
## B31/B34. Clusters with no observed y get u_j ~ N(0, tau^2) from
## .countimp_predict_2l() rather than the NA that glmmPQL's own predict()
## returns and that would later abort .pmm.match() with "'partial' index outside
## bounds". That fix was made for the consolidated engines and 2l.pmm inherited
## it, but no test exercised it through this method.
##
## The donor mechanism is what distinguishes 2l.pmm from every other method in
## the package: imputations must be drawn from the OBSERVED y values, so the
## strongest available check is set membership, which no parametric method can
## satisfy by accident.
if (!exists("sim_2l_gauss", inherits = FALSE)) {
  sim_2l_gauss <- function(N = 240L, J = 24L, seed = 1L, p = 2L, fehlend = 0.25) {
    set.seed(seed)
    grp <- rep(seq_len(J), each = N %/% J)
    X <- matrix(rnorm(N * p), N, p, dimnames = list(NULL, paste0("x", seq_len(p))))
    u <- rep(rnorm(J, 0, 0.5), each = N %/% J)
    y <- 2 + X %*% rep(0.5, p) + u + rnorm(N, 0, 0.6)
    ry <- rbinom(N, 1L, 1 - fehlend) == 0L
    list(y = as.vector(y), ry = !ry,
         x = data.frame(grp = grp, X),
         type = c(-2L, rep(1L, p)))
  }
}

testthat::test_that("B50: 2l.pmm draws only from observed y values", {
  d <- sim_2l_gauss(seed = 21L)
  im <- suppressWarnings(mice.impute.2l.pmm(y = d$y, ry = d$ry, x = d$x, type = d$type))
  testthat::expect_length(im, sum(!d$ry))
  testthat::expect_true(all(is.finite(im)))
  ## the defining property of predictive mean matching: every imputation is a
  ## real observed value, never an interpolation or a model prediction
  testthat::expect_true(all(im %in% d$y[d$ry]))
  ## and the draw must vary -- a single donor for all cases would also satisfy
  ## set membership
  testthat::expect_gt(length(unique(im)), 2L)
})

testthat::test_that("B50: 2l.pmm is silent (B51)", {
  d <- sim_2l_gauss(seed = 22L)
  ## glmmPQL emits its progress through message(gettextf("iteration %d", i)) --
  ## lower case, on stderr, and translated under a non-English locale. Matching
  ## the literal English string would make this test pass on a German machine
  ## whatever verbose is set to, so collect the CONDITIONS instead of the text.
  nachrichten <- character(0)
  withCallingHandlers(
    suppressWarnings(mice.impute.2l.pmm(y = d$y, ry = d$ry, x = d$x, type = d$type)),
    message = function(m) {
      nachrichten <<- c(nachrichten, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  testthat::expect_length(nachrichten, 0L)
  ## and nothing on stdout either (cat() would bypass the handler above)
  aus <- utils::capture.output(
    invisible(suppressWarnings(suppressMessages(
      mice.impute.2l.pmm(y = d$y, ry = d$ry, x = d$x, type = d$type)))))
  testthat::expect_length(aus, 0L)
})

testthat::test_that("B50: a cluster with no observed y is imputed, not aborted (B31/B34)", {
  d <- sim_2l_gauss(seed = 23L)
  ## blank out one entire cluster -- glmmPQL cannot know its random effect
  leer <- d$x$grp == 3L
  d$ry[leer] <- FALSE
  testthat::expect_true(sum(d$ry[leer]) == 0L)
  im <- suppressWarnings(mice.impute.2l.pmm(y = d$y, ry = d$ry, x = d$x, type = d$type))
  testthat::expect_length(im, sum(!d$ry))
  testthat::expect_true(all(is.finite(im)), info = "NA aus unbeobachtetem Cluster")
  testthat::expect_true(all(im %in% d$y[d$ry]))
})

testthat::test_that("B50: the class variable is required, and pmm takes one", {
  d <- sim_2l_gauss(seed = 24L)
  t2 <- d$type; t2[2] <- -2L   # two class variables
  ## Two -2 codes are a three-level model since B87, and the model-based
  ## two-level methods fit it. pmm does NOT: it goes through glmmPQL (nlme),
  ## which takes a single grouping factor. The refusal stays; what changed is
  ## that the message now says which backend imposes the limit instead of
  ## "only one class allowed!".
  err <- tryCatch(suppressWarnings(
    mice.impute.2l.pmm(y = d$y, ry = d$ry, x = d$x, type = t2)),
    error = function(e) conditionMessage(e))
  testthat::expect_type(err, "character")
  testthat::expect_match(err, "one grouping level", fixed = TRUE)
  testthat::expect_match(err, "glmmPQL", fixed = TRUE)
})

testthat::test_that("B50: an empty random part is explained, not passed to nlme (B53)", {
  d <- sim_2l_gauss(seed = 25L)
  ## intercept = FALSE with no predictor coded 2 leaves nothing random. The old
  ## inline construction handed nlme "~0 | grp"; nlme aborted inside
  ## model.matrix.reStruct() with "'data' must be of a vector type, was 'NULL'".
  testthat::expect_false(any(d$type == 2L))
  testthat::expect_error(
    suppressWarnings(mice.impute.2l.pmm(y = d$y, ry = d$ry, x = d$x,
                                        type = d$type, intercept = FALSE)),
    "Random part of the .* model is empty")
  ## the message must name the way out, not just the failure
  msg <- tryCatch(mice.impute.2l.pmm(y = d$y, ry = d$ry, x = d$x,
                                     type = d$type, intercept = FALSE),
                  error = conditionMessage)
  testthat::expect_match(msg, "code a predictor as 2", fixed = TRUE)
  testthat::expect_false(grepl("vector type", msg, fixed = TRUE))
})

testthat::test_that("B50: intercept = FALSE works when a random slope is coded", {
  testthat::skip_on_cran()
  d <- sim_2l_gauss(seed = 25L)
  t2 <- d$type; t2[2] <- 2L        # x1 becomes fixed + random slope
  set.seed(5L)
  a <- suppressWarnings(mice.impute.2l.pmm(y = d$y, ry = d$ry, x = d$x,
                                           type = t2, intercept = TRUE))
  set.seed(5L)
  b <- suppressWarnings(mice.impute.2l.pmm(y = d$y, ry = d$ry, x = d$x,
                                           type = t2, intercept = FALSE))
  ## same seed, genuinely different random structure: draws must not coincide
  testthat::expect_length(a, sum(!d$ry))
  testthat::expect_length(b, sum(!d$ry))
  testthat::expect_true(all(c(a, b) %in% d$y[d$ry]))
  testthat::expect_false(isTRUE(all.equal(a, b)),
    info = "intercept = FALSE does not change the model")
})

testthat::test_that("B50: a missing class variable is named", {
  d <- sim_2l_gauss(seed = 25L)
  t3 <- d$type; t3[1] <- 1L        # no -2 anywhere
  testthat::expect_error(
    suppressWarnings(mice.impute.2l.pmm(y = d$y, ry = d$ry, x = d$x, type = t3)),
    "class variable", fixed = TRUE)
})

testthat::test_that("B50: donors controls the size of the donor pool", {
  testthat::skip_on_cran()
  d <- sim_2l_gauss(seed = 26L, N = 300L, J = 30L)
  ## with a single donor the imputation is the nearest neighbour and therefore
  ## reproducible without any random draw beyond the beta draw
  set.seed(7L); e1 <- suppressWarnings(mice.impute.2l.pmm(
    y = d$y, ry = d$ry, x = d$x, type = d$type, donors = 1L))
  set.seed(7L); e2 <- suppressWarnings(mice.impute.2l.pmm(
    y = d$y, ry = d$ry, x = d$x, type = d$type, donors = 1L))
  testthat::expect_equal(e1, e2)
  ## a wider pool must still stay inside the observed values
  set.seed(7L); e9 <- suppressWarnings(mice.impute.2l.pmm(
    y = d$y, ry = d$ry, x = d$x, type = d$type, donors = 9L))
  testthat::expect_true(all(e9 %in% d$y[d$ry]))
  testthat::expect_false(isTRUE(all.equal(e1, e9)),
    info = "donors hat keinen Effekt")
})
