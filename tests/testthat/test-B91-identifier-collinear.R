## B91: a grouping identifier must not be dropped for collinearity
##
## Found while building the three-level simulations, reported from there.
##
## The mechanism. In a school-class data set both identifiers are numbered
## consecutively, so they correlate -- classes 41 to 60 simply sit in the
## higher-numbered schools. Measured: 0.998827 at 20 schools, 0.999479 at 30.
## The default threshold of find.collinear() is 0.999. From 30 schools upwards
## the class identifier therefore vanished from the predictorMatrix before any
## method got to see it.
##
## Both consequences were silent and both were wrong:
##   * a model meant as three-level ran on as two-level, without a warning --
##     the variance component of the inner level is then simply missing;
##   * a 2lonly method saw type = (1, 1) and asked the user to "code one
##     predictor as -2" -- an instruction they had already followed.
##
## The second case is the more instructive one: the error message was correct
## about the state the method saw, and misleading nonetheless, because that
## state was not the one the user had created.
##
## The rule against it fits on one line: an identifier is a label, not a
## covariate. Its correlation with another identifier is a property of the
## numbering, not of the model.

skip_if_not_installed("glmmTMB")

## Exactly the setup where it tipped over: the number of schools decides.
.b91_daten <- function(S = 30L, K = 4L, nsub = 20L, seed = 91L) {
  set.seed(seed)
  school <- rep(seq_len(S), each = K * nsub)
  classroom <- rep(seq_len(S * K), each = nsub)
  n <- length(school)
  csize <- rep(stats::rnorm(S * K), each = nsub)
  x1 <- stats::rnorm(n)
  y <- stats::rpois(n, exp(0.6 + 0.4 * x1 + 0.5 * csize))
  d <- data.frame(y = y, x1 = x1, csize = csize, school = school, classroom = classroom)
  d$y[sample.int(n, round(0.2 * n))] <- NA
  d$csize[d$classroom %in% sample(unique(classroom), round(0.2 * S * K))] <- NA
  d
}

.b91_pred <- function() {
  p <- matrix(0L, 5L, 5L,
              dimnames = list(c("y", "x1", "csize", "school", "classroom"),
                              c("y", "x1", "csize", "school", "classroom")))
  p["y", ]   <- c(0L, 1L, 1L, -2L, -2L)     # three-level count model
  p["csize", ] <- c(0L, 1L, 0L, -2L, -2L)     # class-level variable
  p
}


test_that("B91: the identifiers correlate strongly enough to cross the threshold", {
  ## Without this check the test below would say nothing: it would confirm a
  ## rule whose trigger never occurs.
  for (S in c(20L, 30L)) {
    school <- rep(seq_len(S), each = 4L * 20L)
    classroom <- rep(seq_len(S * 4L), each = 20L)
    r <- stats::cor(school, classroom)
    if (S == 20L) expect_lt(r, 0.999)       # just below
    if (S == 30L) expect_gt(r, 0.999)       # just above -- this is where it tipped
  }
})

test_that("B91: the class identifier survives the collinearity check", {
  d <- .b91_daten(S = 30L)
  meth <- c(y = "2l.poisson", x1 = "", csize = "2lonly.pmm",
            school = "", classroom = "")
  imp <- suppressWarnings(countimp(d, method = meth,
           predictorMatrix = .b91_pred(), m = 1, maxit = 1, printFlag = FALSE))
  ## both variables were imputed
  expect_false(anyNA(unlist(imp$imp$y)))
  expect_false(anyNA(unlist(imp$imp$csize)))
  ## and the cluster variable stayed constant within each class
  cd <- countimp_complete(imp, 1L)
  expect_identical(max(tapply(cd$csize, cd$classroom, function(z) length(unique(z)))),
                   1L)
})

test_that("B91: the same at 20 schools, where the threshold did not bite", {
  ## The case that worked before. It has to keep working -- a fix must not
  ## repair one half and break the other.
  d <- .b91_daten(S = 20L)
  meth <- c(y = "2l.poisson", x1 = "", csize = "2lonly.pmm",
            school = "", classroom = "")
  imp <- suppressWarnings(countimp(d, method = meth,
           predictorMatrix = .b91_pred(), m = 1, maxit = 1, printFlag = FALSE))
  expect_false(anyNA(unlist(imp$imp$csize)))
})

test_that("B91: the method sees both identifiers, not just one", {
  ## The property at issue, checked directly: the type vector arriving at the
  ## method carries both -2 entries. A spy on the method is the right tool
  ## here, because the imputed values do not reveal how many levels were
  ## used to produce them.
  seen <- NULL
  ci_spy("mice.impute.2l.poisson", function(y, ry, x, type, ...) {
    seen <<- type
    rep(1, sum(!ry))
  })
  d <- .b91_daten(S = 30L)
  suppressWarnings(countimp(d,
    method = c(y = "2l.poisson", x1 = "", csize = "", school = "", classroom = ""),
    predictorMatrix = .b91_pred(), m = 1, maxit = 1, printFlag = FALSE))
  expect_false(is.null(seen))
  expect_identical(sum(seen == -2L), 2L)
})

test_that("B91: ordinary predictors are still checked for collinearity", {
  ## The exemption covers identifiers only. Two identical covariates must
  ## still be caught, or the fix would have disabled the check.
  set.seed(912); n <- 300L
  x1 <- stats::rnorm(n)
  d <- data.frame(y = stats::rpois(n, exp(0.5 + 0.4 * x1)),
                  x1 = x1, x2 = x1,          # exactly the same column
                  id = rep(seq_len(30), each = 10L))
  d$y[sample.int(n, 60L)] <- NA
  p <- matrix(0L, 4L, 4L, dimnames = list(names(d), names(d)))
  p["y", ] <- c(0L, 1L, 1L, -2L)
  ## countimp reports the dropped column through its log; here it is enough
  ## that the run completes and does not fail on the singularity
  imp <- suppressWarnings(countimp(d, method = c(y = "2l.poisson", x1 = "",
                                                 x2 = "", id = ""),
           predictorMatrix = p, m = 1, maxit = 1, printFlag = FALSE))
  expect_false(anyNA(unlist(imp$imp$y)))
})


## --- B92: the same defect in the second filter ------------------------------
##
## check.data() is not the only place that discards columns: remove.lindep()
## does it again inside the sampler, per variable and per iteration, and it
## works on the OBSERVED rows through an eigendecomposition. That is why it
## depended on sample and seed which of the two identifiers fell -- measured on
## the same setup, only the number of schools differing:
##
##   S = 50   drops 'school'  -> 'class' remains  -> later error
##   S = 60   drops 'class'   -> 'school' remains -> RUNS, silently two-level
##
## A call meant as three-level became two-level, without a word. Reported from
## the simulation side as a "sporadic 2l.pmm case".

test_that("B92: remove.lindep drops no grouping identifier", {
  ## Kept small: the property is what arrives at the method, not the result of
  ## an expensive multilevel fit.
  ## The sizes are measured, not chosen: of seven setups tried, this is the
  ## smallest where remove.lindep() actually discards an identifier
  ## (n = 1200). At S = 25, K = 4 it discards nothing and the test would have
  ## said nothing -- the mutation probe showed that before this version
  ## existed.
  ##
  ##   S=40 K=4  cor 0.999707  ->  nothing
  ##   S=50 K=4  cor 0.999812  ->  'class' falls
  ##   S=60 K=2  cor 0.999896  ->  'class' falls        <- used here
  set.seed(921)
  S <- 60L; K <- 2L; nsub <- 10L
  school <- rep(seq_len(S), each = K * nsub)
  classroom <- rep(seq_len(S * K), each = nsub)
  n <- length(school)
  d <- data.frame(y = stats::rpois(n, 3), x1 = stats::rnorm(n),
                  school = school, classroom = classroom)
  d$y[sample.int(n, round(0.2 * n))] <- NA
  expect_gt(stats::cor(school, classroom), 0.9998)

  p <- matrix(0L, 4L, 4L, dimnames = list(names(d), names(d)))
  p["y", ] <- c(0L, 1L, -2L, -2L)

  seen <- NULL
  ci_spy("mice.impute.2l.poisson", function(y, ry, x, type, ...) {
    seen <<- type
    rep(1, sum(!ry))
  })
  suppressWarnings(countimp(d, method = c(y = "2l.poisson", x1 = "",
                                          school = "", classroom = ""),
    predictorMatrix = p, m = 1, maxit = 1, printFlag = FALSE))
  expect_false(is.null(seen))
  expect_identical(sum(seen == -2L), 2L)
})
