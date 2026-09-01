## B90: variables that live at a cluster level
##
## Two halves, and the first is a behaviour change:
##
##  1. THE REFUSAL. A variable that is constant within a cluster and imputed row
##     by row produces a completed data set the design cannot contain -- a class
##     whose class size differs between its own pupils. Measured before this
##     block existed: 12 of 60 classes lost their constancy, while the
##     regression coefficient looked almost unharmed (0.486 against 0.507 on
##     complete data, true 0.500). That combination is what makes it dangerous:
##     the damage is in the structure, not in the number one would look at.
##  2. THE METHODS. 2lonly.* aggregate to one row per cluster, impute there and
##     hand the value back to every row. The property to defend is exactly the
##     one that was broken: one value per cluster, always.
##
## The trap in building this, worth naming because it is invisible: inside a run
## the engine has already filled every missing cell with a hot-deck starting
## value, so `is.na(y)` finds nothing and a cluster-level variable arrives
## row-wise varying. Everything here must go through `ry`.

.b90_daten <- function(nclass = 40L, nsub = 10L, seed = 90L, kl_na = 8L,
                       zaehlung = FALSE) {
  set.seed(seed)
  classroom <- rep(seq_len(nclass), each = nsub)
  csize <- rep(stats::rnorm(nclass), each = nsub)
  if (zaehlung) csize <- round(exp(csize + 1.5))
  x1 <- stats::rnorm(nclass * nsub)
  y <- stats::rpois(nclass * nsub, exp(0.6 + 0.4 * x1 + 0.5 * scale(csize)[, 1]))
  d <- data.frame(y = y, x1 = x1, csize = csize, classroom = classroom)
  d$csize[d$classroom %in% sample.int(nclass, kl_na)] <- NA
  d
}

.b90_pred <- function() {
  p <- matrix(0L, 4L, 4L, dimnames = list(c("y", "x1", "csize", "classroom"),
                                          c("y", "x1", "csize", "classroom")))
  p["csize", ] <- c(1L, 1L, 0L, -2L)
  p
}

.b90_je_cluster <- function(cd)
  max(tapply(cd$csize, cd$classroom, function(z) length(unique(z))))


## --- 1: the refusal ----------------------------------------------------------

test_that("B90: row-wise imputation of a cluster-level variable is refused", {
  d <- .b90_daten()
  err <- tryCatch(suppressWarnings(countimp(d,
           method = c(y = "", x1 = "", csize = "pmm", classroom = ""),
           predictorMatrix = .b90_pred(), m = 1, maxit = 1, printFlag = FALSE)),
         error = function(e) conditionMessage(e))
  expect_type(err, "character")
  ## the message must name the variable, the level, and the way out
  expect_match(err, "'csize' is constant within 'classroom'", fixed = TRUE)
  expect_match(err, "2lonly.pmm", fixed = TRUE)
  expect_match(err, "countimp.check.levels", fixed = TRUE)
})

test_that("B90: the refusal can be switched off, and then the damage is real", {
  d <- .b90_daten()
  ## base R rather than withr::local_options(): withr is not a declared
  ## dependency of this package, and one test is not worth adding one.
  alt <- options(countimp.check.levels = FALSE)
  on.exit(options(alt), add = TRUE)
  imp <- suppressWarnings(countimp(d,
           method = c(y = "", x1 = "", csize = "pmm", classroom = ""),
           predictorMatrix = .b90_pred(), m = 1, maxit = 1, printFlag = FALSE))
  ## this is what the check exists to prevent: classes with several values
  expect_gt(.b90_je_cluster(countimp_complete(imp, 1L)), 1L)
})

test_that("B90: variables that genuinely vary within clusters are untouched", {
  ## `x1` varies inside every class, so nothing here concerns it
  d <- .b90_daten()
  d$x1[sample.int(nrow(d), 60L)] <- NA
  p <- .b90_pred()
  p["x1", ] <- c(1L, 0L, 1L, -2L)
  p["csize", ] <- 0L
  expect_silent(suppressWarnings(countimp(d,
    method = c(y = "", x1 = "pmm", csize = "", classroom = ""),
    predictorMatrix = p, m = 1, maxit = 1, printFlag = FALSE)))
})

test_that("B90: too few observed clusters means no verdict", {
  ## With one or two clusters, constancy is not evidence of a level. The check
  ## must stay quiet rather than guess.
  d <- .b90_daten(nclass = 2L, nsub = 20L, kl_na = 0L)
  d$csize[1:5] <- NA
  expect_silent(suppressWarnings(countimp(d,
    method = c(y = "", x1 = "", csize = "pmm", classroom = ""),
    predictorMatrix = .b90_pred(), m = 1, maxit = 1, printFlag = FALSE)))
})


## --- 2: the methods ----------------------------------------------------------

test_that("B90: every 2lonly method returns one value per cluster", {
  for (meth in c("2lonly.pmm", "2lonly.norm", "2lonly.pois", "2lonly.nb")) {
    zaehl <- meth %in% c("2lonly.pois", "2lonly.nb")
    d <- .b90_daten(zaehlung = zaehl)
    imp <- suppressWarnings(countimp(d,
             method = c(y = "", x1 = "", csize = meth, classroom = ""),
             predictorMatrix = .b90_pred(), m = 2, maxit = 1, printFlag = FALSE))
    for (i in 1:2) {
      cd <- countimp_complete(imp, i)
      expect_identical(.b90_je_cluster(cd), 1L, info = meth)
      expect_false(anyNA(cd$csize), info = meth)
    }
    ## the count variants must return counts
    if (zaehl) {
      v <- unlist(imp$imp$csize)
      expect_true(all(v >= 0 & v == round(v)), info = meth)
    }
  }
})

test_that("B90: the imputations differ between m, and stay in range", {
  d <- .b90_daten()
  imp <- suppressWarnings(countimp(d,
           method = c(y = "", x1 = "", csize = "2lonly.pmm", classroom = ""),
           predictorMatrix = .b90_pred(), m = 5, maxit = 1, printFlag = FALSE))
  M <- vapply(1:5, function(i) mean(countimp_complete(imp, i)$csize), numeric(1))
  ## proper imputation: not the same value every time
  expect_gt(stats::sd(M), 0)
  ## pmm draws from observed donors, so nothing outside their range
  obs_ <- range(d$csize, na.rm = TRUE)
  expect_true(all(unlist(imp$imp$csize) >= obs_[1L] &
                  unlist(imp$imp$csize) <= obs_[2L]))
})

test_that("B90: a variable that is not cluster-level is refused by the method", {
  ## the mirror image of the engine check: asking for 2lonly on a variable that
  ## varies within clusters
  d <- .b90_daten()
  d$x1[sample.int(nrow(d), 60L)] <- NA
  p <- .b90_pred()
  p["x1", ] <- c(1L, 0L, 1L, -2L)
  p["csize", ] <- 0L
  err <- tryCatch(suppressWarnings(countimp(d,
           method = c(y = "", x1 = "2lonly.pmm", csize = "", classroom = ""),
           predictorMatrix = p, m = 1, maxit = 1, printFlag = FALSE)),
         error = function(e) conditionMessage(e))
  expect_type(err, "character")
  expect_match(err, "not constant within any of the grouping")
})

test_that("B90: the level check reads observed cells only", {
  ## The engine fills missing cells with starting values before the first draw.
  ## A check on is.na() would see a row-wise varying variable and refuse its own
  ## data -- measured while building this: the first call failed with "not
  ## constant within any of the grouping variable(s)".
  kk <- ci(".countimp_constant_within")
  y <- c(1, 1, 99, 2, 2, 98)        # 99/98 are starting values, not data
  g <- c("a", "a", "a", "b", "b", "b")
  ry <- c(TRUE, TRUE, FALSE, TRUE, TRUE, FALSE)
  expect_true(kk(y, g, ry))         # constant where observed
  expect_false(kk(y, g))            # without ry the starting values count
})


## The unnamed method vector -------------------------------------------------
##
## Reported from the simulation side on 30 August 2026, while building study 3:
## every 2l.* method was unreachable through the documented route. The check
## above indexes the method vector BY NAME (method[[v]], v from the column
## names of the predictorMatrix), but mice's own documentation writes the
## vector unnamed -- method = c("pmm", "", "", "") -- and mice accepts it. So
## the guard against a real error aborted every legitimate call with "subscript
## out of bounds" before it ever ran, and the only way out was to switch the
## guard off entirely.
##
## It went unseen because the function returns early when no column is coded
## -2: no single-level method reaches the loop, and this project's own studies
## build their method vector through a helper that names it.

test_that("B90: an unnamed method vector reaches the check and imputes", {
  d <- .b90_daten()
  d$y[sample.int(nrow(d), 60L)] <- NA
  p <- .b90_pred()
  p["y", ] <- c(0L, 1L, 1L, -2L)

  ## Both variables are imputed, which puts the two method[[v]] lookups of the
  ## check on the same call: the level test for y, and the "^2lonly" skip for
  ## csize. Unnamed -- the form mice documents -- and named must agree
  ## completely, not merely both run: the vector is read positionally, so a
  ## wrong mapping would silently impute with the wrong method rather than
  ## fail.
  unben <- c("2l.nb2", "", "2lonly.pmm", "")
  ben <- c(y = "2l.nb2", x1 = "", csize = "2lonly.pmm", classroom = "")
  set.seed(4L)
  unbenannt <- suppressWarnings(countimp(d, method = unben,
                 predictorMatrix = p, m = 1L, maxit = 1L, printFlag = FALSE))
  set.seed(4L)
  benannt <- suppressWarnings(countimp(d, method = ben, predictorMatrix = p,
               m = 1L, maxit = 1L, printFlag = FALSE))
  expect_identical(unbenannt$imp$y, benannt$imp$y)
  expect_identical(unbenannt$imp$csize, benannt$imp$csize)

  ## and the imputation is usable: every missing y filled with a count
  yv <- unbenannt$imp$y[[1L]]
  expect_equal(length(yv), 60L)
  expect_false(anyNA(yv))
  expect_true(all(yv >= 0 & yv == round(yv)))
})

test_that("B90: a named partial method vector selects the right variables", {
  ## mice allows method = c(y = "2l.nb2") alone. Before the fix the positional
  ## read (vn[nzchar(method)]) picked position 1 -- which happens to be y here,
  ## so this case needs a variable that is NOT first to have any force.
  kk <- ci(".countimp_check_levels")
  d <- .b90_daten()
  p <- .b90_pred()

  ## csize is cluster-level and sits at position 3; a partial vector naming it
  ## must reach the refusal, not check y by mistake
  err <- tryCatch(kk(d, c(csize = "2l.pmm"), p), error = conditionMessage)
  expect_type(err, "character")
  expect_match(err, "is constant within")

  ## and a partial vector naming a different variable must NOT refuse
  expect_null(kk(d, c(x1 = "2l.pmm"), p))
})

test_that("B90: the level check refuses an unnamed vector of the wrong length", {
  ## Defensive, and deliberately tested at function level: through countimp()
  ## mice rejects the mismatch first ("The length of method (2) does not match
  ## the number of columns in the data (5)"), so this branch is reachable only
  ## by calling the check directly -- which is what the 2l draw functions do.
  kk <- ci(".countimp_check_levels")
  d <- .b90_daten()
  p <- .b90_pred()
  err <- tryCatch(kk(d, c("2l.nb2", ""), p), error = conditionMessage)
  expect_match(err, "read positionally")
  expect_no_match(err, "subscript out of bounds")
})


## Too FINE a level is refused as well ---------------------------------------
##
## The mirror image of the block at the top of this file, and the more
## treacherous of the two. A school budget is constant within its school and
## therefore within every class of that school. Aggregating to the class level
## passes every check the aggregation itself can make -- one value per class,
## all consistent -- and then draws one value PER CLASS. The completed data
## then holds a school with several budgets.
##
## Measured before this check existed, 12 schools with 4 classes each: 4 of 12
## schools came out with more than one budget. No error, no warning.
##
## The core has refused it since 25 August (coarser_level_available()); this is
## the R side catching up. It sits in .countimp_check_levels() and not in the
## method, because a mice.impute.* function only sees what the predictorMatrix
## lets through -- with the school coded 0 it is not in `x` at all, and the
## check would be blind to exactly the level it must warn about.

.b90_drei <- function(nS = 12L, nK = 4L, nP = 15L, seed = 3L) {
  set.seed(seed)
  school <- rep(seq_len(nS), each = nK * nP)
  classroom <- rep(seq_len(nS * nK), each = nP)
  n <- length(school)
  budget <- rep(round(exp(stats::rnorm(nS, 2.5, 0.4))), each = nK * nP)
  x1 <- stats::rnorm(n)
  d <- data.frame(y = stats::rpois(n, exp(0.4 + 0.3 * x1)), x1 = x1,
                  budget = budget, classroom = classroom, school = school)
  d$budget[d$school %in% c(2L, 5L, 9L, 11L)] <- NA
  d
}

test_that("B90: a 2lonly method pointed at too fine a level is refused", {
  d <- .b90_drei()
  p <- matrix(0L, 5L, 5L, dimnames = list(names(d), names(d)))
  p["budget", ] <- c(1L, 1L, 0L, -2L, 0L)      # class as grouping, school not
  err <- tryCatch(suppressWarnings(countimp(d,
           method = c(budget = "2lonly.pois"), predictorMatrix = p,
           m = 1L, maxit = 1L, printFlag = FALSE)), error = conditionMessage)
  expect_type(err, "character")
  expect_match(err, "COARSER level", fixed = TRUE)
  expect_match(err, "'school'", fixed = TRUE)   # names the level to use
  ## and the school is not even a predictor here -- the check must see the
  ## whole data frame, not just the columns the method receives
  expect_equal(unname(p["budget", "school"]), 0L)
})

test_that("B90: grouping by the right level still works", {
  testthat::skip_on_cran()
  d <- .b90_drei()
  p <- matrix(0L, 5L, 5L, dimnames = list(names(d), names(d)))
  p["budget", ] <- c(1L, 1L, 0L, -2L, -2L)     # BOTH levels offered
  set.seed(11)
  imp <- suppressWarnings(countimp(d, method = c(budget = "2lonly.pois"),
           predictorMatrix = p, m = 2L, maxit = 1L, printFlag = FALSE))
  cd <- mice::complete(imp, 1L)
  ## the outermost level wins, so every school keeps ONE budget
  je_schule <- tapply(cd$budget, cd$school, function(z) length(unique(z)))
  expect_equal(max(je_schule), 1L)
})

test_that("B90: the two conditions do not cover for each other", {
  ## The trap the core's own mutation probe fell into: "nested" and "strictly
  ## coarser" can both be false at once, and then removing either changes
  ## nothing. These two cases separate them.
  vn <- ci(".countimp_is_nested")

  ## (a) a mere RENAMING of the same grouping: nested, but not coarser
  classroom <- rep(1:8, each = 5L)
  neu <- classroom + 100L
  expect_true(vn(classroom, neu))                      # nested
  expect_equal(length(unique(neu)), length(unique(classroom)))   # not coarser

  ## (b) a coarser column that lies ACROSS the grouping: coarser, not nested
  quer <- rep(1:2, times = 20L)
  expect_lt(length(unique(quer)), length(unique(classroom)))     # coarser
  expect_false(vn(classroom, quer))                              # not nested
})

## Each of the three conditions alone -----------------------------------------
##
## The refusal rests on three: the candidate must be NESTED around the chosen
## grouping, STRICTLY COARSER than it, and the variable must be constant within
## it too. A test where one candidate satisfies all three cannot tell them
## apart -- drop any one condition and the other two still refuse. A mutation
## probe said exactly that: removing the nesting test, the size test or the
## constancy test changed nothing.
##
## So each case below violates exactly ONE condition, and each must therefore
## RUN. That is what makes dropping that condition visible: the refusal would
## then fire where it must not.
##
## This is the same trap the core fell into on 25 August, recorded in
## ZUSTAND.md 3p -- and this time it was found by probing rather than by luck.

.b90_eine_ebene <- function(nK = 20L, nP = 12L, seed = 5L) {
  set.seed(seed)
  classroom <- rep(seq_len(nK), each = nP)
  n <- length(classroom)
  csize <- rep(round(exp(stats::rnorm(nK, 2.2, 0.4))), each = nP)
  x1 <- stats::rnorm(n)
  d <- data.frame(y = stats::rpois(n, exp(0.4 + 0.3 * x1)), x1 = x1,
                  csize = csize, classroom = classroom)
  d$csize[d$classroom %in% c(2L, 7L, 13L)] <- NA
  d
}

.b90_lauf <- function(d, pred) {
  set.seed(11)
  suppressWarnings(countimp(d, method = c(csize = "2lonly.pois"),
    predictorMatrix = pred, m = 1L, maxit = 1L, printFlag = FALSE))
}

test_that("B90: a mere RENAMING of the grouping does not trigger the refusal", {
  ## Violates SIZE only: nested, constant within, but exactly as many distinct
  ## values. Without the size test this would be refused, and then every data
  ## set carrying a second identifier column would be.
  d <- .b90_eine_ebene()
  d$klasse_id <- d$classroom + 1000L
  p <- matrix(0L, 5L, 5L, dimnames = list(names(d), names(d)))
  p["csize", ] <- c(1L, 1L, 0L, -2L, 0L)
  expect_s3_class(.b90_lauf(d, p), "mids")
})

test_that("B90: a coarser column ACROSS the grouping does not trigger it", {
  ## Violates NESTING only -- and getting there needs care, which is itself
  ## worth recording.
  ##
  ## In the ordinary case the nesting test cannot be separated from the
  ## constancy test: if a column lies ACROSS the classes, then each of its
  ## levels contains rows from several classes, and a class-level variable is
  ## therefore not constant within it. Constancy already refuses. Measured --
  ## a mutation removing the nesting test changed nothing.
  ##
  ## The case where they come apart is the degenerate one: a variable whose
  ## observed values happen to be all alike. Then it is constant within EVERY
  ## column, including one that lies across the grouping, and only the nesting
  ## test still stands between the run and a wrong refusal.
  ##
  ## HONEST LIMIT: this block does not separate them either. Calling
  ## .countimp_check_levels() directly on such data does -- measured, `quer`
  ## comes out coarser and constant, so without the nesting test it refuses --
  ## but through countimp() the mutation still passes. Whatever swallows it
  ## sits between the two, and it has not been chased down. So: the nesting
  ## test is NOT covered by a test that would fail without it. It stays in as
  ## a cheap pre-filter and because the core carries the same condition
  ## (is_nested_in), not because anything here proves it earns its place.
  d <- .b90_eine_ebene()
  d$csize <- ifelse(is.na(d$csize), NA_real_, 9)   # every observed value alike
  d$quer <- rep(1:2, length.out = nrow(d))     # alternates WITHIN each class
  p <- matrix(0L, 5L, 5L, dimnames = list(names(d), names(d)))
  p["csize", ] <- c(1L, 1L, 0L, -2L, 0L)
  expect_lt(length(unique(d$quer)), length(unique(d$classroom)))
  ## constancy alone cannot refuse here: csize IS constant within quer
  expect_equal(length(unique(d$csize[!is.na(d$csize)])), 1L)
  expect_s3_class(.b90_lauf(d, p), "mids")
})

test_that("B90: a coarser level the variable VARIES within does not trigger it", {
  ## Violates CONSTANCY only: a genuine level above the class, but the class
  ## size differs between the classes of one school -- so the school is not
  ## the level of this variable and grouping by the class is right.
  d <- .b90_eine_ebene()
  d$school <- ((d$classroom - 1L) %/% 4L) + 1L    # 5 schools of 4 classes
  p <- matrix(0L, 5L, 5L, dimnames = list(names(d), names(d)))
  p["csize", ] <- c(1L, 1L, 0L, -2L, 0L)
  ## the variable is genuinely NOT constant within school
  const_ <- tapply(d$csize, d$school, function(z) length(unique(z[!is.na(z)])))
  expect_gt(max(const_), 1L)
  expect_s3_class(.b90_lauf(d, p), "mids")
})
