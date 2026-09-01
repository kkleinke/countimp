## B43 -- the encapsulation layer for foreign packages.
##
## Two failure modes this guards against:
##
##  1. An optional package is not installed and the user gets a message about
##     their data instead of about the package. Before the fix, calling
##     mice.impute.zip() without pscl reported "failed to fit the imputation
##     model ... The original message was: there is no package called 'pscl'"
##     -- the cause was at the end of a message whose beginning was wrong.
##
##  2. A foreign package renames a field countimp reads. .countimp_feld() and
##     .countimp_theta() turn that into a message naming the class, the fields
##     looked for, and the package version.

testthat::test_that("B43: a missing package is reported as a missing package", {
  err <- tryCatch(ci(".countimp_need_pkg")("einPaketDasEsNichtGibt"),
                  error = function(e) conditionMessage(e))
  testthat::expect_match(err, "is required for")
  testthat::expect_match(err, "not installed")
  testthat::expect_match(err, "install.packages")
  ## an installed package passes silently
  testthat::expect_true(ci(".countimp_need_pkg")("stats"))
})

testthat::test_that("B43: the message names the alternative when there is one", {
  err <- tryCatch(ci(".countimp_need_pkg")("keinPaket", what = "method X",
                                           ersatz = "method Y"),
                  error = function(e) conditionMessage(e))
  testthat::expect_match(err, "method X", fixed = TRUE)
  testthat::expect_match(err, "or use method Y", fixed = TRUE)
})

testthat::test_that("B43: the pscl guard sits before the pscl call", {
  ## pscl is in Suggests, not in Imports: the one-level zero-inflation and
  ## hurdle methods need it, the rest of the package does not. Without this
  ## skip the test fails in an environment without Suggests -- and it fails on
  ## the missing dependency, not on the property under test.
  skip_if_not_installed("pscl")
  ## The guard must be the first statement of .countimp_zi_fit, so that a
  ## missing package cannot be reported as a fitting failure.
  src <- paste(deparse(body(ci(".countimp_zi_fit"))), collapse = "\n")
  pos_wache <- regexpr(".countimp_need_pkg", src, fixed = TRUE)
  pos_pscl  <- regexpr("pscl::", src, fixed = TRUE)
  testthat::expect_gt(pos_wache, 0)
  testthat::expect_gt(pos_pscl, 0)
  testthat::expect_lt(pos_wache, pos_pscl)
})

testthat::test_that("B43: .countimp_have_pkg reports availability, not errors", {
  testthat::expect_true(ci(".countimp_have_pkg")("stats"))
  testthat::expect_false(ci(".countimp_have_pkg")("einPaketDasEsNichtGibt"))
})

testthat::test_that("B43: theta is read from all three foreign conventions", {
  set.seed(3); n <- 300
  x1 <- stats::rnorm(n)
  y <- MASS::rnegbin(n, mu = exp(1 + 0.5 * x1), theta = 1.4)
  d <- data.frame(y = y, x1 = x1)

  ## MASS::glm.nb -- theta on the fit and in the summary
  f1 <- MASS::glm.nb(y ~ x1, data = d)
  s1 <- summary(f1)
  testthat::expect_equal(ci(".countimp_theta")(f1, s1), s1$theta)
  testthat::expect_equal(ci(".countimp_theta")(f1), f1$theta)

  testthat::skip_if_not_installed("pscl")
  yz <- y; yz[stats::runif(n) < 0.35] <- 0
  dz <- data.frame(y = yz, x1 = x1)
  f2 <- pscl::zeroinfl(y ~ x1 | x1, data = dz, dist = "negbin")
  testthat::expect_equal(as.numeric(ci(".countimp_theta")(f2)), as.numeric(f2$theta))
  f3 <- pscl::hurdle(y ~ x1 | x1, data = dz, dist = "negbin")
  testthat::expect_equal(as.numeric(ci(".countimp_theta")(f3))[1],
                         as.numeric(f3$theta)[1])
})

testthat::test_that("B43: a renamed theta slot gives a message naming the package", {
  testthat::skip_if_not_installed("pscl")
  set.seed(4); n <- 200
  x1 <- stats::rnorm(n)
  y <- MASS::rnegbin(n, mu = exp(1 + 0.4 * x1), theta = 1.6)
  y[stats::runif(n) < 0.3] <- 0
  f <- pscl::zeroinfl(y ~ x1 | x1, data = data.frame(y = y, x1 = x1),
                      dist = "negbin")
  ## simulate the rename pscl could plausibly make
  f$theta <- NULL
  err <- tryCatch(ci(".countimp_theta")(f), error = function(e) conditionMessage(e))
  testthat::expect_match(err, "cannot read theta")
  testthat::expect_match(err, "pscl")
  testthat::expect_match(err, "\\$theta, \\$logtheta")

  ## the documented fallback: a logtheta slot is honoured
  f$logtheta <- log(2.5)
  testthat::expect_equal(ci(".countimp_theta")(f), 2.5)
})

testthat::test_that("B43: sigma() is used only where it IS a dispersion", {
  ## For an lm, stats::sigma() returns the residual standard deviation. Reading
  ## that as a dispersion would make .countimp_check_fit() flag a healthy fit,
  ## so the accessor must refuse the class rather than return a number.
  set.seed(9); n <- 200
  x <- stats::rnorm(n)
  d <- data.frame(y = MASS::rnegbin(n, mu = exp(1 + 0.4 * x), theta = 2), x = x)
  f_lm <- stats::lm(y ~ x, data = d)
  testthat::expect_null(ci(".countimp_theta")(f_lm, pflicht = FALSE))
  err <- tryCatch(ci(".countimp_theta")(f_lm), error = function(e) conditionMessage(e))
  testthat::expect_match(err, "cannot read a dispersion parameter")
  testthat::expect_match(err, "class lm")
  ## and the diagnostics stay silent on it
  testthat::expect_length(ci(".countimp_check_fit")(f_lm, method = "lm"), 0L)
})

testthat::test_that("B43: degenerate dispersion is still detected after the rewrite", {
  set.seed(9); n <- 200
  x <- stats::rnorm(n)
  d <- data.frame(y = MASS::rnegbin(n, mu = exp(1 + 0.4 * x), theta = 2), x = x)
  f <- MASS::glm.nb(y ~ x, data = d)
  g <- f; g$theta <- 2e6
  testthat::expect_true("theta_degenerate_high" %in%
                          ci(".countimp_check_fit")(g, method = "nb"))
  g$theta <- 1e-8
  testthat::expect_true("theta_degenerate_low" %in%
                          ci(".countimp_check_fit")(g, method = "nb"))
  g$theta <- 2
  testthat::expect_length(ci(".countimp_check_fit")(g, method = "nb"), 0L)
})

testthat::test_that("B43: glmmTMB dispersion goes through the same accessor", {
  testthat::skip_on_cran()
  set.seed(9); G <- 25; nj <- 8; N <- G * nj
  grp <- rep(seq_len(G), each = nj)
  u <- stats::rnorm(G, 0, 0.4); x <- stats::rnorm(N)
  d <- data.frame(y = MASS::rnegbin(N, mu = exp(1 + 0.4 * x + u[grp]), theta = 2),
                  x = x, g = factor(grp))
  f <- glmmTMB::glmmTMB(y ~ x + (1 | g), family = glmmTMB::nbinom2, data = d)
  testthat::expect_equal(ci(".countimp_theta")(f), glmmTMB::sigma(f))
  testthat::expect_length(ci(".countimp_check_fit")(f, method = "2l.nb"), 0L)
})

testthat::test_that("B43: call sites keep a length-one theta", {
  ## The accessor returns NULL when no convention applies, and as.numeric(NULL)
  ## is numeric(0) -- which makes `if (is.finite(th) && th > 0)` raise
  ## "argument is of length zero" instead of taking the fallback branch. Every
  ## call site must therefore coalesce to NA_real_. Guarded by inspection of
  ## the sources, since the failure only shows on an unusual fit object.
  pfad <- if (dir.exists("../../R")) "../../R" else "R"
  testthat::skip_if_not(dir.exists(pfad), "package sources not locatable")
  for (f in list.files(pfad, pattern = "[.]R$", full.names = TRUE)) {
    if (basename(f) == "foreign-packages.R") next
    zeilen <- sub("#.*$", "", readLines(f, warn = FALSE))
    idx <- grep("as\\.numeric\\(\\.countimp_theta\\(", zeilen)
    for (i in idx)
      testthat::expect_match(zeilen[i], "%\\|\\|% *NA_real_",
                             info = paste(basename(f), i))
  }
  ## and the behaviour itself
  obj <- structure(list(), class = "keineKonvention")
  th <- suppressWarnings(as.numeric(ci(".countimp_theta")(obj, pflicht = FALSE) %||% NA_real_))
  testthat::expect_length(th, 1L)
  testthat::expect_true(is.na(th))
})

testthat::test_that("B43: glmmTMB::sigma is no longer called outside the layer", {
  pfad <- if (dir.exists("../../R")) "../../R" else "R"
  testthat::skip_if_not(dir.exists(pfad), "package sources not locatable")
  ## sigma() on a non-glmmTMB fit returns the residual SD, not a dispersion --
  ## the accessor refuses those classes, a direct call would not.
  for (f in list.files(pfad, pattern = "[.]R$", full.names = TRUE)) {
    zeilen <- sub("#.*$", "", readLines(f, warn = FALSE))
    testthat::expect_length(grep("glmmTMB::sigma", zeilen, fixed = TRUE), 0L)
  }
})

testthat::test_that("B43: .countimp_feld reports class, fields and package", {
  ## pscl is in Suggests, not in Imports: the one-level zero-inflation and
  ## hurdle methods need it, the rest of the package does not. Without this
  ## skip the test fails in an environment without Suggests -- and it fails on
  ## the missing dependency, not on the property under test.
  skip_if_not_installed("pscl")
  obj <- structure(list(a = 1), class = c("fremdklasse", "list"))
  testthat::expect_equal(ci(".countimp_feld")(obj, "a", "the a value"), 1)
  ## first present candidate wins -- how renames are absorbed
  testthat::expect_equal(ci(".countimp_feld")(obj, c("neu", "a"), "the a value"), 1)
  err <- tryCatch(ci(".countimp_feld")(obj, c("x", "y"), "the x value", pkg = "pscl"),
                  error = function(e) conditionMessage(e))
  testthat::expect_match(err, "fremdklasse")
  testthat::expect_match(err, "x, y")
  testthat::expect_match(err, "pscl")
  ## optional fields return NULL instead of stopping
  testthat::expect_null(ci(".countimp_feld")(obj, "x", "the x value", pflicht = FALSE))
})

testthat::test_that("B43: every Suggests package is guarded before use", {
  ## Packages in Suggests may be absent. Any :: call to them must be preceded
  ## by a guard, or sit inside a class-conditional branch that cannot be
  ## reached without the package.
  pfad <- if (dir.exists("../../R")) "../../R" else "R"
  testthat::skip_if_not(dir.exists(pfad), "package sources not locatable")
  desc <- if (file.exists("../../DESCRIPTION")) "../../DESCRIPTION" else "DESCRIPTION"
  testthat::skip_if_not(file.exists(desc), "DESCRIPTION not locatable")

  sug <- read.dcf(desc, fields = "Suggests")[1, 1]
  pkgs <- trimws(gsub("\\(.*?\\)", "", strsplit(sug, ",")[[1]]))
  pkgs <- setdiff(pkgs, c("testthat", ""))

  for (p in pkgs) {
    treffer <- character(0)
    for (f in list.files(pfad, pattern = "[.]R$", full.names = TRUE)) {
      zeilen <- sub("#.*$", "", readLines(f, warn = FALSE))
      idx <- grep(paste0(p, "::"), zeilen, fixed = TRUE)
      if (length(idx)) treffer <- c(treffer, paste0(basename(f), ":", idx))
    }
    if (!length(treffer)) next
    ## either a guard exists for this package, or all call sites are inside a
    ## branch conditional on an object class only that package produces
    alle <- unlist(lapply(list.files(pfad, pattern = "[.]R$", full.names = TRUE),
                          function(f) sub("#.*$", "", readLines(f, warn = FALSE))))
    hat_wache <- any(grepl(paste0("need_pkg\\(\"", p, "\"|have_pkg\\(\"", p,
                                  "\"|requireNamespace\\(\"", p, "\""), alle))
    ## A call site needs no guard when it is unreachable without the package:
    ## the only nlme call sits in an inherits(fit, "lme") branch, and objects
    ## of class lme exist only if nlme is installed (MASS::glmmPQL returns
    ## class glmmPQL/lme and imports nlme itself). Additionally wrapped in
    ## try(). Verified by inspection, not assumed.
    unerreichbar <- identical(p, "nlme") &&
      any(grepl("inherits\\(fit, *\"lme\"\\)", alle)) &&
      any(grepl("try\\(nlme::VarCorr", alle))
    testthat::expect_true(hat_wache || unerreichbar,
                          info = paste(p, "unguarded at", paste(treffer, collapse = ", ")))
  }
})

testthat::test_that("B61: no internal function is defined in two files", {
  ## Found the hard way: R/fitdiag.R defined .countimp_zi_formula, which already
  ## existed in R/impute1lzi.R with a different signature. Collation order
  ## decided which one survived, and the loser's callers failed with "argument
  ## \"type\" is missing" -- a message pointing at the wrong file entirely.
  ## R gives no warning for this, so it needs a test.
  pfad <- if (dir.exists("../../R")) "../../R" else "R"
  testthat::skip_if_not(dir.exists(pfad), "package sources not locatable")
  defs <- list()
  for (f in list.files(pfad, pattern = "[.]R$", full.names = TRUE)) {
    zeilen <- sub("#.*$", "", readLines(f, warn = FALSE))
    ## Top-level definitions only -- anchored at column 1 with no leading
    ## whitespace. An indented `ein_zug <- function()` is a local helper inside
    ## another function; two files may each have one without any collision,
    ## since neither is visible outside its parent. Only unindented definitions
    ## share the package namespace and can overwrite one another.
    hit <- regmatches(zeilen, regexpr("^[.]?[a-zA-Z][a-zA-Z0-9._]*\\s*<-\\s*function",
                                      zeilen))
    nam <- trimws(sub("\\s*<-\\s*function.*$", "", hit))
    for (n in nam) defs[[n]] <- c(defs[[n]], basename(f))
  }
  doppelt <- names(defs)[vapply(defs, function(x) length(unique(x)) > 1L, logical(1))]
  testthat::expect_length(doppelt, 0L)
  if (length(doppelt))
    testthat::fail(paste0("defined in several files: ",
                          paste(sprintf("%s (%s)", doppelt,
                                        vapply(defs[doppelt],
                                               function(x) paste(unique(x), collapse = ", "),
                                               character(1))),
                                collapse = "; ")))
})
