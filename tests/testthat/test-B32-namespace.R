## B32 -- every call into a foreign package must name its namespace.
##
## Up to 3.0.0 the eight two-level NB methods called rnegbin() unqualified. In
## the installed package import(MASS) hides that: the function sits in the
## package's imports frame. It surfaces the moment the code is not run as an
## installed package, or the moment the import form changes -- and moving from
## import(MASS) to importFrom() is exactly what a lean dependency surface wants.
##
## This test is a standing guard, not a one-off check: it scans every file in
## R/ and fails as soon as a new unqualified call appears. That is the point --
## the user's stated priority is that an update to a foreign package must not
## be able to disturb countimp.

testthat::test_that("no foreign function is called without its namespace", {
  ## pscl is in Suggests, not in Imports: the one-level zero-inflation and
  ## hurdle methods need it, the rest of the package does not. Without this
  ## skip the test fails in an environment without Suggests -- and it fails on
  ## the missing dependency, not on the property under test.
  skip_if_not_installed("pscl")
  ## An installed package has no R/*.R files (they are byte-compiled into
  ## R/countimp.rdb), so the sources must be found in the tree. Candidates
  ## cover: testthat's own wd (tests/testthat), R CMD check, and the repo root.
  cand <- c("../../R", "../../../R", "R", "paket/R",
            "../../paket/R", "../../../paket/R")
  hit <- cand[vapply(cand, function(p)
    dir.exists(p) && length(list.files(p, pattern = "[.]R$")) > 0L, logical(1))]
  testthat::skip_if(length(hit) == 0L, "package sources not locatable")
  files <- list.files(hit[1], pattern = "[.]R$", full.names = TRUE)
  ## guard against the guard: the scan is worthless if it sees no engines
  testthat::expect_true(any(grepl("impute2l", basename(files))))

  ## Functions countimp uses from foreign packages. stats/utils/methods are
  ## base-R and always attached, so they are not at risk here.
  foreign <- c("rnegbin", "glm.nb", "ginv", "mvrnorm", "polr",       # MASS
               "zeroinfl", "hurdle",                                  # pscl
               "glmmTMB", "fixef", "ranef", "VarCorr",                # glmmTMB
               "sdreport",                                            # TMB
               "multinom",                                            # nnet
               "lme", "lmeControl")                                    # nlme

  offenders <- character(0)
  for (f in files) {
    src <- readLines(f, warn = FALSE)
    ## drop comments (roxygen and plain) -- prose may name a function freely
    src <- sub("#.*$", "", src)
    ## drop string literals: a function name inside an error message or a
    ## fallback reason is text, not a call. Both quote styles, non-greedy.
    src <- gsub('"[^"]*"', '""', src)
    src <- gsub("'[^']*'", "''", src)
    for (fn in foreign) {
      pat <- paste0("(^|[^[:alnum:]._:])", gsub("[.]", "[.]", fn), "[[:space:]]*\\(")
      hit <- grep(pat, src)
      for (i in hit) {
        ## qualified calls carry "pkg::" immediately before the name
        if (!grepl(paste0("::[[:space:]]*", gsub("[.]", "[.]", fn)), src[i])) {
          offenders <- c(offenders,
                         sprintf("%s:%d  %s", basename(f), i, trimws(src[i])))
        }
      }
    }
  }
  testthat::expect_equal(offenders, character(0),
                         info = paste(offenders, collapse = "\n"))
})

testthat::test_that("stats::sigma is used, not glmmTMB::sigma", {
  ## sigma() has been a generic in stats since R 3.3. Calling the stats generic
  ## keeps working through S3 dispatch even if glmmTMB stops re-exporting its
  ## own sigma(); both return the same value for an nbinom2 fit.
  src <- paste(deparse(ci(".countimp_stabilize")), collapse = "\n")
  testthat::expect_match(src, "stats::sigma", fixed = TRUE)
})

testthat::test_that("the declared dependencies are the ones actually used", {
  ## A dependency that is declared but unused is dead weight; one that is used
  ## but undeclared is a latent failure. Both are drift risks.
  d <- system.file("DESCRIPTION", package = "countimp")
  testthat::skip_if(!nzchar(d) || !file.exists(d), "DESCRIPTION not available")
  imp <- read.dcf(d, fields = "Imports")[1, 1]
  testthat::skip_if(is.na(imp), "no Imports field")
  imp <- trimws(unlist(strsplit(imp, ",")))
  imp <- sub("[[:space:]]*\\(.*$", "", imp)
  ## mice must NOT be an import: independence from it is a project goal
  testthat::expect_false("mice" %in% imp)
  ## the packages the engines need must be declared
  for (p in c("glmmTMB", "MASS", "stats")) testthat::expect_true(p %in% imp, info = p)
})
