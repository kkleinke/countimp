## B38 -- countimp must not resolve its own internals through a hard-coded
## asNamespace("countimp").
##
## .countimp_find_method() used to do exactly that. Under the development load
## path (source() of R/*.R into the global environment) it therefore returned
## functions from the INSTALLED build: edited sources were never what ran, and
## two copies of .countimp_state existed side by side, so a state flag set by
## one was invisible to the other. The bug hid real failures in this suite --
## after the fix, three previously "passing" tests turned out to be asserting
## nothing.

testthat::test_that("B38: .countimp_home() finds the copy that is running", {
  home <- ci(".countimp_home")()
  testthat::expect_true(is.environment(home))
  ## it must be the environment countimp's own state lives in
  testthat::expect_true(exists(".countimp_state", envir = home, inherits = FALSE))

  ## Asserted as SAME STATE, not as identical environments. Under R CMD check
  ## testthat runs the tests in test_env("countimp"), an environment that copies
  ## the namespace's bindings and reports environmentName() == "countimp" itself.
  ## .countimp_home() then legitimately hands back that copy rather than the
  ## namespace, and expect_identical(home, ci_home()) failed with the useless
  ## message that two <env:namespace:countimp> differ -- while the property the
  ## fix is about held throughout. The bindings are copied, the objects they
  ## point to are not: one .countimp_state serves both, which is exactly what
  ## B38 exists to guarantee.
  st_home <- get(".countimp_state", envir = home,      inherits = FALSE)
  st_ci   <- get(".countimp_state", envir = ci_home(), inherits = FALSE)
  testthat::expect_identical(st_home, st_ci)
  ## and a write through one hull is visible through the other
  if (is.environment(st_home)) {
    st_home$.b38_probe <- 4711L
    on.exit(suppressWarnings(rm(".b38_probe", envir = st_home)), add = TRUE)
    testthat::expect_identical(st_ci$.b38_probe, 4711L)
  }
})

testthat::test_that("B38: method lookup and state agree on one copy", {
  ## The decisive property: the function the resolver hands out must write to
  ## the same state environment the caller can read. Under the old code these
  ## were different environments whenever the sources were sourced.
  f <- ci(".countimp_find_method")("nb")
  testthat::expect_true(is.function(f))
  testthat::expect_identical(environment(f), ci_home())

  st <- get(".countimp_state", envir = ci_home(), inherits = FALSE)
  st2 <- get(".countimp_state", envir = environment(f), inherits = FALSE)
  testthat::expect_identical(st, st2)
})

testthat::test_that("B38: a state flag set by the package is visible to a test", {
  ## This is what broke B04: the warning fired, the flag went up -- in the other
  ## copy. Set it here, have package code read it, and read it back.
  st <- get(".countimp_state", envir = ci_home(), inherits = FALSE)
  old <- st$underdisp_warned
  on.exit(st$underdisp_warned <- old, add = TRUE)

  st$underdisp_warned <- FALSE
  first <- ci(".countimp_warn_underdisp")
  testthat::expect_warning(first(0.5), "underdispers")
  testthat::expect_true(st$underdisp_warned)
  ## and the throttle holds: no second warning
  testthat::expect_silent(first(0.5))
})

testthat::test_that("B38: no source file hard-codes asNamespace(\"countimp\")", {
  cand <- c("../../R", "../../../R", "R", "paket/R")
  hit <- cand[vapply(cand, function(p)
    dir.exists(p) && length(list.files(p, pattern = "[.]R$")) > 0L, logical(1))]
  testthat::skip_if(length(hit) == 0L, "package sources not locatable")
  files <- list.files(hit[1], pattern = "[.]R$", full.names = TRUE)
  sites <- character(0)
  for (f in files) {
    src <- sub("#.*$", "", readLines(f, warn = FALSE))
    for (i in grep('asNamespace\\("countimp"\\)', src))
      sites <- c(sites, sprintf("%s:%d", basename(f), i))
  }
  ## the single legitimate use is the fallback inside .countimp_home()
  testthat::expect_true(all(grepl("^families[.]R:", sites)),
                        info = paste(sites, collapse = ", "))
  testthat::expect_true(length(sites) <= 1L,
                        info = paste(sites, collapse = ", "))
})

testthat::test_that("B38: no test file hard-codes countimp::: or the env name", {
  ## countimp:::.foo() in a test has the same defect as the package bug: it
  ## silently tests the installed build. ci() replaces it.
  ## This file necessarily names both patterns (in its title and in the search
  ## itself), so it excludes itself from the scan.
  fs <- setdiff(list.files(".", pattern = "^test-.*[.]R$"),
                "test-B38-load-path.R")
  testthat::skip_if(length(fs) == 0L, "tests not locatable")
  bad_ns <- character(0)
  bad_nm <- character(0)
  bad_bare <- character(0)
  for (f in fs) {
    src <- sub("#.*$", "", readLines(f, warn = FALSE))
    if (any(grepl("countimp:::", src, fixed = TRUE))) bad_ns <- c(bad_ns, f)
    ## environmentName(...) compared against the literal "countimp"
    if (any(grepl('environmentName\\([^)]*\\)[^"]*"countimp"', src)))
      bad_nm <- c(bad_nm, f)
    ## A BARE call to a PACKAGE INTERNAL is the same defect and fails worse:
    ## under R CMD check nothing of the package is in scope, so the test ERRORS
    ## with "could not find function" instead of failing on its assertion, and
    ## the message points at the wrong thing. Found in test-B69 line 109 by a
    ## full-suite run, where the countimp::: scan had nothing to grep for.
    ##
    ## The decision is by NAMESPACE MEMBERSHIP, not by the .countimp_ prefix:
    ## a prefix rule produced two false positives -- test-B34 carries the
    ## pattern inside a search string, and .countimp_quelle() is a helper
    ## defined in helper-internals.R that borrows the prefix. Neither is a
    ## package internal, so neither is the bug this guards against.
    ## String literals are blanked first so a pattern inside quotes is not a
    ## call, and ci(".countimp_x") is untouched for the same reason.
    ohne_str <- gsub('"[^"]*"', '""', src)
    kand <- unlist(regmatches(ohne_str, gregexpr(
      "[.]countimp_[[:alnum:]._]*[[:space:]]*\\(", ohne_str)))
    kand <- unique(sub("[[:space:]]*\\($", "", kand))
    ## No file-wide exemption for files that define functions: it would excuse
    ## a real bare call in any file that happens to define a helper. A local
    ## definition that SHADOWS a package internal would be flagged here, and
    ## should be -- shadowing an internal in a test hides which one ran.
    if (length(kand) &&
        any(vapply(kand, exists, logical(1),
                   where = asNamespace("countimp"), inherits = FALSE)))
      bad_bare <- c(bad_bare, f)
  }
  testthat::expect_equal(bad_ns, character(0))
  testthat::expect_equal(bad_nm, character(0))
  testthat::expect_equal(bad_bare, character(0))
})
