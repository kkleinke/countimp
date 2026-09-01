## B33 -- the manual page must name the distribution the code actually fits.
##
## Three of the four two-level Poisson methods were documented as "NB
## regression variant" while fitting family = "poisson". A user reading the
## help page would pick .poisson.boot believing they get an NB model -- that
## is, believing they model the overdispersion that motivated the choice.
##
## The check runs against the function bodies, so it holds for the consolidated
## wrappers and would also have caught the defect in the old copies.

testthat::test_that("B33: family in the body matches the family in the name", {
  poisson.methods <- c("mice.impute.2l.poisson", "mice.impute.2l.poisson.noint",
                       "mice.impute.2l.poisson.boot",
                       "mice.impute.2l.poisson.noint.boot")
  nb.methods <- c("mice.impute.2l.nb", "mice.impute.2l.nb.noint",
                  "mice.impute.2l.nb.boot", "mice.impute.2l.nb.noint.boot",
                  "mice.impute.2l.nb2", "mice.impute.2l.nb2.noint",
                  "mice.impute.2l.nb2.boot", "mice.impute.2l.nb2.noint.boot")
  for (m in poisson.methods) {
    bod <- paste(deparse(body(get(m))), collapse = " ")
    testthat::expect_match(bod, '"poisson"', fixed = TRUE, info = m)
    testthat::expect_false(grepl('"nbinom2"', bod, fixed = TRUE), info = m)
  }
  for (m in nb.methods) {
    bod <- paste(deparse(body(get(m))), collapse = " ")
    testthat::expect_match(bod, '"nbinom2"', fixed = TRUE, info = m)
    testthat::expect_false(grepl('"poisson"', bod, fixed = TRUE), info = m)
  }
})

testthat::test_that("B33: the nb2 aliases ARE their originals", {
  ## Not "behave the same" -- literally the same function object. Up to 3.0.0
  ## they were full copies of the bodies, which is how copies drift apart.
  testthat::expect_identical(mice.impute.2l.nb2,             mice.impute.2l.nb)
  testthat::expect_identical(mice.impute.2l.nb2.noint,       mice.impute.2l.nb.noint)
  testthat::expect_identical(mice.impute.2l.nb2.boot,        mice.impute.2l.nb.boot)
  testthat::expect_identical(mice.impute.2l.nb2.noint.boot,  mice.impute.2l.nb.noint.boot)
})

testthat::test_that("B33: the manual page text agrees with the family", {
  ## Read the generated Rd rather than the sources: that is what the user sees.
  rd <- system.file("man", "mice.impute.2l.poisson.Rd", package = "countimp")
  if (!nzchar(rd) || !file.exists(rd)) {
    cand <- c("man/mice.impute.2l.poisson.Rd", "../../man/mice.impute.2l.poisson.Rd",
              "paket/man/mice.impute.2l.poisson.Rd")
    rd <- cand[file.exists(cand)][1]
  }
  testthat::skip_if(is.na(rd) || !file.exists(rd), "Rd file not available")
  txt <- paste(readLines(rd, warn = FALSE), collapse = "\n")

  ## In the \section{Methods} / \describe block every \item names a method and
  ## its description. A Poisson method must not be called an NB variant.
  items <- regmatches(txt, gregexpr("\\\\item\\{[^}]*\\}\\{[^}]*\\}", txt))[[1]]
  testthat::skip_if(length(items) == 0L, "no \\item entries in the Rd")
  bad <- character(0)
  for (it in items) {
    nm <- sub("^\\\\item\\{([^}]*)\\}.*$", "\\1", it)
    ds <- sub("^\\\\item\\{[^}]*\\}\\{([^}]*)\\}$", "\\1", it)
    if (!grepl("poisson", nm, fixed = TRUE)) next
    if (grepl("nb2", nm, fixed = TRUE)) next
    ## a Poisson entry may not be described as NB
    if (grepl("\\bNB\\b|negative binomial", ds) && !grepl("Poisson", ds))
      bad <- c(bad, paste(nm, "->", ds))
  }
  testthat::expect_equal(bad, character(0), info = paste(bad, collapse = "\n"))
})

testthat::test_that("B33: the header prose is not cut off mid-sentence", {
  ## Two sentences in the old header ended abruptly ("... either using a
  ## Bayesian" then a bare "#'", and "... rather than" then a new \enumerate).
  ## A line ENDING in a conjunction is normal in wrapped prose -- the defect is
  ## a dangling line whose sentence never continues, i.e. the next roxygen line
  ## is empty or starts a new block. That is what this checks.
  cand <- c("R/mice.impute.2l.poisson.R", "../../R/mice.impute.2l.poisson.R",
            "paket/R/mice.impute.2l.poisson.R")
  f <- cand[file.exists(cand)][1]
  testthat::skip_if(is.na(f), "source file not locatable")
  ln <- grep("^#'", readLines(f, warn = FALSE), value = TRUE)
  txt <- sub("^#'[[:space:]]?", "", ln)

  dangling <- character(0)
  ## line 1 is the Rd title and line 2 the blank separator: a title carries no
  ## terminal punctuation by convention, so the scan starts at line 3.
  for (i in seq.int(3L, length(txt) - 1L)) {
    ## does this line end mid-sentence? (no terminal punctuation, no brace)
    if (!grepl("[[:alnum:],]$", txt[i])) next
    if (grepl("^[[:space:]]*$", txt[i])) next
    ## roxygen TAG lines are not prose: @describeIn/@param/@aliases text
    ## legitimately ends on a word and is followed by the next tag.
    if (grepl("^[[:space:]]*@", txt[i])) next
    nxt <- txt[i + 1L]
    ## a continuation must be prose; an empty line or a new roxygen tag or a
    ## new \enumerate/\itemize block means the sentence was abandoned
    broken <- grepl("^[[:space:]]*$", nxt) ||
              grepl("^[[:space:]]*@", nxt) ||
              grepl("^[[:space:]]*\\\\(enumerate|itemize|describe)", nxt)
    if (broken) dangling <- c(dangling, sprintf("line %d: %s", i, txt[i]))
  }
  testthat::expect_equal(dangling, character(0),
                         info = paste(dangling, collapse = "\n"))
})


## The DESCRIPTION carries a Date -------------------------------------------
##
## Reported from the simulation side on 1 September 2026: citation() rendered
## the package as "Kleinke K (????)" and warned "could not determine year for
## 'countimp' from package DESCRIPTION file". Without a Date field, citation()
## has no year to use -- so anyone citing the package in a paper gets a
## question mark where the year belongs.
##
## Kept as a test because the field is easy to lose in a version bump, and the
## symptom only shows when someone actually calls citation().

test_that("B33: citation() can determine the year", {
  d <- utils::packageDescription("countimp")
  expect_false(is.null(d$Date))
  expect_match(d$Date, "^[0-9]{4}-[0-9]{2}-[0-9]{2}$")
  ## and it must actually reach citation(): no warning, and a year in the text
  zit <- testthat::expect_no_warning(utils::citation("countimp"))
  expect_match(format(zit, style = "text"), "\\(20[0-9]{2}\\)")
})
