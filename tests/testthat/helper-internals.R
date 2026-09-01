## Resolve package internals against the copy of countimp actually under test.
##
## The suite runs under two load paths:
##
##   1. R CMD check / test_check("countimp") -- the package is installed and
##      loaded, and the internals live in its namespace.
##   2. A development runner that source()s R/*.R into the global environment,
##      which is how a fix is checked before a rebuild.
##
## Writing countimp:::.foo() hard-codes path 1. Under path 2 it silently
## resolves against the INSTALLED build, so an edited source file is not what
## gets tested -- the tests pass while testing yesterday's code, and a new
## internal that does not exist in the installed namespace makes them error out
## for the wrong reason. Worse, a test that mixes countimp:::.foo() with an
## unqualified bar() ends up with two copies of the package state.
##
## ci("name") prefers the sourced definition when there is one and falls back to
## the namespace otherwise, so a test says what it means under both paths.
## Functions carry their own enclosing environment, so state objects such as
## .countimp_diag stay consistent with the functions that write to them.
ci <- function(name) {
  if (exists(name, envir = globalenv(), inherits = FALSE))
    return(get(name, envir = globalenv(), inherits = FALSE))
  ns <- tryCatch(asNamespace("countimp"), error = function(e) NULL)
  if (is.null(ns))
    stop("countimp is neither sourced nor installed, cannot resolve '", name,
         "'", call. = FALSE)
  if (!exists(name, envir = ns, inherits = FALSE))
    stop("'", name, "' exists in neither the sourced sources nor the installed ",
         "namespace. If it is new, rebuild the package or source R/*.R.",
         call. = FALSE)
  get(name, envir = ns, inherits = FALSE)
}

## Which path are we on? Useful in skip conditions and failure messages.
##
## Decided on a FUNCTION, not on the state object .countimp_diag: under
## R CMD check the tests run with the package attached but nothing of it in the
## global environment, and .countimp_diag is an environment that only lands
## there when R/*.R is sourced. Keying on it was right for path 2 but made
## path 1 depend on load order -- if any earlier test had put a name of that
## spelling into globalenv(), ci_home() switched to the wrong environment and
## B38 failed with the confusing message that two identical-looking namespaces
## differ. A function that every load path defines is the stable marker.
ci_path <- function() {
  if (exists(".countimp_home", envir = globalenv(), inherits = FALSE) &&
      is.function(get(".countimp_home", envir = globalenv(), inherits = FALSE)))
    "sourced" else "installed"
}

## The environment countimp's own functions live in on the current load path:
## the namespace when installed, the global environment when sourced. Use this
## instead of the literal "countimp" when a test asserts that a function came
## from countimp rather than from mice -- the literal only holds on one of the
## two paths (B38).
ci_home <- function() {
  if (identical(ci_path(), "sourced")) globalenv() else asNamespace("countimp")
}

## Replace an internal with a spy for the duration of one test, on both load
## paths. assign(..., envir = ci_home()) alone is not enough: under R CMD check
## the namespace is SEALED and every binding locked, so the assignment fails
## with "cannot change value of locked binding" -- B40 errored out for exactly
## that reason while passing in the sourced runner. unlockBinding() opens it,
## and the original is put back on exit either way.
##
## Returns invisibly; call it inside test_that() and the restore happens when
## the block leaves. `env` defaults to the caller's frame so that on.exit is
## registered in the test block, not in this helper.
ci_spy <- function(name, spy, env = parent.frame()) {
  home <- ci_home()
  orig <- get(name, envir = home, inherits = FALSE)
  gesperrt <- bindingIsLocked(name, home)
  if (gesperrt) unlockBinding(name, home)
  assign(name, spy, envir = home)
  ## The restore expression must carry its VALUES, not the names `name`/`home`/
  ## `orig` -- those live in this helper's frame and are gone by the time
  ## on.exit fires in the caller. bquote(.(x)) splices them in.
  restore <- bquote({
    if (bindingIsLocked(.(name), .(home))) unlockBinding(.(name), .(home))
    assign(.(name), .(orig), envir = .(home))
    if (.(gesperrt)) lockBinding(.(name), .(home))
  })
  do.call(base::on.exit, list(restore, add = TRUE), envir = env)
  invisible(orig)
}

## Count data with a missing outcome, for the B55 dispatch tests. Kept here so
## the four test blocks measure the same data and no block silently changes it.
.b55_daten <- function(seed, n = 150L) {
  set.seed(seed)
  x <- stats::rnorm(n)
  d <- data.frame(x = x, y = stats::rpois(n, exp(0.6 + 0.3 * x)))
  d$y[sample.int(n, size = round(0.2 * n))] <- NA
  d
}

## Read one file from the package tree, or return NULL if the tree is not
## reachable. An installed package has no R/*.R and no man/*.Rd next to the
## tests, so every source-inspecting check must be able to skip. Candidates
## cover testthat's own wd, R CMD check's test dir, and the repo root -- the
## same list the older tests spell out inline (test-B32).
.countimp_quelle <- function(relpfad) {
  wurzeln <- c("../..", "../../..", ".", "paket", "../../paket", "../../../paket")
  for (w in wurzeln) {
    p <- file.path(w, relpfad)
    if (file.exists(p)) return(readLines(p, warn = FALSE))
  }
  NULL
}

## Count data with a clean Poisson signal and gaps in y, for B58.
.b58_daten <- function(n = 80L, seed = 4L) {
  set.seed(seed)
  d <- data.frame(x = stats::rnorm(n), y = stats::rpois(n, 3))
  d$y[sample.int(n, 20L)] <- NA
  d
}

## Overdispersed count data with gaps in y. Shared by B40 (engine switches)
## and B56 (retry loop): a helper defined inside a test file is invisible to
## every other test file, so it belongs here.
zaehl_daten <- function(n = 200, theta = 1.5, seed = 5) {
  set.seed(seed)
  x1 <- stats::rnorm(n); x2 <- stats::rnorm(n)
  y <- MASS::rnegbin(n, mu = exp(1 + 0.6 * x1 - 0.4 * x2), theta = theta)
  d <- data.frame(y = y, x1 = x1, x2 = x2)
  d$y[sample.int(n, round(0.25 * n))] <- NA
  d
}
