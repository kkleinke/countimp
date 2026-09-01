## zzz.R -- load-time checks on the dependency stack
##
## Motivation. During the 2.1.0 work half of the multilevel methods crashed
## the R session with a segfault, not an error -- no message, no traceback,
## no way to catch it with try(). The cause was outside countimp: the
## installed glmmTMB binary had been compiled against TMB 1.9.19 while
## TMB 1.9.23 was on the library path. glmmTMB detects this and warns, but
## the warning is easy to miss, and the consequence is a hard crash inside
## compiled code the moment a model with random effects is fitted.
##
## A user hitting this sees countimp die mid-imputation and reasonably
## blames countimp. The check below moves the diagnosis to load time and
## says what to do about it.

## Implementation note. Up to 2.6.0 this called glmmTMB:::check_dep_version(),
## which is exactly the kind of dependence on a foreign internal that this
## package is trying to shed: a rename in glmmTMB turns the safety check into
## an error, and R CMD check flags the ':::' as a NOTE. There is no exported
## equivalent -- neither glmmTMB nor TMB exports the version the binary was
## built against -- so the comparison is done here from public data where
## possible, and the one value that is only available internally is read
## defensively, never with ':::'.
##
## Returns NA when the build version cannot be determined. That is the honest
## answer: absence of evidence, not evidence of a match. Callers stay silent
## in that case rather than claim the backend is sound.

.countimp_tmb_build_version <- function() {
  ns <- tryCatch(asNamespace("glmmTMB"), error = function(e) NULL)
  if (is.null(ns)) return(NULL)
  ## mget() with ifnotfound cannot fail if glmmTMB drops or renames the object
  v <- mget(".TMB.build.version", envir = ns,
            ifnotfound = list(NULL), inherits = FALSE)[[1]]
  if (is.null(v)) return(NULL)
  tryCatch(as.character(v), error = function(e) NULL)
}

.countimp_check_tmb <- function() {
  if (!requireNamespace("glmmTMB", quietly = TRUE)) return(invisible(NA))
  built <- .countimp_tmb_build_version()
  if (is.null(built)) return(invisible(NA))     # cannot tell; say nothing
  mismatch <- tryCatch(
    utils::packageVersion("TMB") != built,
    error = function(e) NA)
  if (is.na(mismatch)) return(invisible(NA))
  if (isTRUE(mismatch)) {
    packageStartupMessage(sprintf(paste0(
      "countimp: glmmTMB %s was compiled against a different version of TMB\n",
      "  than the one installed (TMB %s). Multilevel imputation can CRASH the\n",
      "  R session (segfault) rather than raise a catchable error.\n",
      "  Fix: install.packages('glmmTMB', type = 'source')\n",
      "  See ?glmmTMB::reinstalling"),
      as.character(utils::packageVersion("glmmTMB")),
      as.character(utils::packageVersion("TMB"))))
    return(invisible(FALSE))
  }
  invisible(TRUE)
}

.onAttach <- function(libname, pkgname) {
  .countimp_check_tmb()        # ABI consistency of the compiled backend
  .countimp_check_backend()    # family constructors still present (families.R)
}
