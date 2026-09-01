## Encapsulation layer for foreign packages -----------------------------------
##
## countimp reads fields out of objects returned by MASS, pscl, glmmTMB and
## stats -- fields that are part of no documented interface. A rename in any of
## those packages breaks countimp without countimp having changed. Two earlier
## findings were exactly this: B35 (pscl's theta slot names) and B42
## (summary.glm on the plain list glm.fit returns).
##
## This file collects the two kinds of foreign contact:
##
##   1. availability -- is the package installed at all (.countimp_need_pkg)
##   2. field access -- read one documented-by-us quantity out of a foreign
##      object, with an explicit error when the field is gone
##
## Nothing here fits a model. Fitting stays in the family cores; they call
## these accessors instead of reaching into the objects themselves.


## Which optional package does what, for error messages.
.countimp_pkg_zweck <- c(
  mice = "the non-count imputation methods (pmm, logreg, polyreg, ...)",
  pscl = "the single-level zero-inflated and hurdle methods (zip, zinb, hp, hnb)",
  nlme = "reading variance components from lme fits",
  nnet = "the multinomial method for unordered factors")


## Require an optional package, or stop with an actionable message.
##
## The message names the package, what it is needed for, how to install it, and
## -- where one exists -- the alternative that does not need it. Called before
## the first foreign call, not after a failure: a missing package must not
## surface as "failed to fit the imputation model" (B43).
.countimp_need_pkg <- function(pkg, what = NULL, ersatz = NULL) {
  if (requireNamespace(pkg, quietly = TRUE)) return(invisible(TRUE))
  if (is.null(what)) {
    what <- if (pkg %in% names(.countimp_pkg_zweck)) .countimp_pkg_zweck[[pkg]]
            else "this imputation method"
  }
  msg <- paste0("countimp: package '", pkg, "' is required for ", what,
                ", but it is not installed.\n",
                "  Install it with install.packages(\"", pkg, "\")")
  if (!is.null(ersatz)) msg <- paste0(msg, ", or use ", ersatz)
  stop(msg, ".", call. = FALSE)
}


## Is an optional package usable? For code that degrades instead of stopping.
.countimp_have_pkg <- function(pkg) {
  isTRUE(requireNamespace(pkg, quietly = TRUE))
}


## Read a field out of a foreign fit object.
##
## `feld` may name several candidates; the first one present wins, which is how
## renames across package versions are absorbed. When none is present the error
## names the object's class, the field(s) looked for, and the package that owns
## the object -- enough for a bug report without a debugger session.
.countimp_feld <- function(obj, feld, was, pkg = NULL, pflicht = TRUE) {
  for (f in feld) {
    if (!is.null(obj[[f]])) return(obj[[f]])
  }
  if (!pflicht) return(NULL)
  stop("countimp: cannot read ", was, " from an object of class ",
       paste(class(obj), collapse = "/"), ".\n",
       "  Looked for: ", paste(feld, collapse = ", "), ".\n",
       if (!is.null(pkg))
         paste0("  This object comes from package '", pkg,
                "' (version ", .countimp_pkg_version(pkg), "); a rename there ",
                "is the likely cause.\n"),
       "  Please report this, naming the package version.", call. = FALSE)
}


## Package version as a string, or "unknown".
.countimp_pkg_version <- function(pkg) {
  v <- try(as.character(utils::packageVersion(pkg)), silent = TRUE)
  if (inherits(v, "try-error")) "unknown" else v
}


## The dispersion parameter of a negative binomial fit.
##
## MASS::glm.nb calls it theta and reports it on the fit and in the summary;
## pscl::zeroinfl and pscl::hurdle carry it as a coefficient on the log scale
## (B35). glmmTMB reports it through sigma(). One accessor, three conventions.
## `pflicht = FALSE` returns NULL instead of stopping, for callers that carry
## their own message -- .countimp_zi_theta() does, and keeps it.
.countimp_theta <- function(fit, summ = NULL, pflicht = TRUE) {
  if (inherits(fit, "negbin"))
    return(.countimp_feld(if (is.null(summ)) fit else summ, "theta",
                          "the negative binomial dispersion theta", "MASS",
                          pflicht = pflicht))
  if (inherits(fit, c("zeroinfl", "hurdle"))) {
    ## pscl 1.5.9 carries theta on the fit itself and its standard error as
    ## SE.logtheta -- verified by inspection, NOT as a log-scale coefficient:
    ## coef() returns only count_* and zero_* entries. So the only fallback
    ## worth having is the log-scale slot pair, not a coefficient search.
    th <- fit[["theta"]]
    if (!is.null(th)) return(th)
    lg <- fit[["logtheta"]]
    if (!is.null(lg)) return(exp(unname(lg)))
    if (!pflicht) return(NULL)
    stop("countimp: cannot read theta from a pscl ",
         paste(class(fit), collapse = "/"), " fit (pscl ",
         .countimp_pkg_version("pscl"), ").\n",
         "  Looked for: $theta, $logtheta.\n",
         "  A rename in pscl is the likely cause; please report this, ",
         "naming the pscl version.", call. = FALSE)
  }
  ## sigma() only for classes where it IS the count dispersion. For an lm it
  ## returns the residual standard deviation, which is not a dispersion at all
  ## -- returning it would make the diagnostics flag a healthy fit.
  if (inherits(fit, "glmmTMB")) {
    s <- try(stats::sigma(fit), silent = TRUE)
    if (!inherits(s, "try-error") && length(s) == 1L && is.finite(s))
      return(as.numeric(s))
  }
  if (!pflicht) return(NULL)
  stop("countimp: cannot read a dispersion parameter from an object of class ",
       paste(class(fit), collapse = "/"), ".\n",
       "  Known conventions: MASS $theta, pscl $theta, glmmTMB sigma().",
       call. = FALSE)
}
