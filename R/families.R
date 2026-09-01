## families.R -- single point of contact for distribution families
##
## WHY THIS FILE EXISTS
##
## Up to version 2.0.8 countimp specified families inline in 76 places, in
## three different notations:
##
##   family = list(family = "nbinom2", link = "log")   (44x, deprecated)
##   family = "binomial"                               (32x)
##   family = "gaussian"                               ( 1x)
##
## The list notation has been deprecated in glmmTMB since 1.1.5. Every fit
## emitted two warnings:
##
##   "specifying 'family' as a plain list is deprecated"
##   "some components missing from 'family': downstream methods may fail"
##
## Two consequences. First, the warnings drowned the diagnostic signal: in a
## multilevel validation run 57 of 57 draws carried warnings, all of them
## self-inflicted. Second, and more serious, the second message is a promise
## of future breakage. glmmTMB 1.1.14 still repairs the incomplete family
## object internally -- verified: family()$linkinv, $variance, predict(),
## residuals(), simulate(), sigma() and confint() all return identical values
## for both notations. A later release need not keep doing so.
##
## Routing every family through one function means the next glmmTMB API change
## is a single-line fix here rather than a search across 76 call sites. That
## is the point: countimp should not have to be edited in dozens of places
## because a dependency changed its mind.
##
## Note that "poisson", "binomial" and "gaussian" come from stats, while the
## count-specific families come from glmmTMB. Qualifying each explicitly means
## the package no longer relies on import(glmmTMB) in NAMESPACE putting the
## right names on the search path -- and it makes visible which dependency
## each family actually belongs to.

.countimp_family <- function(name, link = "log") {
  switch(name,
    ## --- stats ---------------------------------------------------------
    poisson  = stats::poisson(link = link),
    ## the hurdle/zero part is a logit model; accept the default link call
    ## with link = "log" so callers need not special-case it
    binomial = stats::binomial(link = if (link == "log") "logit" else link),
    gaussian = stats::gaussian(),
    ## --- glmmTMB -------------------------------------------------------
    nbinom1           = glmmTMB::nbinom1(link = link),
    nbinom2           = glmmTMB::nbinom2(link = link),
    truncated_poisson = glmmTMB::truncated_poisson(link = link),
    truncated_nbinom1 = glmmTMB::truncated_nbinom1(link = link),
    truncated_nbinom2 = glmmTMB::truncated_nbinom2(link = link),
    ## compois is the two-level route only. Single-level COM-Poisson is fitted
    ## by countimp's own core (R/compois.R) so that it survives on machines
    ## without a C++ toolchain; random effects genuinely need glmmTMB.
    compois           = glmmTMB::compois(link = link),
    stop(sprintf(
      "countimp: unknown distribution family '%s'. Available: %s", name,
      paste(c("poisson", "binomial", "gaussian", "nbinom1", "nbinom2",
              "truncated_poisson", "truncated_nbinom1", "truncated_nbinom2",
              "compois"),
            collapse = ", ")),
      call. = FALSE))
}

## Startup check. glmmTMB renamed the family constructors once already
## (nbinom2 was previously reached through the list notation only). If a
## future version drops or renames one, fail loudly at load time with a
## message that names the problem, instead of failing deep inside an
## imputation run with a missing-object error.

.countimp_check_backend <- function() {
  needed <- c("nbinom1", "nbinom2", "truncated_poisson",
              "truncated_nbinom1", "truncated_nbinom2", "compois")
  missing <- needed[!vapply(needed, function(f)
    exists(f, envir = asNamespace("glmmTMB"), inherits = FALSE), logical(1))]
  if (length(missing)) {
    warning(sprintf(paste0(
      "countimp: glmmTMB %s does not provide the family constructor(s) %s.\n",
      "  Multilevel count imputation will fail. Please report this at\n",
      "  https://github.com/kkleinke/countimp/issues -- countimp needs to be\n",
      "  adapted to the new glmmTMB API."),
      as.character(utils::packageVersion("glmmTMB")),
      paste(sQuote(missing), collapse = ", ")), call. = FALSE)
    return(invisible(FALSE))
  }
  invisible(TRUE)
}

## ---------------------------------------------------------------------------
## Starting values (since 2.6.0)
##
## The FCS loop needs an initial fill-in before the first real draw. Up to
## 2.5.0 this was mice::mice.impute.sample(). That single call was the last
## HARD dependency on mice inside the imputation engine -- a package that is
## otherwise only needed for non-count methods and for pooling. Hot-deck
## sampling from the observed values is four lines, so the dependency is not
## worth carrying: mice moved to Suggests in 2.6.0.
##
## Semantics are those of mice.impute.sample(): draw with replacement from the
## observed part of y. Identical in distribution, not bit-identical (a
## different number of RNG draws is consumed).
## ---------------------------------------------------------------------------

.countimp_impute_sample <- function(y, ry, wy = NULL) {
  if (is.null(wy)) wy <- !ry
  obs <- y[ry]
  if (!length(obs))
    stop("Cannot draw starting values: variable has no observed cases.",
         call. = FALSE)
  if (length(obs) == 1L) rep(obs, sum(wy))
  else sample(obs, size = sum(wy), replace = TRUE)
}

## Guard for the optional mice dependency -------------------------------------
.countimp_need_mice <- function(what) {
  if (!requireNamespace("mice", quietly = TRUE))
    stop("Package 'mice' is required for ", what, ", but is not installed.\n",
         "  Install it with install.packages(\"mice\"), or use countimp's own\n",
         "  count-data methods, which no longer depend on mice.", call. = FALSE)
  invisible(TRUE)
}

## Resolving imputation methods to functions ----------------------------------
##
## Methods are named by a string ("nb", "hp", "pmm", ...) and turned into a
## function name ("mice.impute.nb"). Three places can supply that function:
##
##   1. countimp's own namespace   -- the count-data methods
##   2. the user's workspace       -- own methods, and methods from packages
##                                    the user has attached
##   3. mice's namespace           -- pmm, logreg, polyreg, ... for the
##                                    non-count variables in the data
##
## Earlier versions used exists(mode = "function") followed by do.call() on
## the *name*. Both resolve through the search path, so a method living in
## mice was found only when mice happened to be attached -- installed was not
## enough. Under R CMD check nothing is attached, which made the documented
## examples fail. Resolving the function object explicitly removes that
## dependence on the search path.

## Resolution order (since 3.0.0)
##
##   1. countimp's own count methods       mice.impute.<m>     in countimp
##   2. countimp's own base methods        countimp.impute.<m> in countimp
##   3. the user's own function            mice.impute.<m>     in globalenv
##   4. mice's method                      mice.impute.<m>     in mice
##
## Rationale for this order. (1) before (4) is what makes countimp's count
## models take precedence over same-named mice methods. (2) before (3) and (4)
## means a plain method = "pmm" works with mice absent -- which is the point of
## the base methods in basemethods.R. (3) before (4) lets a user override a
## method by defining mice.impute.<m> in their session, the convention mice
## itself established.
##
## The user can reverse (2) and (4) with options(countimp.methods = "mice"),
## for instance to reproduce results obtained with mice's implementations.
## Two of countimp's base methods differ from mice's on purpose -- polyreg and
## polr draw the coefficients from their posterior, mice's do not -- so this
## switch is a real reproducibility control, not cosmetic.
## Where does this copy of countimp live?
##
## .countimp_find_method() has to look up sibling functions by name, and the
## obvious asNamespace("countimp") is wrong in one important case: while the
## package is being developed the sources are source()d into the global
## environment, and asNamespace() then silently returns the INSTALLED build.
## Method lookup would resolve to the installed code while every other call in
## the session runs the edited sources -- two copies of the package, each with
## its own .countimp_state, so a fix appears to have no effect and a state flag
## set by one is invisible to the other.
##
## environment(sys.function()) answers the question correctly on both paths: it
## is the namespace when the package is installed and loaded, and the global
## environment when the sources were source()d. The asNamespace() fallback
## covers the case of being called from an environment whose parent chain has
## been detached from either.
.countimp_home <- function() {
  e <- environment(sys.function(-1L))
  if (is.environment(e) && exists(".countimp_state", envir = e, inherits = FALSE))
    return(e)
  e <- topenv(parent.frame())
  if (is.environment(e) && exists(".countimp_state", envir = e, inherits = FALSE))
    return(e)
  tryCatch(asNamespace("countimp"), error = function(err) globalenv())
}

.countimp_find_method <- function(method) {
  fname <- paste0("mice.impute.", method)
  cname <- paste0("countimp.impute.", method)
  ns    <- .countimp_home()

  ## (1) own count methods
  if (exists(fname, envir = ns, mode = "function", inherits = FALSE))
    return(get(fname, envir = ns, mode = "function"))

  prefer_mice <- identical(getOption("countimp.methods", "countimp"), "mice")

  own_base <- function() {
    if (exists(cname, envir = ns, mode = "function", inherits = FALSE))
      get(cname, envir = ns, mode = "function") else NULL
  }
  from_mice <- function() {
    if (requireNamespace("mice", quietly = TRUE)) {
      mns <- asNamespace("mice")
      if (exists(fname, envir = mns, mode = "function", inherits = FALSE))
        return(get(fname, envir = mns, mode = "function"))
    }
    NULL
  }

  if (prefer_mice) {
    f <- from_mice()
    if (!is.null(f)) return(f)
  } else {
    ## (2) own base methods
    f <- own_base()
    if (!is.null(f)) return(f)
  }

  ## (3) the user's own function
  if (exists(fname, envir = globalenv(), mode = "function", inherits = TRUE))
    return(get(fname, envir = globalenv(), mode = "function"))

  ## (4) whichever of the two was not tried above
  f <- if (prefer_mice) own_base() else from_mice()
  if (!is.null(f)) return(f)

  NULL
}

## Argument checking for countimp() -----------------------------------------
##
## The engine signature ends in `...`, so an unrecognised argument is absorbed
## without a word. That is tolerable for arguments passed on to a method, and
## dangerous for arguments the caller believes are steering the imputation:
## a typo (`maxti = 20`) or a mice-only argument (`ignore =`) would silently
## produce imputations under different settings than requested.
##
## The check therefore distinguishes three cases:
##   * known engine arguments        -> accepted
##   * known mice-only arguments     -> error naming the countimp equivalent
##   * anything else                 -> error listing the closest known name
##
## Method arguments are NOT enumerated by hand. The first version of this check
## carried a hand-written list of engine arguments only, and rejected `EV`,
## `bounds`, `donors` and `theta` -- documented arguments of countimp's own
## methods -- with "Unknown argument in countimp(): `EV` (did you mean `m`?)"
## (B76). A hand-written list of everything the package accepts is guaranteed to
## drift behind the package.
##
## The method arguments are therefore read from the formals of every
## mice.impute.* function in the namespace at check time, so a new method with a
## new argument is accepted the moment it exists. Arguments consumed inside the
## fitting helpers rather than declared by the method (they travel in `...`) are
## the one part that still needs naming, in .countimp_dots_args below.
##
## Genuinely unknown names remain an error with a means of escape:
## options(countimp.check.args = FALSE).
.countimp_engine_args <- c(
  "data", "m", "method", "predictorMatrix", "where", "visitSequence",
  "form", "post", "defaultMethod", "maxit", "diagnostics", "printFlag",
  "seed", "imputationMethod", "defaultImputationMethod", "data.init",
  ## formula interface
  "formulas", "family", "zero", "draw")


## Arguments consumed by the fitting helpers rather than declared in a method
## signature: they arrive in the method's `...` and are read further down. These
## cannot be discovered from formals(), so they are named here -- the only
## hand-maintained part of the allowlist.
.countimp_dots_args <- c("quiet", "theta", "theta.min", "se.max", "EV")


## Arguments accepted by countimp(): engine arguments, plus every argument
## declared by any mice.impute.* method in the namespace, plus the pass-through
## arguments above. Computed on each call rather than cached, because a user can
## define a method in their own workspace between two calls.
.countimp_known_args <- function() {
  ## .countimp_home(), never asNamespace("countimp"): under the development load
  ## path the latter reads the INSTALLED build, so a method added to the sources
  ## would not be seen and its arguments would be rejected (B38).
  meth <- tryCatch({
    ns <- .countimp_home()
    ms <- grep("^mice\\.impute\\.", ls(ns), value = TRUE)
    unlist(lapply(ms, function(m) names(formals(get(m, envir = ns)))))
  }, error = function(e) character(0))
  meth <- setdiff(unique(meth), c("y", "ry", "x", "wy", "...", "type"))
  unique(c(.countimp_engine_args, meth, .countimp_dots_args))
}

## mice-only arguments, with the countimp route to the same goal.
.countimp_mice_only <- c(
  ignore   = paste("use `where` to control which cells are imputed;",
                   "`ignore` selects rows used for fitting and has no",
                   "countimp equivalent yet"),
  blocks   = "countimp imputes variable by variable; pass `visitSequence` instead",
  blots    = "pass method arguments directly, e.g. countimp(..., donors = 3)",
  calltype = "not applicable: countimp reads models from `formulas` or `predictorMatrix`")

.countimp_check_args <- function(arguments) {
  if (isFALSE(getOption("countimp.check.args", TRUE))) return(invisible(NULL))
  nms <- names(arguments)
  nms <- nms[nzchar(nms)]
  known <- .countimp_known_args()
  unknown <- setdiff(nms, known)
  if (!length(unknown)) return(invisible(NULL))

  ## R matches arguments before `...` partially, and callers rely on it
  ## (`print = FALSE` for `printFlag`). Accept an unambiguous abbreviation --
  ## the engine's own matching resolves it -- and only report what is left.
  unknown <- unknown[vapply(unknown, function(u) {
    hit <- startsWith(known, u)
    sum(hit) != 1L
  }, logical(1L))]
  if (!length(unknown)) return(invisible(NULL))

  mo <- intersect(unknown, names(.countimp_mice_only))
  if (length(mo))
    stop("Argument", if (length(mo) > 1L) "s" else "", " ",
         paste0("`", mo, "`", collapse = ", "),
         " belong", if (length(mo) == 1L) "s" else "",
         " to mice::mice() and ha", if (length(mo) == 1L) "s" else "ve",
         " no effect in countimp():\n",
         paste0("  `", mo, "`: ", .countimp_mice_only[mo], collapse = "\n"),
         call. = FALSE)

  ## Not a mice argument and not ours: most likely a typo. Name the nearest
  ## known argument, which turns the message into a fix.
  near <- vapply(unknown, function(u) {
    d <- utils::adist(u, known, ignore.case = TRUE)[1L, ]
    if (min(d) <= max(2L, nchar(u) %/% 3L)) known[which.min(d)] else NA_character_
  }, character(1L))
  hint <- ifelse(is.na(near), "", paste0(" (did you mean `", near, "`?)"))
  stop("Unknown argument", if (length(unknown) > 1L) "s" else "", " in countimp(): ",
       paste0("`", unknown, "`", hint, collapse = ", "),
       "\n  If this argument is meant for an imputation method, set ",
       "options(countimp.check.args = FALSE) to pass it through.",
       call. = FALSE)
}

## Refill of implausible draws in the outlier arm (EV = TRUE) ----------------
##
## When EV = TRUE, draws flagged as outliers (and NaNs) are set to NA and
## redrawn by predictive mean matching. The PMM variant used is midastouch,
## which lives in mice -- an optional dependency here. Routing the call
## through .countimp_find_method() keeps it working when mice is installed
## but not attached, and produces a message naming the cause when mice is
## absent, instead of "could not find function mice.impute.midastouch".
##
## `wy` must name the positions to be refilled, i.e. exactly the flagged
## draws. It cannot be left implicit: midastouch defaults to wy = !ry, and
## !ry is not the flagged set whenever the imputation targets differ from the
## missing-data pattern (mice's `where` argument) or whenever only a subset of
## the draws was flagged. Getting a vector of the wrong length back here used
## to be recycled silently into the caller's imputation vector.

.countimp_ev_refill <- function(y, ry, x, wy = NULL) {
  f <- .countimp_find_method("midastouch")
  if (is.null(f))
    stop("EV = TRUE redraws outlying imputations by predictive mean matching\n",
         "  (mice.impute.midastouch), which is provided by package 'mice'.\n",
         "  Either install mice (install.packages(\"mice\")) or call the\n",
         "  imputation method with EV = FALSE.", call. = FALSE)
  out <- f(y, ry, x, wy = wy)
  n.want <- if (is.null(wy)) sum(!ry) else sum(wy)
  if (length(out) != n.want)
    stop(sprintf(paste0("countimp: the donor method returned %d value(s) for %d ",
                        "position(s) to refill.\n  This is a bug in countimp -- ",
                        "please report it. Use EV = FALSE to avoid the outlier arm."),
                 length(out), n.want), call. = FALSE)
  out
}

## The whole EV = TRUE arm, in one place ---------------------------------------
##
## This block was copied into all 56 imputation methods: flag outliers, set
## them NA, write the retained draws into y[wy], build the flagged-position
## mask, refill. Copied code is where the two 2.6.0 defects came from, and the
## block has a subtlety that must not be re-derived per copy: `wy.ev` indexes
## into y, so it is which(wy)[idx], NOT idx -- using idx directly refills the
## wrong rows whenever the imputation targets are not the first positions.
##
## EV = TRUE remains deprecated and not recommended; see .countimp_warn_ev().
.countimp_ev_handle <- function(imp, y, x, wy) {
  outl <- getOutliers(imp, rho = c(0.3, 0.3), FLim = c(0.15, 0.85))
  idx  <- c(outl$iLeft, outl$iRight, which(is.nan(imp)))
  if (!length(idx)) return(imp)
  imp[idx] <- NA
  y[wy] <- imp
  wy.ev <- logical(length(y))
  wy.ev[which(wy)[idx]] <- TRUE
  imp[idx] <- .countimp_ev_refill(y, !is.na(y), x, wy = wy.ev)
  imp
}

## Which of these method names cannot be resolved? Used by check.method() to
## report all missing methods at once instead of failing on the first.
.countimp_missing_methods <- function(methods) {
  methods <- unique(methods[nzchar(methods)])
  methods[vapply(methods,
                 function(m) is.null(.countimp_find_method(m)),
                 logical(1))]
}

## ---------------------------------------------------------------------------
## The univariate-method contract (since 2.7.0)
##
## countimp's own FCS engine calls imputation methods that may come from
## countimp, from the user's workspace, or from mice. That makes the calling
## convention an INTERFACE, and it needs to be written down in one place
## instead of being implicit in the sampler loop.
##
## The convention is mice's, because that is what third-party methods are
## written against:
##
##   f(y, ry, x, wy = NULL, type = NULL, ...)
##     y     numeric/factor vector, the target, length n
##     ry    logical, TRUE where y is observed
##     x     NUMERIC MATRIX of predictors, n rows, no intercept column,
##           factors already expanded to dummies
##     wy    logical, TRUE where an imputation is wanted
##     type  integer vector, one entry per COLUMN OF x
##   returns a vector of length sum(wy)
##
## Why this wrapper exists -- a real incompatibility it fixes:
##
## The sampler inherited from mice 2.46 passed `x` as a data.frame and passed
## y/ry/x positionally. mice >= 3.0 passes a matrix. Most methods do not care;
## mice.impute.lasso.norm() does, because it hands x to glmnet, which requires
## a matrix. So that method worked when called from mice and failed when
## called from countimp -- with an error message ("'list' object cannot be
## coerced to type 'double'") that names neither the method nor the cause.
##
## The deeper point: passing a data.frame worked for 33 of mice's 34 methods
## BY LUCK, not by contract. Any future method that requires a matrix would
## have broken the same way. Converting once, here, removes a whole class of
## breakage instead of one instance of it -- which is the point of having a
## contract at all.
##
## Passing y/ry/x by name rather than by position closes the second half of
## the same gap: a method whose first arguments are ordered differently (or
## that takes only some of them) is called correctly.
## ---------------------------------------------------------------------------

## Resolve one per-variable argument out of its named-list form.
##
## Shared by `bounds` and `censor`; a third scale argument gets one line rather
## than a third copy of the same twenty. Returns `dots` unchanged when the
## argument is absent or not given per variable.
## `braucht` says whether the METHOD about to be called declares this argument.
## It has to, because a call may mix scale families: imputing one variable with
## `cp` (censor) and another with `bp` (bounds) passes BOTH arguments to both
## methods, and without this the bounds list -- correctly naming only the
## bounded variable -- was reported as incomplete while `cp` was being called
## ("'bounds' does not cover variable 'y', which is imputed by a bounded
## method", on a `y` that was imputed by a censored one). Measured on a
## two-variable example, 22 August 2026: the run aborted although both
## specifications were complete.
##
## An argument the method does not declare belongs to another variable's
## method, so it is dropped here rather than checked. The B58 rule -- a named
## list that omits the current variable is an error, not a silent fallback --
## keeps applying to the method that actually reads the argument.
.countimp_pro_variable <- function(dots, arg, vname, beispiel, wofuer,
                                   braucht = TRUE) {
  if (!(arg %in% names(dots)) || !is.list(dots[[arg]])) return(dots)
  if (!isTRUE(braucht)) { dots[[arg]] <- NULL; return(dots) }
  bl <- dots[[arg]]
  if (is.null(names(bl)) || any(!nzchar(names(bl))))
    stop("when '", arg, "' is a list, every element must be named after a ",
         "variable, ", beispiel, ".", call. = FALSE)
  if (is.null(vname))
    stop("Internal error: '", arg, "' given per variable but the variable name ",
         "did not reach the method.", call. = FALSE)
  if (!vname %in% names(bl))
    stop("'", arg, "' does not cover variable '", vname, "', which is imputed ",
         "by ", wofuer, ". Named entries: ", paste(names(bl), collapse = ", "),
         ". Add an entry for '", vname, "' or use a different method.",
         call. = FALSE)
  dots[[arg]] <- bl[[vname]]
  dots
}


.countimp_call_method <- function(f, y, ry, x, wy = NULL, type = NULL,
                                  vname = NULL, ...) {
  if (is.null(wy)) wy <- !ry

  ## `bounds` and `censor` describe the SCALE of one variable, so a single
  ## countimp() call may need a different value for each of them. Two spellings
  ## are accepted for both:
  ##
  ##   bounds = c(0, 8)                  one specification for every variable
  ##   bounds = list(sym = c(0, 8),      per variable, by name
  ##                 days = c(0, 5))
  ##
  ## The list form is resolved here, where the variable name is known, so the
  ## methods themselves only ever see one variable's value. A named list that
  ## does not mention the current variable is an error rather than a silent
  ## fallback to unbounded/uncensored imputation -- the user asked for the
  ## restriction and would not learn that one variable was left out (B58).
  ##
  ## `censor` shares the rule but not the shape: it is a single limit or one
  ## limit per case, so is.list() is what separates "per variable" from "per
  ## case" -- a per-case vector is numeric, never a list.
  dots <- list(...)
  fmls_f <- names(formals(f))
  dots <- .countimp_pro_variable(dots, "bounds", vname,
                                 "e.g. bounds = list(y = c(0, 8))",
                                 "a bounded method",
                                 braucht = "bounds" %in% fmls_f)
  dots <- .countimp_pro_variable(dots, "censor", vname,
                                 "e.g. censor = list(y = 10)",
                                 "a censored method",
                                 braucht = "censor" %in% fmls_f)

  ## Predictors as a numeric matrix -- the mice >= 3.0 convention.
  ## data.matrix() is the right tool: it keeps a matrix untouched and encodes
  ## the columns of a data.frame numerically. Factors should already be
  ## expanded by padModel(), but a residual one would silently become its
  ## integer codes, so it is caught rather than mis-modelled.
  if (!is.matrix(x)) {
    fac <- vapply(as.data.frame(x), function(z) is.factor(z) || is.character(z),
                  logical(1))
    if (any(fac))
      stop("Predictor(s) ", paste(names(fac)[fac], collapse = ", "),
           " reached the imputation method unexpanded.\n",
           "  Imputation methods require numeric predictors; a factor would\n",
           "  otherwise be modelled as its integer codes.", call. = FALSE)
    x <- data.matrix(x)
  }
  if (is.null(colnames(x)))
    colnames(x) <- paste0("V", seq_len(ncol(x)))

  ## `type` must align with the columns of x, or a method that reads it (the
  ## two-level methods do) silently uses the wrong roles.
  if (!is.null(type) && length(type) != ncol(x))
    stop("Internal error: `type` has ", length(type), " entries but `x` has ",
         ncol(x), " columns.", call. = FALSE)

  ## Only pass `type` to methods that can take it. Methods without `type` and
  ## without `...` would fail on an unused argument; and a method that passes
  ## `...` on to a fitting function (glmnet, ranger) would receive `type` as a
  ## stray argument there.
  fmls <- names(formals(f))
  args <- list(y = y, ry = ry, x = x, wy = wy)
  if (!is.null(type) && ("type" %in% fmls || "..." %in% fmls))
    args$type <- type
  if (length(dots) && !("..." %in% fmls))
    dots <- dots[names(dots) %in% fmls]
  do.call(f, c(args, dots))
}
