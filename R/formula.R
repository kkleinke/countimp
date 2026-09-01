## formula.R -- construction and validation of the random-effects term
##
## Background. The `.noint` methods drop the intercept from the RANDOM part
## of the model (it stays in the fixed part -- see the help pages). The code
## therefore assembles a term of the form
##
##     (0 + <random slopes> | <class variable>)
##
## When the user supplies no random slope at all (no predictor coded 2 or 6
## in `type`), the slope list is empty and the term collapses to
##
##     (0 | <class variable>)
##
## which asks for a random-effects block with zero parameters. glmmTMB does
## not reject this: it passes the empty structure down to TMB and the R
## session dies with a SEGFAULT inside compiled code -- no error, no
## traceback, nothing try() can catch. (lme4 raises a plain error here:
## "length(theta) > 0L is not TRUE".)
##
## Until glmmTMB 1.1.14 this crash was masked in some installations by an
## unrelated TMB ABI mismatch that aborted earlier. It is a real defect in
## countimp: a documented method combination kills the session.
##
## `.countimp_randeff()` is the single place where the term is built. It
## refuses the degenerate case with a message that names the offending
## method and tells the user what to change.

## ---------------------------------------------------------------------------
## PART 2 -- the lme4-style formula interface (since 2.6.0)
##
## Motivation. Up to 2.5.0 an imputation model was specified through three
## coupled objects: a `predictorMatrix` row, a `type` vector whose codes 1-6
## encode two orthogonal facts at once (fixed vs. fixed+random, and zero vs.
## count vs. both part), and a method NAME that encodes the distribution, the
## multilevel structure, the random-intercept handling and the draw method
## (63 exported names). Nothing in that scheme is visible at the call site:
## `type = c(2, 5, -2)` with method "2l.hnb.noint.zero" is a model definition
## written in three places at once.
##
## The formula interface writes the same model as one expression per model
## part, in the notation users already know from lme4/glmmTMB:
##
##   countimp(data = d,
##            formulas = list(y ~ x1 + x2 + (1 + x1 | id)),
##            zero     = list(y ~ x1 + (1 | id)),
##            family   = hurdle_nb(), draw = "bayes")
##
## DESIGN DECISIONS (settled with the maintainer, see NEWS 2.6.0)
##
##   1. Two arguments for the two model parts (`formulas` and `zero`) rather
##      than a pscl-style `y ~ count | zero` pipe. `|` already means "grouping"
##      inside lme4 terms; overloading it forces a hand-written splitter that
##      has to reach inside parentheses, and a typo silently yields a different
##      model. glmmTMB (ziformula) and brms use separate arguments too.
##   2. `family` selects the distribution. The `2l.` prefix and the
##      `.noint.zero/.count/.both` suffixes become properties of the FORMULA
##      -- (1|g) present or not, (0+x|g) vs (1+x|g) -- so they disappear from
##      the name. 58 of the 81 method names collapse to 7 families x `draw`
##      (counted, not estimated: test-B53 generates every name the constructors
##      can produce and checks that each one exists as an export). The
##      remaining 23 are the families that need a scale argument -- bounds,
##      censor -- which a formula has no place for.
##   3. `.` on the right-hand side expands to all other variables in `data`.
##      Incomplete variables with no formula keep the default method for their
##      measurement level; countimp reports them at startup.
##   4. `draw = "bayes" | "boot"` replaces the `.boot` name suffix.
##
## This is a COMPILER, not a second engine: it emits exactly the `type` vector,
## `predictorMatrix` row, method name and intercept flags that the existing
## methods already consume. The old interface is untouched and stays valid.
## ---------------------------------------------------------------------------

## Family constructors ------------------------------------------------------
## Each returns a small list; `.countimp_method_name()` maps it plus the
## formula properties to one of the existing method names (58 of the 81 are
## reachable this way; the count is checked in test-B53, not estimated here).

## `at_levels` says at which levels the family EXISTS as an exported method, and
## `args` carries a scale argument that belongs to the family rather than to
## the formula (a bound, a censoring limit).
##
## Both were added for the scale families. `at_levels` closes a trap that was
## measured, not imagined: .countimp_method_name() below composes "2l." with
## the stem, so a one-level-only family in a two-level formula produced the
## name "2l.cp" -- which does not exist. The user learned that from
## "Imputation method '2l.cp' was not found" in the middle of their run,
## instead of from the translation.
.countimp_fam <- function(stem, twopart, label, at_levels = c("1l", "2l"),
                          args = list()) {
  structure(list(stem = stem, twopart = twopart, label = label,
                 at_levels = at_levels, args = args),
            class = "countimp_family_spec")
}


## A scale argument is MANDATORY in the family that declares it.
##
## The alternative -- allowing bounded_poisson() and taking the bound from
## countimp(bounds = ...) -- would leave the scale in two places, which is the
## coupling this interface exists to remove. Whoever needs the number to come
## from elsewhere has the classic interface, unchanged.
.countimp_skala_pflicht <- function(fehlt, arg, fname, beispiel) {
  if (!fehlt) return(invisible(TRUE))
  stop("`", fname, "()` needs its `", arg, "` argument: ", beispiel, ".\n",
       "  The scale belongs to the family, so that it stands next to the ",
       "variable it\n  describes. To supply it separately instead, use the ",
       "classic interface with\n  `method` and `", arg, "`.", call. = FALSE)
}

#' Distribution families for the countimp formula interface
#'
#' @description Select the imputation model for a count variable. Used as the
#'   \code{family} argument of \code{\link{countimp}} when models are specified
#'   through formulas. The multilevel structure and the random-intercept
#'   handling are taken from the formula, not from the family.
#'
#' @details
#' \strong{Scale families.} \code{zerotrunc_*}, \code{bounded_*} and
#' \code{censored_*} describe counts whose \emph{support} is restricted, and
#' the restriction belongs to the family rather than to the formula:
#'
#' \preformatted{
#' countimp(d,
#'          formulas = list(visits ~ x1, days ~ x1),
#'          family   = list(visits = censored_poisson(censor = 10),
#'                          days   = bounded_poisson(bounds = c(0, 5))),
#'          m = 5)
#' }
#'
#' \strong{Not a family: quasipoisson.} \code{quasipoisson} appears in the
#' method names (\code{method = c(y = "quasipoisson")}) but deliberately has
#' \emph{no} family constructor. It is not a distribution but a variance
#' assumption on a Poisson core: quasi-likelihood has no likelihood, and
#' therefore no posterior for \code{draw = "bayes"} to draw from. The method
#' remains available through the classic interface for whoever wants it; the
#' formula interface does not lead anyone there. Writing
#' \code{family = quasipoisson()} -- the spelling from \code{glm()} -- reports
#' this rather than failing obscurely.
#'
#' The scale then stands next to the variable it describes, instead of in a
#' second list that has to be kept in step with the method names. Three
#' consequences, all deliberate:
#' \itemize{
#'   \item \code{bounds} and \code{censor} are \strong{mandatory} in the
#'     families that take them. Supplying the number separately is still
#'     possible through the classic interface
#'     (\code{method = c(y = "bp")} plus \code{bounds = }), but not by mixing
#'     the two.
#'   \item Giving the same scale \strong{both} ways is an error, not a
#'     precedence rule: the two can differ, and imputations drawn under a bound
#'     that was not meant still look like plausible counts.
#'   \item These families exist \strong{only for single-level data}. There is
#'     no \code{2l.bp} or \code{2l.cp}; a grouping term in the formula is
#'     refused by name rather than compiled into a method that does not exist.
#' }
#'
#' @param bounds Numeric \code{c(lower, upper)} for the bounded families; see
#'   \code{\link{mice.impute.bp}}.
#' @param censor The censoring limit for the censored families -- one number,
#'   or one per case; see \code{\link{mice.impute.cp}}.
#'
#' @return An object of class \code{countimp_family_spec}.
#' @seealso \code{\link{countimp_spec}} to inspect the translation, including
#'   where the scale argument ended up (\code{$args}).
#' @examples
#' hurdle_nb()
#' zi_poisson()
#' compois()
#' bounded_poisson(bounds = c(0, 8))
#' censored_poisson(censor = 10)
#' zerotrunc_nb()
#' @name countimp_families
#' @export
hurdle_nb      <- function() .countimp_fam("hnb",  TRUE,  "hurdle negative binomial")
#' @rdname countimp_families
#' @export
hurdle_poisson <- function() .countimp_fam("hp",   TRUE,  "hurdle Poisson")
#' @rdname countimp_families
#' @export
zi_nb          <- function() .countimp_fam("zinb", TRUE,  "zero-inflated negative binomial")
#' @rdname countimp_families
#' @export
zi_poisson     <- function() .countimp_fam("zip",  TRUE,  "zero-inflated Poisson")
#' @rdname countimp_families
#' @export
nb             <- function() .countimp_fam("nb",   FALSE, "negative binomial")
#' @rdname countimp_families
#' @export
poisson_count  <- function() .countimp_fam("poisson", FALSE, "Poisson")
## compois is the one later family that fits this interface without an extra
## argument, and every name it can generate exists: cmp, cmp.boot, 2l.cmp,
## 2l.cmp.noint and their .boot variants. The bounded, zero-truncated and
## censored families are deliberately NOT here -- each needs a scale argument
## (`bounds`, `censor`) that the formula carries no place for, and two of them
## have no two-level counterpart, so a family constructor could produce a name
## like "2l.cp" that does not exist. test-B53 checks that property for every
## constructor rather than trusting this note.
#' @rdname countimp_families
#' @export
compois       <- function() .countimp_fam("cmp", FALSE, "Conway-Maxwell-Poisson")

## The scale families. Each exists at ONE level only -- there is no 2l.ztp,
## 2l.bp or 2l.cp -- which is what the `at_levels` field records.
#' @rdname countimp_families
#' @export
zerotrunc_poisson <- function()
  .countimp_fam("ztp", FALSE, "zero-truncated Poisson", "1l")
#' @rdname countimp_families
#' @export
zerotrunc_nb      <- function()
  .countimp_fam("ztnb", FALSE, "zero-truncated negative binomial", "1l")
#' @rdname countimp_families
#' @export
bounded_poisson   <- function(bounds) {
  .countimp_skala_pflicht(missing(bounds), "bounds", "bounded_poisson",
                          "bounded_poisson(bounds = c(0, 8))")
  .countimp_fam("bp", FALSE, "bounded Poisson", "1l", list(bounds = bounds))
}
#' @rdname countimp_families
#' @export
bounded_nb        <- function(bounds) {
  .countimp_skala_pflicht(missing(bounds), "bounds", "bounded_nb",
                          "bounded_nb(bounds = c(0, 8))")
  .countimp_fam("bnb", FALSE, "bounded negative binomial", "1l",
                list(bounds = bounds))
}
#' @rdname countimp_families
#' @export
censored_poisson  <- function(censor) {
  .countimp_skala_pflicht(missing(censor), "censor", "censored_poisson",
                          "censored_poisson(censor = 10)")
  .countimp_fam("cp", FALSE, "right-censored Poisson", "1l",
                list(censor = censor))
}
#' @rdname countimp_families
#' @export
censored_nb       <- function(censor) {
  .countimp_skala_pflicht(missing(censor), "censor", "censored_nb",
                          "censored_nb(censor = 10)")
  .countimp_fam("cnb", FALSE, "right-censored negative binomial", "1l",
                list(censor = censor))
}

#' @export
print.countimp_family_spec <- function(x, ...) {
  cat("<countimp family: ", x$label, ">\n", sep = "")
  invisible(x)
}

## Formula dissection --------------------------------------------------------
## Splits an lme4 formula into fixed terms, random slopes, grouping factor and
## the random-intercept flag. Deliberately hand-rolled on terms() rather than
## calling lme4::findbars(): countimp must not gain a dependency on lme4 just
## to read a formula, and the bar syntax we accept is a small subset.


## Variables appearing inside `( ... | g )` terms -- excluded from `.` expansion
## so that a random slope is not silently duplicated as a fixed-only term.
.countimp_bar_vars <- function(f) {
  tl <- attr(stats::terms(f, keep.order = TRUE, allowDotAsName = TRUE),
             "term.labels")
  bars <- tl[grepl("\\|", tl)]
  if (!length(bars)) return(character(0))
  unique(unlist(lapply(bars, function(b) all.vars(stats::as.formula(paste("~", b))))))
}

.countimp_parse_formula <- function(f, data, part = "count") {
  if (!inherits(f, "formula"))
    stop("`", part, "` must be a formula, got ", class(f)[1L], call. = FALSE)

  lhs <- if (length(f) == 3L) all.vars(f[[2L]]) else character(0)
  if (part == "count" && length(lhs) != 1L)
    stop("Formula must have exactly one variable on the left-hand side: ",
         deparse1(f), call. = FALSE)

  ## `.` must be resolved BEFORE terms(): terms.formula() rejects a dot unless
  ## it is given `data`, and passing `data` there would also expand it inside
  ## the random-effects bars, which is not what `(1 | id)` means.
  ftxt <- paste(deparse(f), collapse = " ")
  if (grepl("(^|[^[:alnum:]._])\\.([^[:alnum:]._]|$)", ftxt)) {
    others <- setdiff(names(data), c(lhs, .countimp_bar_vars(f)))
    repl <- if (length(others)) paste(others, collapse = " + ") else "1"
    ftxt <- sub("(^|[^[:alnum:]._])\\.([^[:alnum:]._]|$)",
                paste0("\\1", repl, "\\2"), ftxt)
    f <- stats::as.formula(ftxt, env = environment(f))
  }

  tl <- attr(stats::terms(f, keep.order = TRUE), "term.labels")
  bars <- grepl("\\|", tl)
  fixed <- tl[!bars]

  ## Several `( ... | g )` terms are allowed since the three-level extension:
  ## (1 | school) + (1 | class) is a three-level model, and glmmTMB fits it
  ## natively. `re` carries one entry per level in the order written; `group`,
  ## `slopes` and `rint` keep their old single-level meaning for the FIRST
  ## level, because that is what `type` and the method name can encode.
  ##
  ## The limit is deliberate and named rather than silently absorbed: the
  ## classic interface transports the model as a predictorMatrix row (codes
  ## -2, 0..6) plus a method name. That channel can say "this variable groups"
  ## (-2, now more than once) and "this variable has a random slope" (2/4/6),
  ## but it cannot say WHICH level a slope belongs to, and it has only one
  ## .noint flag. So levels after the first take `(1 | g)` and nothing else.
  grp <- character(0); slopes <- character(0); rint <- FALSE
  re <- list()
  if (any(bars)) {
    for (b in tl[bars]) {
      parts <- strsplit(b, "\\|")[[1L]]
      g <- trimws(parts[2L])
      toks <- trimws(strsplit(trimws(parts[1L]), "\\+")[[1L]])
      ri <- !("0" %in% toks) && !("-1" %in% toks)
      sl <- setdiff(toks, c("0", "1", "-1", ""))
      re[[length(re) + 1L]] <- list(group = g, slopes = sl, rint = ri)
    }
    if (anyDuplicated(vapply(re, `[[`, character(1), "group")))
      stop("The same grouping factor appears in more than one random-effect ",
           "term: ", deparse1(f), ".\n  Write one term per grouping factor, ",
           "e.g. (1 + x | id) instead of (1 | id) + (0 + x | id).",
           call. = FALSE)
    for (k in seq_along(re)[-1L]) {
      if (length(re[[k]]$slopes) || !re[[k]]$rint)
        stop("Only the first random-effect term may carry random slopes or ",
             "drop its intercept.\n  Term ", k, " is (",
             paste(c(if (!re[[k]]$rint) "0", re[[k]]$slopes), collapse = " + "),
             " | ", re[[k]]$group, ") in: ", deparse1(f),
             ".\n  Reason: countimp passes the model on as type codes plus a ",
             "method name, and that\n  channel carries one slope set and one ",
             "intercept flag. Put the slopes on the\n  first term, or impute ",
             "this variable with the classic interface.", call. = FALSE)
    }
    grp    <- vapply(re, `[[`, character(1), "group")
    slopes <- re[[1L]]$slopes
    rint   <- re[[1L]]$rint
    ## a random slope is a fixed effect too (lme4 semantics)
    fixed <- unique(c(fixed, slopes))
  }

  ## drop the grouping factors from the fixed part if the user listed them
  if (length(grp)) fixed <- setdiff(fixed, grp)

  unknown <- setdiff(c(fixed, slopes, grp), names(data))
  if (length(unknown))
    stop("Variable(s) not in `data`: ", paste(unknown, collapse = ", "),
         call. = FALSE)

  list(y = if (length(lhs)) lhs else NA_character_,
       fixed = fixed, slopes = slopes, group = grp,
       rint = rint, re = re, twolevel = length(grp) > 0L)
}

## type / predictorMatrix construction ---------------------------------------
## Codes (see ?mice.impute.2l.hnb):
##   1 fixed both parts        2 fixed+random both parts
##   3 fixed count only        4 fixed+random count only
##   5 fixed zero only         6 fixed+random zero only
##  -2 class (grouping) variable

.countimp_build_type <- function(cnt, zro, varnames) {
  type <- setNames(integer(length(varnames)), varnames)

  in.c  <- function(v) v %in% cnt$fixed
  in.z  <- function(v) !is.null(zro) && v %in% zro$fixed
  ran.c <- function(v) v %in% cnt$slopes
  ran.z <- function(v) !is.null(zro) && v %in% zro$slopes

  for (v in varnames) {
    if (identical(v, cnt$y)) next
    both <- in.c(v) && (is.null(zro) || in.z(v))
    only.c <- in.c(v) && !in.z(v)
    only.z <- !in.c(v) && in.z(v)
    rand <- ran.c(v) || ran.z(v)
    type[v] <- if (both)   { if (rand) 2L else 1L }
          else if (only.c) { if (rand) 4L else 3L }
          else if (only.z) { if (rand) 6L else 5L }
          else 0L
  }
  if (length(cnt$group)) type[cnt$group] <- -2L
  type
}

## Method-name resolution ----------------------------------------------------
## Maps family + formula properties + draw onto one of the exported names.
## Returns the name AND the intercept flags, because for the .noint variants
## the flags are baked into the name (they are formals with fixed defaults).

.countimp_method_name <- function(fam, cnt, zro, draw = c("bayes", "boot")) {
  draw <- match.arg(draw)
  stem <- fam$stem
  twol <- cnt$twolevel

  ## Refuse a level the family does not have, instead of composing a name that
  ## does not exist. Families built before this field carry both levels.
  at_levels <- if (is.null(fam$at_levels)) c("1l", "2l") else fam$at_levels
  if (isTRUE(twol) && !("2l" %in% at_levels))
    stop("The ", fam$label, " family has no two-level counterpart, but the ",
         "formula carries a\n  grouping term. There is no method \"2l.",
         stem, "\". Drop the (1 | g) term, or use a\n  family that exists at ",
         "two levels.", call. = FALSE)
  if (!isTRUE(twol) && !("1l" %in% at_levels))
    stop("The ", fam$label, " family exists only for two-level models, but the ",
         "formula has no\n  grouping term. Add one, e.g. + (1 | id).",
         call. = FALSE)

  if (!twol) {
    ## single-level: no random intercepts, no .noint variants
    ## NOTE: the engine prepends "mice.impute." itself (see zisampler.R:82),
    ## so the bare stem is what belongs in `method`.
    nm <- if (stem == "poisson") "poisson" else stem
    if (draw == "boot") nm <- paste0(nm, ".boot")
    return(list(method = nm, intercept.c = TRUE, intercept.z = TRUE))
  }

  ic <- cnt$rint
  iz <- if (fam$twopart && !is.null(zro)) zro$rint else TRUE

  suffix <- if (fam$twopart) {
    if      ( ic &&  iz) ""
    else if (!ic &&  iz) ".noint.count"
    else if ( ic && !iz) ".noint.zero"
    else                 ".noint.both"
  } else {
    if (ic) "" else ".noint"
  }
  stem2l <- if (stem == "poisson") "poisson" else stem
  nm <- paste0("2l.", stem2l, suffix)
  if (draw == "boot") nm <- paste0(nm, ".boot")
  list(method = nm, intercept.c = ic, intercept.z = iz)
}

## Scale arguments of the families, as one named list per argument:
##   list(bounds = list(days = c(0, 5)), censor = list(visits = 10))
## which is exactly the per-variable form countimp()'s methods read.
.countimp_skalen_args <- function(fam.list) {
  out <- list()
  for (v in names(fam.list)) {
    a <- fam.list[[v]]$args
    if (!length(a)) next
    for (nm in names(a)) out[[nm]][[v]] <- a[[nm]]
  }
  out
}


#' Compile a formula specification into countimp's classic arguments
#'
#' @description Internal workhorse of the formula interface. Translates
#'   \code{formulas}/\code{zero}/\code{family}/\code{draw} into the
#'   \code{method} vector, \code{predictorMatrix} and per-variable \code{type}
#'   vectors that the imputation methods consume. Exposed so that users can
#'   inspect -- and report -- what was actually fitted.
#' @param data data frame to be imputed.
#' @param formulas list of count-part formulas (or a single formula). Random
#'   effects are written lme4-style. Since 3.0.0 a formula may carry SEVERAL
#'   grouping terms, which is how a three-level model is specified:
#'   \code{y ~ x1 + (1 | school) + (1 | class)}. The first term carries the
#'   random slopes and decides the \code{.noint} variant; further terms enter
#'   as plain random intercepts, and anything else in them is an error rather
#'   than a silent omission (the model travels on as type codes, which carry
#'   one slope set and one intercept flag).
#' @param zero zero-part formula(s), or \code{NULL}. Which target a zero-part
#'   formula belongs to may be stated in three ways: as a one-sided formula
#'   (\code{zero = ~ z1 + (1 | id)}) when \code{formulas} names a single
#'   target; as a two-sided formula (\code{zero = y ~ z1 + (1 | id)}); or as a
#'   named list (\code{zero = list(y = ~ z1 + (1 | id))}). With more than one
#'   target, a one-sided formula cannot be attributed and is an error rather
#'   than a guess. The random-intercept term of the zero part is read
#'   separately from that of the count part and selects the \code{.noint.zero}
#'   / \code{.noint.count} / \code{.noint.both} variant accordingly; omitting
#'   \code{zero} for a two-part family reuses the count-part formula.
#' @param family a family from \code{\link{countimp_families}}, or a named list
#'   mapping target variables to families.
#' @param draw \code{"bayes"} (default) or \code{"boot"}.
#' @return A list with \code{method}, \code{predictorMatrix}, \code{type} and
#'   \code{blots}, ready to be passed to \code{\link{countimp}}.
#' @examples
#' ## --- flat file (wide format, one row per person) ------------------------
#' data(crim4w)
#' spec <- countimp_spec(crim4w,
#'           formulas = list(ACRIM ~ BCRIM + FEMALE + RE),
#'           family = hurdle_poisson())
#' spec$method[spec$method != ""]
#'
#' ## --- multilevel (long format, several occasions per person) -------------
#' ## the term (1 | ID) requests a random intercept and selects the
#' ## corresponding two-level method
#' data(crim4l)
#' spec2 <- countimp_spec(crim4l,
#'            formulas = list(DELINQ ~ TIME + FEMALE + (1 | ID)),
#'            family = poisson_count())
#' spec2$method[spec2$method != ""]
#'
#' ## --- separate zero part --------------------------------------------------
#' ## the count part keeps its random intercept, the zero part has a random
#' ## slope only -- which selects the ".noint.zero" variant
#' spec3 <- countimp_spec(crim4l,
#'            formulas = list(DELINQ ~ TIME + FEMALE + (1 | ID)),
#'            zero     = ~ TIME + (0 + TIME | ID),
#'            family   = hurdle_poisson())
#' spec3$method[spec3$method != ""]
#' @export
countimp_spec <- function(data, formulas, zero = NULL,
                          family = hurdle_nb(), draw = "bayes") {
  if (inherits(formulas, "formula")) formulas <- list(formulas)
  if (inherits(zero, "formula"))     zero     <- list(zero)
  if (!is.list(formulas) || !length(formulas))
    stop("`formulas` must be a formula or a non-empty list of formulas.",
         call. = FALSE)

  varnames <- names(data)
  nvar <- length(varnames)
  method <- setNames(rep("", nvar), varnames)
  pred <- matrix(0L, nvar, nvar, dimnames = list(varnames, varnames))
  types <- list()

  ## family may be one spec for all, or a named list per target
  fam.of <- function(y) {
    if (inherits(family, "countimp_family_spec")) return(family)
    if (is.list(family) && !is.null(family[[y]])) return(family[[y]])
    stop("No family given for target variable '", y, "'.", call. = FALSE)
  }
  ## Zero-part formulas are matched to targets. Three notations are accepted:
  ##   (a) named list       zero = list(y = ~ x + (1 | id))
  ##   (b) two-sided        zero = y ~ x + (1 | id)
  ##   (c) one-sided        zero = ~ x + (1 | id)
  ## (c) is the notation the documentation and every example use, and it was
  ## silently DROPPED before: the old matcher only looked at the LHS, a
  ## one-sided formula has none, so `zro` stayed NULL and fell back to `zro <-
  ## cnt`. Consequence -- the zero part inherited the COUNT part's random
  ## structure and predictors: `zero = ~ z1 + (0 + z1 | g)` produced method
  ## "2l.hnb" instead of "2l.hnb.noint.zero", and the type vector described the
  ## count part twice. No error, no warning, a different model than requested.
  ziele <- vapply(formulas, function(f) {
    if (!inherits(f, "formula") || length(f) != 3L)
      stop("Each entry of `formulas` must be a two-sided formula ",
           "(y ~ predictors).", call. = FALSE)
    v <- all.vars(f[[2L]])
    if (length(v) != 1L)
      stop("Formula must have exactly one variable on the left-hand side: ",
           deparse(f), call. = FALSE)
    v
  }, character(1))

  zmap <- list()
  if (!is.null(zero)) {
    if (!is.list(zero) || !all(vapply(zero, inherits, logical(1), "formula")))
      stop("`zero` must be a formula or a list of formulas.", call. = FALSE)
    nz <- names(zero)
    if (is.null(nz)) nz <- rep("", length(zero))
    frei <- ziele
    for (i in seq_along(zero)) {
      z <- zero[[i]]
      y <- if (nzchar(nz[i])) nz[i]
           else if (length(z) == 3L) all.vars(z[[2L]])[1L]
           else if (length(frei) == 1L) frei
           else if (length(zero) == length(ziele)) ziele[i]
           else NA_character_
      if (is.na(y))
        stop("Cannot tell which target the one-sided `zero` formula ",
             deparse(z), " belongs to. With more than one entry in ",
             "`formulas`, name the zero part explicitly -- either as ",
             "`zero = list(", ziele[1L], " = ~ ...)` or as a two-sided ",
             "formula `", ziele[1L], " ~ ...`.", call. = FALSE)
      if (!y %in% ziele)
        stop("`zero` formula ", deparse(z), " refers to '", y,
             "', which has no count-part formula in `formulas`.",
             call. = FALSE)
      if (!is.null(zmap[[y]]))
        stop("More than one `zero` formula given for '", y, "'.",
             call. = FALSE)
      ## A one-sided zero part gets the target as LHS so that the parser --
      ## which insists on a response -- sees a well-formed formula. Built on
      ## the call, not via deparse/reformulate: deparse() breaks long formulas
      ## into several strings and reformulate() would then drop all but the
      ## first. environment() is carried over so that data-masking still works.
      if (length(z) == 3L) {
        zmap[[y]] <- z
      } else {
        zz <- call("~", as.name(y), z[[2L]])
        zz <- eval(call("as.formula", zz), envir = environment(z))
        environment(zz) <- environment(z)
        zmap[[y]] <- zz
      }
      frei <- setdiff(frei, y)
    }
  }
  zero.of <- function(y) zmap[[y]]

  for (f in formulas) {
    cnt <- .countimp_parse_formula(f, data, "count")
    y <- cnt$y
    if (!y %in% varnames)
      stop("Target variable '", y, "' is not a column of `data`.", call. = FALSE)
    fam <- fam.of(y)
    zf <- zero.of(y)
    zro <- if (!is.null(zf)) .countimp_parse_formula(zf, data, "zero") else NULL

    if (!fam$twopart && !is.null(zro))
      warning("Family '", fam$label, "' has no zero part; `zero` formula for '",
              y, "' is ignored.", call. = FALSE)

    ## a two-part family with no explicit zero formula: reuse the count part
    if (fam$twopart && is.null(zro)) zro <- cnt

    ## Both model parts must name the SAME grouping levels. The type codes
    ## carry -2 per variable, not per model part, so a zero formula with fewer
    ## (or other) grouping terms cannot be expressed -- and was silently given
    ## the count part's levels instead. Measured: `zero = ~ z1 + (1 | id)`
    ## beside `y ~ x1 + (1 | id) + (1 | kl)` produced the zero-part formula
    ## `Y ~ z1 + (1 | id) + (1 | kl)`, i.e. a random effect the user did not
    ## write. Invisible in the output, and exactly the shape of B53.
    ##
    ## Possible only since the three-level extension: with one grouping factor
    ## the two parts could not differ.
    if (fam$twopart && !setequal(cnt$group, zro$group))
      stop("The count and zero parts name different grouping levels: (",
           paste(cnt$group, collapse = ", "), ") against (",
           paste(zro$group, collapse = ", "), ") for '", y, "'.\n",
           "  countimp passes the grouping on as type codes, which are per ",
           "VARIABLE and not per\n  model part, so the two parts cannot ",
           "differ in their levels. Give both parts the\n  same (1 | g) ",
           "terms, or drop `zero` to let the zero part reuse the count part.",
           call. = FALSE)

    res <- .countimp_method_name(fam, cnt, zro, draw)
    method[y] <- res$method
    ty <- .countimp_build_type(cnt, if (fam$twopart) zro else NULL, varnames)
    types[[y]] <- ty
    pred[y, ] <- ifelse(ty != 0L, 1L, 0L)
    pred[y, y] <- 0L
    if (length(cnt$group)) pred[y, cnt$group] <- -2L
  }

  ## incomplete variables without a formula: default method, reported
  incomplete <- varnames[colSums(is.na(data)) > 0L]
  unmodelled <- setdiff(incomplete, names(types))
  if (length(unmodelled)) {
    for (v in unmodelled) {
      x <- data[[v]]
      method[v] <- if (is.numeric(x)) "pmm"
                   else if (is.logical(x) || (is.factor(x) && nlevels(x) == 2L)) "logreg"
                   else if (is.ordered(x)) "polr"
                   else "polyreg"
      pred[v, ] <- 1L; pred[v, v] <- 0L
    }
    message("countimp: no formula given for incomplete variable(s) ",
            paste(unmodelled, collapse = ", "),
            " -- using default method(s) ",
            paste(unique(method[unmodelled]), collapse = ", "), ".")
  }

  blots <- list()
  for (y in names(types)) {
    ty <- types[[y]]
    blots[[y]] <- list(type = ty[ty != 0L])
  }

  ## Scale arguments carried by the families, collected per variable -- the
  ## shape the methods already understand (see .countimp_pro_variable). Part of
  ## the spec so that countimp_spec() shows the WHOLE translation: a user who
  ## inspects it should not have to guess where the bound went.
  structure(list(method = method, predictorMatrix = pred,
                 type = types, blots = blots, draw = draw,
                 args = .countimp_skalen_args(
                   stats::setNames(lapply(names(types), fam.of), names(types)))),
            class = "countimp_spec")
}

#' @export
print.countimp_spec <- function(x, ...) {
  cat("<countimp specification>\n")
  m <- x$method[x$method != ""]
  for (v in names(m)) {
    cat("  ", v, ": ", m[[v]], sep = "")
    if (!is.null(x$type[[v]])) {
      ty <- x$type[[v]]; ty <- ty[ty != 0L]
      cat("  [type: ", paste(names(ty), ty, sep = "=", collapse = ", "), "]", sep = "")
    }
    cat("\n")
  }
  cat("  draw: ", x$draw, "\n", sep = "")
  invisible(x)
}

## `style` selects the dialect of the returned string:
##   "glmmTMB"  (1+x|grp)   -- a term inside a one-sided model formula
##   "nlme"     ~1+x|grp    -- a standalone formula for lme()/glmmPQL()
## Both dialects share the empty-random-part guard below, which is the reason
## this helper is used by mice.impute.2l.pmm() as well: nlme reacts to
## `~0 | grp` with "'data' must be of a vector type, was 'NULL'" from deep
## inside model.matrix.reStruct(), which names neither the intercept switch nor
## the missing slope (B53).
.countimp_randeff <- function(grp.name, slopes = character(0),
                              intercept = TRUE, part = "count",
                              style = c("glmmTMB", "nlme")) {
  style    <- match.arg(style)
  grp.name <- as.character(grp.name)
  if (!length(grp.name))
    stop("countimp: no grouping factor for the random part.", call. = FALSE)
  slopes   <- slopes[nzchar(slopes)]

  if (!intercept && length(slopes) == 0L) {
    stop(sprintf(paste0(
      "Random part of the %s model is empty.\n",
      "  This method removes the random intercept (intercept.%s = FALSE), so at\n",
      "  least one random slope is required, but `type` codes none.\n",
      "  Fix one of:\n",
      "    - code a predictor as 2 (fixed + random slope) or 6 in `type`, or\n",
      "    - use the method WITHOUT the '.noint' suffix, which keeps the random\n",
      "      intercept.\n",
      "  (Passing (0 | %s) to glmmTMB crashes the R session rather than\n",
      "  raising an error, so countimp stops here instead.)"),
      part, substr(part, 1L, 1L), grp.name[1L]), call. = FALSE)
  }

  ## First level: intercept flag and slopes as given. Further levels: plain
  ## random intercepts -- the only thing the type codes can express for them
  ## (.countimp_parse_formula refuses anything else with a named error).
  lead <- if (intercept) "1" else "0"
  kern <- if (length(slopes))
    paste0(lead, "+", paste(slopes, collapse = "+"), "|", grp.name[1L])
  else
    paste0(lead, "|", grp.name[1L])
  if (style == "nlme") {
    ## nlme takes one grouping level per call; the two-level hurdle path is the
    ## only caller and passes a single factor. Refuse the rest rather than
    ## silently dropping a level.
    if (length(grp.name) > 1L)
      stop("countimp: the nlme backend handles one grouping level, got ",
           length(grp.name), " (", paste(grp.name, collapse = ", "), ").",
           call. = FALSE)
    return(paste0("~", kern))
  }
  weitere <- if (length(grp.name) > 1L)
    paste0(" + (1|", grp.name[-1L], ")", collapse = "") else ""
  paste0("(", kern, ")", weitere)
}

## Decode `type` into the two model parts -------------------------------------
## The inverse of .countimp_build_type(): given the integer `type` vector that
## mice derives from the predictorMatrix, return the fixed and random predictor
## names for the count part and for the zero part, plus the grouping variable.
##
## This existed inline in 68 places across the two-level methods, once per
## model part per method. Two defects fixed in 2.6.0 both came from that
## duplication: intercept.z was applied to the count part, and the count-part
## fixed set was overwritten with nam[-grp]. A single decoder removes the
## opportunity to reintroduce either.
##
## Codes: 1 fixed both | 2 fixed+random both | 3 fixed count | 4 fixed+random
## count | 5 fixed zero | 6 fixed+random zero | -2 class variable.
.countimp_decode_type <- function(type, varnames) {
  type <- as.integer(type)
  if (length(type) != length(varnames))
    stop("countimp: length(type) != length(varnames) -- ", length(type),
         " vs ", length(varnames), ".", call. = FALSE)

  ## More than one -2 means more than one grouping level: (1 | school) +
  ## (1 | class) is coded as two -2 entries, in column order. Up to 3.0.0 a
  ## second -2 was an error ("only one class allowed!"); it is now the
  ## three-level case. Random slopes stay with the FIRST grouping variable --
  ## the type codes carry no level for them (see .countimp_parse_formula).
  grp <- which(type == -2L)
  if (!length(grp))
    stop("countimp: no class variable in `type`.\n",
         "  A multilevel method needs at least one predictor coded -2 (the\n",
         "  grouping variable). Code the cluster identifier as -2 in the\n",
         "  predictorMatrix row for this target.", call. = FALSE)

  pick <- function(codes) varnames[sort(unique(which(type %in% codes)))]

  list(
    ## a character VECTOR since the three-level extension; group[1] is the
    ## level that carries the random slopes and the .noint flag
    group   = varnames[grp],
    ## count part: shared (1, 2) plus count-only (3, 4)
    c.fixed = pick(c(1L, 2L, 3L, 4L)),
    c.slope = pick(c(2L, 4L)),
    ## zero part: shared (1, 2) plus zero-only (5, 6)
    z.fixed = pick(c(1L, 2L, 5L, 6L)),
    z.slope = pick(c(2L, 6L))
  )
}

## Build the two glmmTMB formulas for a two-level two-part method.
## `part` selects which side; the random part goes through .countimp_randeff(),
## which raises an informative error instead of letting (0 | g) segfault.
.countimp_2l_formula <- function(dec, part = c("count", "zero"),
                                 response, intercept = TRUE) {
  part <- match.arg(part)
  fixed  <- if (identical(part, "count")) dec$c.fixed else dec$z.fixed
  slopes <- if (identical(part, "count")) dec$c.slope else dec$z.slope
  fx <- if (length(fixed)) paste(fixed, collapse = "+") else "1"
  re <- .countimp_randeff(dec$group, slopes, intercept = intercept, part = part)
  stats::as.formula(paste(response, "~", fx, "+", re))
}
