#' Multiple imputation of incomplete count data
#'
#' @description Main entry point of the package. Imputation models can be
#'   specified in two equivalent ways:
#'
#'   \describe{
#'     \item{Formula interface (recommended, since 2.6.0)}{Pass \code{formulas}
#'       and \code{family}. Each model is written in \pkg{lme4} notation, e.g.
#'       \code{y ~ x1 + x2 + (1 | id)}, and the family selects the count
#'       distribution. Predictor roles and the random-intercept structure are
#'       read from the formula; no \code{predictorMatrix} is needed.}
#'     \item{Classic interface}{Pass \code{method} and \code{predictorMatrix}
#'       with \pkg{mice}'s type codes (1 = fixed effect, 2 = random effect,
#'       3 = class variable, -2 = grouping). Retained for backwards
#'       compatibility.}
#'   }
#'
#'   The two interfaces are mutually exclusive. Supplying \code{formulas}
#'   together with \code{method} or \code{predictorMatrix} is an error.
#'
#' @param data A data frame or matrix containing the incomplete data. May be
#'   given positionally, in both the formula and the classic interface.
#' @param ... Arguments passed on to the imputation engine. In addition to the
#'   usual arguments (\code{m}, \code{maxit}, \code{seed}, \code{method},
#'   \code{predictorMatrix}, ...), the following control the formula
#'   interface:
#'   \describe{
#'     \item{\code{formulas}}{A formula or named list of formulas, one per
#'       variable to be imputed, in \pkg{lme4} notation. A term
#'       \code{(1 | id)} requests a random intercept for the cluster variable
#'       \code{id} and selects the corresponding multilevel method.}
#'     \item{\code{family}}{A single family specification or a named list with
#'       one entry per formula. See \code{\link{countimp_families}} for the
#'       available choices: \code{\link{poisson_count}}, \code{\link{nb}},
#'       \code{\link{zi_poisson}}, \code{\link{zi_nb}},
#'       \code{\link{hurdle_poisson}} and \code{\link{hurdle_nb}}.}
#'     \item{\code{zero}}{The second model part of a two-part family --
#'       the inflation part of \code{zi_poisson()}/\code{zi_nb()}, the hurdle
#'       of \code{hurdle_poisson()}/\code{hurdle_nb()}. Omitting it reuses the
#'       count part's right-hand side. Three notations are accepted:
#'       \code{zero = ~ z1} (one-sided, when a single variable is imputed),
#'       \code{zero = y ~ z1} (two-sided), or
#'       \code{zero = list(y = ~ z1)} (named, and the only one that works for
#'       several targets). See the Zero-inflated and hurdle models section.}
#'     \item{\code{draw}}{How regression parameters are drawn:
#'       \code{"bayes"} (default) draws from the approximate posterior,
#'       \code{"boot"} uses a bootstrap of the observed data.}
#'   }
#' @section Variables that live at a cluster level:
#' In multilevel data some variables do not vary within their cluster: class
#' size, school type, ward capacity. They occupy one cell per row, but they
#' carry one value per \emph{cluster}, and imputing them row by row invents
#' variation the design cannot have. countimp refuses that and points at the
#' \code{\link{mice.impute.2lonly.pmm}} family, which imputes at the cluster
#' level:
#'
#' \preformatted{
#' pred["classsize", ] <- c(1, 1, 0, -2)      # -2 marks the cluster variable
#' countimp(d, method = c(y = "2l.poisson", classsize = "2lonly.pmm",
#'                        x1 = "", classroom = ""),
#'          predictorMatrix = pred, m = 5)
#' }
#'
#' Available: \code{2lonly.pmm}, \code{2lonly.norm}, \code{2lonly.pois} and
#' \code{2lonly.nb}. A variable that genuinely varies within its clusters is
#' untouched by the check, and a variable observed in fewer than three clusters
#' is skipped -- constancy across two clusters is not evidence of a level. The
#' check can be switched off with
#' \code{options(countimp.check.levels = FALSE)}, which restores the old
#' behaviour; measured on 60 classes with 20\% of the classes missing their
#' class variable, that behaviour destroyed the constancy in 12 of them.
#'
#' Predictors of any level are \emph{used} without special treatment: a
#' cluster-level covariate is simply a column that happens to be constant
#' within its cluster, so \code{y ~ x1 + classsize + schooltype + (1 | school) +
#' (1 | classroom)} needs nothing beyond the formula.
#'
#' @section Zero-inflated and hurdle models:
#' Two-part families model the zeros separately, so they take a second formula.
#' \code{zero} describes that part; the count part stays in \code{formulas}:
#'
#' \preformatted{
#' ## single level: counts from x1, excess zeros from z1
#' countimp(d, formulas = list(y ~ x1),
#'          zero     = ~ z1,
#'          family   = zi_poisson(), m = 5)
#'
#' ## two levels: each part carries its own random structure
#' countimp(d, formulas = list(y ~ x1 + (1 | id)),
#'          zero     = ~ z1 + (1 | id),
#'          family   = zi_nb(), m = 5)
#'
#' ## several incomplete variables: `zero` must be named
#' countimp(d, formulas = list(y1 ~ x1 + (1 | id), y2 ~ x2 + (1 | id)),
#'          zero     = list(y1 = ~ z1 + (1 | id), y2 = ~ z2 + (1 | id)),
#'          family   = list(y1 = zi_nb(), y2 = hurdle_nb()), m = 5)
#' }
#'
#' \strong{Which family models what.} \code{zi_*} treats a zero as coming from
#' either of two sources -- a structural zero or a count that happened to be
#' zero -- and both linear predictors are moved by one draw from the joint
#' posterior. \code{hurdle_*} treats one gate as deciding zero versus positive,
#' and the positive part cannot be zero. Which one is right is a question about
#' the data, not a statistic; \code{\link{countimp_fit_diag}} shows how far
#' apart the two fit, and says out loud that it cannot settle the choice.
#'
#' \strong{The random parts are read separately.} The two parts may differ in
#' their random \emph{intercepts}, and that is what selects the method variant:
#'
#' \tabular{ll}{
#'   \code{y ~ x1 + (1 | id)}, \code{zero = ~ z1 + (1 | id)}
#'     \tab \code{2l.zinb} \cr
#'   \code{y ~ x1 + (0 + x1 | id)}, \code{zero = ~ z1 + (1 | id)}
#'     \tab \code{2l.zinb.noint.count} \cr
#'   \code{y ~ x1 + (1 | id)}, \code{zero = ~ z1 + (0 + z1 | id)}
#'     \tab \code{2l.zinb.noint.zero} \cr
#'   \code{y ~ x1 + (0 + x1 | id)}, \code{zero = ~ z1 + (0 + z1 | id)}
#'     \tab \code{2l.zinb.noint.both}
#' }
#'
#' They may \strong{not} differ in their grouping levels: the model travels on
#' as type codes, which mark a variable as grouping (\code{-2}) without saying
#' for which part. A \code{zero} formula with fewer \code{(1 | g)} terms than
#' the count part is an error rather than a silent completion.
#'
#' \strong{A zero formula for a one-part family} (\code{poisson_count()},
#' \code{nb()}, \code{compois()}, the scale families) is ignored with a
#' warning -- there is no second part to put it in.
#'
#' @return An object of class \code{c("countimp_mids", "mids")}, whatever
#'   packages are loaded. The inherited \code{mids} class keeps the object
#'   usable with \pkg{mice} (see \code{\link[mice]{mice}}) when that package is
#'   available, while the \code{countimp_mids} class ensures countimp's own
#'   methods never displace mice's methods for mice's objects. Use
#'   \code{\link{countimp_diagnostics}} to retrieve what the imputation run
#'   recorded and \code{\link{miinference}} for pooling.
#' @seealso \code{\link{countimp_spec}} to inspect the translation of a
#'   formula specification into the classic representation;
#'   \code{\link{countimp_families}} for the available count distributions;
#'   \code{\link{countimp_diagnostics}} for post-imputation checks.
#' @author Kristian Kleinke
#' @examples
#' data(crim4w)
#'
#' ## --- formula interface (since 2.6.0) -----------------------------------
#' ## one model per incomplete count variable; predictors are read from the
#' ## right-hand side, so no predictorMatrix is needed
#' imp <- countimp(data = crim4w,
#'                 formulas = list(ACRIM = ACRIM ~ FEMALE + RE + GY + HA),
#'                 family   = list(ACRIM = nb()),
#'                 m = 2, maxit = 2, seed = 1)
#'
#' ## a two-part model: counts from one set of predictors, the zeros from
#' ## another. Omitting `zero` would reuse the count part's right-hand side.
#' countimp_spec(data     = crim4w,
#'               formulas = list(ACRIM = ACRIM ~ FEMALE + RE + GY),
#'               zero     = ~ FEMALE + HA,
#'               family   = list(ACRIM = zi_nb()))
#'
#' ## inspect how the specification is translated before imputing
#' countimp_spec(data     = crim4w,
#'               formulas = list(ACRIM = ACRIM ~ FEMALE + RE),
#'               family   = list(ACRIM = hurdle_nb()))
#'
#' ## --- classic interface -------------------------------------------------
#' ini <- countimp(crim4w, maxit = 0)
#' meth <- ini$method
#' meth[6:7] <- "hp"
#' meth[8:9] <- "pmm"
#' pred <- ini$predictorMatrix
#' pred[, "id"] <- 0
#' pred["ACRIM", ] <- c(0, 1, 3, 2, 0, 3, 3, 2, 1)
#' ## "hp" is a hurdle method and needs pscl, which is a suggested package --
#' ## so the example runs it only where pscl is available. Without the guard
#' ## R CMD check fails on a machine that installs hard dependencies only.
#' if (requireNamespace("pscl", quietly = TRUE)) {
#'   imp2 <- countimp(data = crim4w, method = meth, predictorMatrix = pred)
#'   imp2$method
#' }
## The names of the family constructors, in ONE place.
##
## They used to be listed in the error message itself, which named six of
## thirteen -- the scale families from B89 were missing. The same mistake sat
## in `analyse/k27_b12_stand.R`, where it cost more: that script kept its own
## list of seven families and therefore reported twelve method names as "not
## reachable through formulas" when they had been reachable for weeks. The
## project record carried the wrong figure (58 of 81 instead of 70 of 85) as an
## open item for two days.
##
## test-B96 holds this list against the namespace. A hand-maintained list is
## defensible when a test catches up with it; without that test it is a time
## bomb.
.countimp_familien <- function() {
  c("poisson_count", "nb",
    "zi_poisson", "zi_nb", "hurdle_poisson", "hurdle_nb",
    "zerotrunc_poisson", "zerotrunc_nb", "bounded_poisson", "bounded_nb",
    "censored_poisson", "censored_nb", "compois")
}

## What `stats` carries under the same name, and where it belongs here.
##
## `family = poisson()` is the likeliest slip of all: it is the spelling from
## glm(), it exists, and it returns an object of class "family". The old
## message answered that with "must be a family object ... not an object of
## class family" -- a contradiction no reader can resolve. Nor did it mention
## that for quasipoisson there is deliberately NO family.
.countimp_stats_familie <- list(
  poisson      = "poisson_count()",
  quasipoisson = NA_character_,   # deliberately none: see below
  binomial     = NA_character_,
  gaussian     = NA_character_,
  Gamma        = NA_character_,
  inverse.gaussian = NA_character_
)

## Quasi-Poisson deliberately gets no family constructor.
##
## It is not a distribution but a variance assumption on a Poisson core: there
## is no likelihood, hence no posterior for `draw = "bayes"` to draw from. The
## method remains reachable through `method = c(y = "quasipoisson")` -- whoever
## wants it gets it -- but the formula interface leads nobody there, and
## `family = quasipoisson()` now says why.
## Short label for a stats family object, for the message.
##
## Called `.countimp_short_message` on the first attempt -- a name the package already
## uses (diagnostics.R shortens message texts with it). R kept the existing
## definition, and the message came out truncated: "`family$CCRIM` is poisson"
## instead of the full label. With the opposite load order my version would
## have won and broken message shortening throughout the package. Check names
## first, not afterwards.
.countimp_stats_family_short <- function(fam) {
  nm <- if (is.null(fam$family)) "?" else fam$family
  paste0("`", nm, "()` from stats/glm, not a countimp family.")
}

.countimp_familie_falsch <- function(fam) {
  konstruktoren <- paste0(.countimp_familien(), "()", collapse = ", ")

  ## A family object from stats: take its name and translate it.
  if (inherits(fam, "family")) {
    nm <- fam$family
    ziel <- if (!is.null(nm) && nm %in% names(.countimp_stats_familie))
              .countimp_stats_familie[[nm]] else ""
    kopf <- paste0("`family = ", nm, "()` is the family object of stats/glm, ",
                   "not a countimp family.")
    if (identical(nm, "quasipoisson"))
      return(paste0(kopf, "\n  countimp has no quasipoisson family on purpose: ",
                    "quasi-likelihood has no\n  posterior to draw from. The ",
                    "method is still available through the classic\n  ",
                    "interface: `method = c(<variable> = \"quasipoisson\")`."))
    if (!is.na(ziel) && nzchar(ziel))
      return(paste0(kopf, "\n  Use `family = ", ziel, "` instead."))
    return(paste0(kopf, "\n  countimp imputes counts; available families: ",
                  konstruktoren, "."))
  }

  ## A list holding something else: say WHAT and WHERE.
  if (is.list(fam) && length(fam)) {
    schlecht <- which(!vapply(fam, inherits, logical(1), "countimp_family_spec"))
    if (length(schlecht)) {
      wo <- names(fam)[schlecht[1L]]
      wo <- if (is.null(wo) || !nzchar(wo)) paste0("[[", schlecht[1L], "]]")
            else paste0("$", wo)
      element <- fam[[schlecht[1L]]]
      ## A stats family object inside the list gets the same translation as
      ## `family = poisson()` -- the same slip one level down, and the more
      ## common form when several variables are imputed.
      if (inherits(element, "family"))
        return(paste0("`family", wo, "` is ", .countimp_stats_family_short(element),
                      "\n  ", sub("^[^\n]*\n  ", "",
                                  .countimp_familie_falsch(element))))
      was <- paste(class(element), collapse = "/")
      return(paste0("`family", wo, "` is not a countimp family but an object ",
                    "of class ", was, ".\n  Use one of: ", konstruktoren, "."))
    }
  }

  ## A character string: the method name is meant.
  if (is.character(fam) && length(fam) == 1L)
    return(paste0("`family` must be a family object, not a character string.",
                  "\n  For `family = \"", fam, "\"` you probably want ",
                  "`method = c(<variable> = \"", fam, "\")` instead.",
                  "\n  Family constructors: ", konstruktoren, "."))

  paste0("`family` must be a family object (or a named list of them), not an ",
         "object of class ", paste(class(fam), collapse = "/"),
         ".\n  Use the constructors: ", konstruktoren, ".")
}

#' @export
countimp <- function(data, ...)
{
  ## `data` used to travel inside `...`, which meant the formula interface
  ## only recognised it when it was passed by name -- countimp(df, formulas =
  ## ..., family = ...) failed with "`data` is required", while the classic
  ## interface accepted the same positional call. Naming it here makes both
  ## interfaces behave alike; it is put back into the argument list so that
  ## the downstream do.call()s are unchanged.
  arguments <- list(...)
  if (!missing(data)) arguments$data <- data

  ## ---- formula interface (since 2.6.0) ------------------------------------
  ##
  ## A `formulas`/`family` pair that does not match the pattern below used to
  ## fall through this block without a word: `formulas` and `family` were then
  ## dropped further down as unused arguments and every incomplete variable was
  ## imputed by the default method. countimp(d, formulas = list(b ~ a),
  ## family = "poisson") returned pmm imputations and reported success (B81).
  ## The likely spellings are therefore checked before the dispatch.
  if (!is.null(arguments$formulas)) {
    fam <- arguments$family
    ok <- inherits(fam, "countimp_family_spec") ||
      (is.list(fam) && length(fam) &&
         all(vapply(fam, inherits, logical(1), "countimp_family_spec")))
    if (is.null(fam))
      stop("`formulas` needs a matching `family`. Use one of the family ",
           "constructors, e.g. family = poisson_count(), nb(), zi_poisson(), ",
           "hurdle_nb().", call. = FALSE)
    if (!ok) stop(.countimp_familie_falsch(fam), call. = FALSE)
  } else if (inherits(arguments$family, "countimp_family_spec")) {
    stop("`family` was given without `formulas`. The formula interface needs ",
         "both; the classic interface uses `method` and `predictorMatrix`.",
         call. = FALSE)
  }

  ## `formulas` + `family` are compiled into `method`/`predictorMatrix`/`type`
  ## and then handed to the classic engine. Nothing below this block changed.
  if (!is.null(arguments$formulas) && !is.null(arguments$family) &&
      inherits(arguments$family, "countimp_family_spec") ||
      (!is.null(arguments$formulas) && is.list(arguments$family) &&
       length(arguments$family) &&
       inherits(arguments$family[[1L]], "countimp_family_spec"))) {

    if (!is.null(arguments$predictorMatrix) || !is.null(arguments$method))
      stop("Specify the imputation model EITHER through `formulas`/`family` ",
           "OR through `method`/`predictorMatrix`, not both.", call. = FALSE)
    if (is.null(arguments$data))
      stop("`data` is required when using the formula interface.", call. = FALSE)

    spec <- countimp_spec(data = arguments$data,
                          formulas = arguments$formulas,
                          zero = arguments$zero,
                          family = arguments$family,
                          draw = if (is.null(arguments$draw)) "bayes" else arguments$draw)

    ## Scale arguments travelling in the families (bounds, censor) become
    ## ordinary per-variable arguments here -- the form the methods already
    ## read. Giving the same one BOTH ways is an error, not a precedence rule:
    ## the two would differ silently, and an imputation under a bound that was
    ## not meant still looks like plausible counts. Same decision the interface
    ## already takes for `formulas` together with `method`.
    for (nm in names(spec$args)) {
      if (!is.null(arguments[[nm]]))
        stop("`", nm, "` was given twice: once through the family (e.g. ",
             if (nm == "bounds") "bounded_poisson(bounds = ...)"
             else "censored_poisson(censor = ...)",
             ")\n  and once as `", nm, " = ` in the call. Drop one of them -- ",
             "countimp does not\n  pick a winner, because the two can differ ",
             "and the imputations would not show it.", call. = FALSE)
      arguments[[nm]] <- spec$args[[nm]]
    }

    arguments$formulas <- NULL; arguments$zero <- NULL
    arguments$family <- NULL;   arguments$draw <- NULL
    arguments$method <- spec$method
    arguments$predictorMatrix <- spec$predictorMatrix
    ## `type` travels per variable; the engine reads it off the predictorMatrix
    ## row, so encode the codes there (that is what the classic interface does).
    for (y in names(spec$type)) {
      ty <- spec$type[[y]]
      arguments$predictorMatrix[y, names(ty)] <- ty
    }

    ## Grouping variables as integers -- but only on THIS route. The formula
    ## syntax `(1 | g)` is lme4's and glmmTMB's, where a factor is the normal
    ## case: a user who wrote their analysis model there and reuses the formula
    ## for the imputation almost certainly has one. check.data() refuses a
    ## factor as a class variable, and rightly so on the classic route -- there
    ## the -2 is the user's own and a factor could be expanded into indicator
    ## columns further down, so 50 groups would silently become 49 predictors.
    ## Here the variable is marked as a grouping SYNTACTICALLY, before any
    ## predictorMatrix exists, so the conversion is unambiguous and loses
    ## nothing: an identifier carries neither order nor metric.
    pm <- arguments$predictorMatrix
    grp <- colnames(pm)[apply(pm == -2L, 2L, any)]
    for (g in grp) {
      sp <- arguments$data[[g]]
      if (is.factor(sp) || is.character(sp))
        arguments$data[[g]] <- as.integer(as.factor(sp))
    }

    attr(arguments$predictorMatrix, "countimp_spec") <- spec
  }

  ## ---- engine ------------------------------------------------------------
  ## Since 3.0.0 every call runs on countimp's own engine. Up to 2.6.0 the
  ## call was routed to mice::mice unless a zero-inflated count method was
  ## requested, which made countimp's behaviour depend on mice's architecture:
  ## a change to mice's argument list, to its method dispatch or to its return
  ## object changed what countimp did, for models that are entirely countimp's
  ## own. The engine now handles all methods (it always could -- the routing
  ## was historical), and mice's *methods* remain reachable through
  ## .countimp_find_method(), which is a stable, documented interface rather
  ## than an architectural dependency.
  ##
  ## Argument checking is done here rather than in the engine because the
  ## engine takes `...`, which silently absorbs anything it does not know --
  ## including mice-only arguments and plain typos. Before 3.0.0 a call with
  ## `ignore =` was routed to mice and honoured; routing it to the engine
  ## without this check would have silently ignored it and returned different
  ## imputations with no indication that anything was dropped.
  .countimp_check_args(arguments)
  do.call("zimice", arguments)
}
