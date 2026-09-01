## ---------------------------------------------------------------------------
## countimp 2.1: automated imputation diagnostics
##
## Replaces the outlier handling controlled by EV (deprecated in 2.0.8).
##
## Design principle: FLAG AND REPORT, DO NOT SILENTLY REPAIR.
## The old EV mechanism replaced flagged draws by predictive mean matching.
## That truncated the upper tail of the count distribution -- which for
## overdispersed counts is part of the target distribution, not an artefact --
## and reduced coverage of the nominal 95% CI from 94% to 79%.
##
## The checks below separate three genuinely different situations:
##   (1) HARD FAILURE   - the fit or the draw is not a number. Must be repaired
##                        (redraw), because the imputation is unusable.
##   (2) PATHOLOGY      - the fitted model itself predicts impossible means
##                        (linear predictor extrapolated far beyond the data).
##                        Flagged and reported; optionally repaired.
##   (3) HEAVY TAIL     - the fit is fine and the draw is simply large because
##                        the count distribution has a heavy tail. NOT a defect.
##                        Left untouched.
## The deprecated EV mechanism could not distinguish (2) from (3): on clean
## negative binomial data it flagged values in 98% of all draws.
## ---------------------------------------------------------------------------

.countimp_diag <- new.env(parent = emptyenv())
.countimp_diag$log     <- list()
.countimp_diag$enabled <- TRUE

## Number of times a draw is repeated after a HARD failure (non-finite fit,
## non-finite draw). ImputeRobust uses 3 for the same purpose; we keep that.
## Read by .countimp_draw_retry() at the bottom of this file -- until 3.0.0 the
## constant had no reader at all and action = "redraw" was indistinguishable
## from "silent" (B57).
.countimp_max_redraw <- 3L

#' Control and retrieve countimp imputation diagnostics
#'
#' @param enable logical; switch diagnostic recording on or off.
#' @param reset  logical; discard everything recorded so far.
#' @export
countimp_diagnostics <- function(enable = NULL, reset = FALSE) {
  if (!is.null(enable)) .countimp_diag$enabled <- isTRUE(enable)
  if (isTRUE(reset)) .countimp_diag$log <- list()
  ## A pure setter call -- countimp_diagnostics(enable = TRUE) or
  ## (reset = TRUE) -- returns invisibly, so switching the log on does not
  ## print "no diagnostics recorded" at the console.
  setter.only <- (!is.null(enable) || isTRUE(reset)) &&
    length(.countimp_diag$log) == 0L
  if (length(.countimp_diag$log) == 0L) {
    empty <- structure(data.frame(), class = c("countimp_diag", "data.frame"))
    return(if (setter.only) invisible(empty) else empty)
  }
  out <- do.call(rbind, lapply(.countimp_diag$log, as.data.frame,
                               stringsAsFactors = FALSE))
  rownames(out) <- NULL
  structure(out, class = c("countimp_diag", "data.frame"))
}

#' @export
print.countimp_diag <- function(x, ...) {
  if (nrow(x) == 0L) { cat("countimp: no imputation diagnostics recorded.\n"); return(invisible(x)) }
  nb <- sum(x$status != "ok")
  cat(sprintf("countimp imputation diagnostics: %d draws recorded, %d flagged\n",
              nrow(x), nb))
  if (nb > 0L) {
    cat("\nflagged draws:\n")
    print(utils::head(as.data.frame(x)[x$status != "ok", , drop = FALSE], 20L))
    cat("\nInspect these imputations before pooling. A flag does not automatically\n",
        "invalidate the imputation; it marks a draw for visual inspection.\n", sep = "")
  } else {
    cat("all draws passed the automated checks.\n")
  }
  invisible(x)
}

.countimp_record <- function(...) {
  if (!isTRUE(.countimp_diag$enabled)) return(invisible(NULL))
  .countimp_diag$log[[length(.countimp_diag$log) + 1L]] <- list(...)
  invisible(NULL)
}

## Record a numerical-fallback event from an imputation method.
##
## countimp_diagnostics() rbind()s the log entries, so every row has to carry
## the same columns -- a row with a different set of names makes the whole
## call fail, not just that row. This wrapper is the only supported way to
## note an event outside .countimp_screen(): it fills the schema and leaves
## the draw-level numbers at NA, which is what they are for a fallback that
## concerns the fit rather than the draws.
##
## Such events matter to the user: a ridge penalty or a singular Hessian means
## the model was not identified by the data, and imputations from it deserve a
## look even when nothing failed outright.
## Record a numerical event reported by the prediction/draw helpers.
##
## These call sites pass an identifier plus a free-text description. Before
## 3.0.0 they called .countimp_record() with two unnamed arguments, which
## produced log rows with columns "V1"/"V2" instead of the schema. As long as
## only such rows were recorded, countimp_diagnostics() happened to work; as
## soon as one of them met a row from .countimp_screen(), the rbind() failed
## with "numbers of columns of arguments do not match" and the whole
## diagnostics call aborted -- losing every recorded row, including the ones
## the user wanted to see. The event is now folded into the same schema.
.countimp_note_event <- function(event, detail = "") {
  .countimp_note(paste0(event, if (nzchar(detail)) paste0(": ", detail) else ""),
                 where = event)
}

.countimp_note <- function(problem, where = NA_character_) {
  .countimp_record(method     = where,
                   status     = "flagged",
                   problems   = problem,
                   n_imp      = NA_integer_,
                   obs_max    = NA_real_,
                   draw_max   = NA_real_,
                   draw_ratio = NA_real_,
                   mu_ratio   = NA_real_,
                   warning    = "")
}

## --- (1) fit-level checks ---------------------------------------------------
## Returns a character vector of problems; empty if the fit is sound.
.countimp_check_fit <- function(fit, method = NA_character_,
                                theta_hi = 1e6, theta_lo = 1e-6, se_hi = 1e3) {
  p <- character(0)
  if (inherits(fit, "try-error") || is.null(fit)) return("fit_failed")

  ## countimp's own single-level engine passes the list from .countimp_1l_fit()
  ## -- beta/cov/scale/theta, no class, no $coefficients. coef() and vcov()
  ## return NULL on it, which made the checks below report se_nonfinite for a
  ## perfectly sound fit while silently skipping convergence (B59). Translate
  ## that shape into the two things the checks need.
  if (is.null(attr(fit, "class")) && is.list(fit) &&
      all(c("beta", "cov") %in% names(fit))) {
    return(.countimp_check_coefs(fit$beta, fit$cov, fit$scale %||% 1,
                                 fit$theta %||% NA_real_,
                                 theta_hi = theta_hi, theta_lo = theta_lo,
                                 se_hi = se_hi))
  }

  ## Optimiser convergence code. Dispatch on CLASS, not on the presence of a
  ## list element: `fit$fit` is the TMB optimiser list for glmmTMB but the
  ## numeric vector of fitted values for glm.nb, so probing by name silently
  ## picks up the wrong object.
  cv <- NULL
  if (inherits(fit, "glmmTMB")) {
    cv <- tryCatch(fit$fit$convergence, error = function(e) NULL)
  } else if (inherits(fit, c("zeroinfl", "hurdle"))) {
    cv <- tryCatch(as.integer(!isTRUE(fit$converged)), error = function(e) NULL)
  } else if (inherits(fit, "glm")) {
    cv <- tryCatch(as.integer(!isTRUE(fit$converged)), error = function(e) NULL)
  }
  if (length(cv) == 1L && !is.na(cv) && cv != 0) p <- c(p, "not_converged")

  ## fixed effects finite and of plausible magnitude
  be <- try(if (inherits(fit, "glmmTMB")) unlist(glmmTMB::fixef(fit)) else stats::coef(fit),
            silent = TRUE)
  if (inherits(be, "try-error") || !all(is.finite(be))) p <- c(p, "coef_nonfinite")
  else if (max(abs(be)) > 50)                           p <- c(p, "coef_extreme")

  ## standard errors: non-finite SEs indicate a singular or indefinite Hessian.
  ## vcov() returns a plain matrix for glm/glm.nb/zeroinfl but a LIST of
  ## component matrices for glmmTMB -- extract accordingly.
  se <- try(suppressWarnings({
    V <- stats::vcov(fit)
    if (is.list(V) && !is.matrix(V)) V <- V[[1L]]
    sqrt(diag(as.matrix(V)))
  }), silent = TRUE)
  if (inherits(se, "try-error") || !all(is.finite(se))) p <- c(p, "se_nonfinite")
  else if (max(se, na.rm = TRUE) > se_hi)               p <- c(p, "se_extreme")

  ## dispersion parameter degenerate -> the NB has collapsed to Poisson or worse
  ## One accessor for all four conventions (MASS theta, pscl theta, glmmTMB
  ## sigma, stats sigma) -- see .countimp_theta() in foreign-packages.R. Diagnostics
  ## must never stop, so pflicht = FALSE and the try() stays.
  th <- try(.countimp_theta(fit, pflicht = FALSE), silent = TRUE)
  if (!inherits(th, "try-error") && length(th) == 1L && is.finite(th)) {
    if (th > theta_hi) p <- c(p, "theta_degenerate_high")
    if (th < theta_lo) p <- c(p, "theta_degenerate_low")
  }
  p
}

`%||%` <- function(a, b) if (is.null(a)) b else a

## Fit-level checks for countimp's own fit list (beta, cov, scale, theta).
##
## Same defect vocabulary as .countimp_check_fit() so the caller cannot tell
## the two apart -- only the accessors differ. There is no convergence code in
## this shape: glm.fit()'s $converged is dropped by .countimp_1l_fit(), so a
## non-converged fit shows up here as non-finite or extreme coefficients.
.countimp_check_coefs <- function(beta, cov, scale = 1, theta = NA_real_,
                                  theta_hi = 1e6, theta_lo = 1e-6, se_hi = 1e3) {
  p <- character(0)
  if (!length(beta) || !all(is.finite(beta))) p <- c(p, "coef_nonfinite")
  else if (max(abs(beta)) > 50)               p <- c(p, "coef_extreme")

  ## SEs from the unscaled covariance times the dispersion scale -- the same
  ## product .countimp_1l_draw_beta() draws from, so a scale that cannot be
  ## used for the draw is caught here rather than at the Cholesky.
  if (!is.finite(scale) || scale <= 0) p <- c(p, "se_nonfinite")
  else {
    se <- try(suppressWarnings(sqrt(diag(as.matrix(cov) * scale))), silent = TRUE)
    if (inherits(se, "try-error") || !length(se) || !all(is.finite(se)))
      p <- c(p, "se_nonfinite")
    else if (max(se) > se_hi) p <- c(p, "se_extreme")
  }

  if (length(theta) == 1L && is.finite(theta)) {
    if (theta > theta_hi) p <- c(p, "theta_degenerate_high")
    if (theta < theta_lo) p <- c(p, "theta_degenerate_low")
  }
  p
}

## --- (1b) warning capture and classification --------------------------------
## Fitting a count model emits warnings of very different severity. Recording
## all of them is useless: in a multilevel test run 57 of 57 draws carried the
## two glmmTMB family-specification deprecation warnings, which say nothing
## about the imputation. Only warnings matching CRITICAL are recorded; anything
## matching BENIGN is dropped silently; the rest are recorded as "other".
.countimp_warn_critical <- paste(
  "convergence", "Hessian", "singular", "NaN", "NA/NaN", "non-finite",
  "did not converge", "iteration limit", "step size", "infinite",
  "algorithm did not", "rank.deficient", "fitted rates numerically 0",
  "fitted probabilities numerically", "theta", "boundary",
  sep = "|")
.countimp_warn_benign <- paste(
  "components missing from 'family'",     # glmmTMB, triggered by countimp itself
  "plain list is deprecated",             # glmmTMB, ditto
  "Predicting new random effect levels",  # expected: imputing unseen clusters
  "Disabling conditional",
  sep = "|")

#' Evaluate an expression, capturing and classifying warnings
#' @keywords internal
.countimp_quietly <- function(expr) {
  w <- character(0)
  val <- withCallingHandlers(
    tryCatch(expr, error = function(e) structure(conditionMessage(e), class = "countimp_error")),
    warning = function(c) { w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning") })
  w <- w[!grepl(.countimp_warn_benign, w, ignore.case = TRUE)]
  crit <- w[grepl(.countimp_warn_critical, w, ignore.case = TRUE)]
  list(value = val, warnings = w, critical = crit)
}

## --- (2) prediction-level check ---------------------------------------------
## The discriminating statistic. A pathological fit predicts a MEAN far outside
## anything the observed data support; a heavy-tailed but sound fit does not.
## Empirical calibration (see analysis/diagnostik_AB.csv):
##   sound flat-file NB     max(mu_mis)/max(y_obs) <=   1.3
##   sound multilevel NB    max(mu_mis)/max(y_obs) <=  21.7
##   genuine extrapolation  max(mu_mis)/max(y_obs) >= 4600
## The default threshold of 100 sits inside a gap spanning two orders of
## magnitude on either side.
.countimp_check_mu <- function(mu, y_obs, mu_ratio_max = 100) {
  p <- character(0)
  if (!length(mu) || all(is.na(mu))) return("mu_unavailable")
  if (any(!is.finite(mu)))           p <- c(p, "mu_nonfinite")
  om <- suppressWarnings(max(y_obs, na.rm = TRUE))
  if (!is.finite(om) || om <= 0) om <- 1
  r <- suppressWarnings(max(mu, na.rm = TRUE)) / om
  if (is.finite(r) && r > mu_ratio_max) p <- c(p, "mu_extrapolated")
  attr(p, "mu_ratio") <- r
  p
}

## --- (3) draw-level check ---------------------------------------------------
## Only non-finite draws are HARD defects. Large finite draws are never
## discarded -- they are the tail of the target distribution, and dropping the
## large ones until a small one appears would pull the imputations downwards.
##
## But a draw can be finite and still absurd, and that is worth saying. The
## ratio of the largest imputed value to the largest observed one separates the
## two cases cleanly. Measured across families and data situations, each with a
## model that FITS:
##
##   Poisson data, poisson       1.3        NB data (theta = 2), nb      1.3
##   Poisson data, nb            2.6        NB data (theta = 0.5), nb    1.0
##   zero-inflated data, zip     1.7        zero-inflated data, nb       3.7
##
## Against that, the case this check was written for: the delinquency data in
## crim4w, imputed under nb with the other waves as untransformed predictors --
## ratio 59.7, largest imputed value 955 against 16 observed. The cause is not
## the family but the link: counts as predictors in a log-link model extrapolate
## exponentially (eta reaches 7.2 at the observed maxima, so mu = 1380), and
## drawing beta widens that further.
##
## A threshold of 10 sits an order of magnitude above every sound case and an
## order of magnitude below the broken one. The draw is kept and flagged, which
## is what "flagged" means in this package: usable but unusual, recorded,
## warned about, and left alone.
.countimp_check_draws <- function(imp, y_obs, draw_ratio_max = 10) {
  p <- character(0)
  if (!length(imp)) return("no_draws")
  if (any(is.na(imp) | !is.finite(imp))) p <- c(p, "draw_nonfinite")
  om <- suppressWarnings(max(y_obs, na.rm = TRUE)); if (!is.finite(om) || om <= 0) om <- 1
  r <- suppressWarnings(max(imp, na.rm = TRUE)) / om
  if (is.finite(r) && r > draw_ratio_max) p <- c(p, "draw_far_above_observed")
  attr(p, "draw_ratio") <- r
  p
}

## --- orchestrator -----------------------------------------------------------
#' Run the automated checks for one imputation draw
#'
#' @details
#' Every draw taken by countimp's count and zero-inflated engines passes through
#' this function; you rarely need to call it yourself. It sorts problems into two
#' severities, and that distinction governs what happens next:
#'
#' \describe{
#'   \item{\code{"hard"}}{the draw is unusable -- the model did not fit, or the
#'     coefficients, the predicted mean or the draws themselves are not finite.
#'     countimp repeats the draw up to three times. Only if every attempt fails
#'     does the imputation stop, with a message naming the reason and the number
#'     of attempts. Before version 3.0.0 the first failure aborted the whole run,
#'     discarding every completed draw with it.}
#'   \item{\code{"flagged"}}{the draw is usable but unusual -- large but finite
#'     values, extreme standard errors, a predicted mean far outside the range of
#'     the observed data, a degenerate dispersion parameter. Such draws are
#'     recorded, warned about, and \strong{kept}. They are not repeated:
#'     discarding large draws until a small one appears would pull the imputations
#'     towards the centre of the distribution and understate the between-
#'     imputation variance, which is the opposite of what multiple imputation is
#'     for.}
#' }
#'
#' Use \code{\link{countimp_diagnostics}} to read the log afterwards -- it holds
#' one row per draw, including the attempts that were discarded.
#'
#' @param imp     the imputed values.
#' @param y_obs   the observed part of the target variable.
#' @param fit     the fitted imputation model (optional).
#' @param mu      predicted means for the missing cases (optional).
#' @param method  name of the imputation method, for the diagnostic log.
#' @param action  what to do about the problems found. \code{"warn"} (default)
#'   warns about all of them; \code{"silent"} records them without warning;
#'   \code{"redraw"} is the mode the imputation engines use -- it warns about
#'   flagged draws but stays quiet about hard failures, because those are about
#'   to be repeated by the retry loop and warning about a draw the user never
#'   receives would be noise. All three modes record to the diagnostic log and
#'   return the same \code{status} and \code{redraw} verdict; they differ only
#'   in what they report.
#' @param mu_ratio_max threshold for the extrapolation check.
#' @param warnings character vector of warning messages caught while fitting
#'   the imputation model. A non-empty vector adds \code{"critical_warning"}
#'   to the reported problems, so that a model that converged only nominally
#'   is not treated as clean. Default \code{character(0)}.
#' @return a list with \code{imp}, \code{status}, \code{problems}, \code{redraw}.
#' @export
countimp_check <- function(imp, y_obs, fit = NULL, mu = NULL,
                           method = NA_character_, action = c("warn", "silent", "redraw"),
                           mu_ratio_max = 100, warnings = character(0)) {
  action <- match.arg(action)
  pf <- if (!is.null(fit)) .countimp_check_fit(fit, method) else character(0)
  pm <- if (!is.null(mu))  .countimp_check_mu(mu, y_obs, mu_ratio_max) else character(0)
  pd <- .countimp_check_draws(imp, y_obs)
  pw <- if (length(warnings)) "critical_warning" else character(0)

  problems <- c(pf, as.character(pm), as.character(pd), pw)
  hard     <- any(c("fit_failed", "draw_nonfinite", "mu_nonfinite", "coef_nonfinite") %in% problems)
  status   <- if (length(problems) == 0L) "ok" else if (hard) "hard" else "flagged"

  .countimp_record(
    method     = method,
    status     = status,
    problems   = if (length(problems)) paste(problems, collapse = ",") else "",
    n_imp      = length(imp),
    obs_max    = suppressWarnings(max(y_obs, na.rm = TRUE)),
    draw_max   = suppressWarnings(max(imp,  na.rm = TRUE)),
    draw_ratio = attr(pd, "draw_ratio") %||% NA_real_,
    mu_ratio   = attr(pm, "mu_ratio")   %||% NA_real_,
    warning    = if (length(warnings)) substr(paste(warnings, collapse = " | "), 1, 200) else "")

  ## The three action values must differ, or the documented contract is a lie
  ## (B57). "warn" reports everything; "silent" reports nothing; "redraw" leaves
  ## HARD failures to the caller's retry loop -- warning about a draw that is
  ## about to be repeated would report a problem the user never receives -- but
  ## still warns about flagged draws, which the retry does not repair.
  melden <- switch(action,
    warn   = length(problems) > 0L,
    silent = FALSE,
    redraw = length(problems) > 0L && !hard)
  if (melden) {
    warning(sprintf("countimp: imputation draw flagged (%s): %s. See countimp_diagnostics().",
                    method, paste(problems, collapse = ", ")), call. = FALSE)
  }
  list(imp = imp, status = status, problems = problems, redraw = hard)
}


## Strip countimp's own prefix from a nested message, so the abort message
## reads "... failed on 3 successive draws (the covariance matrix is not
## positive definite)" rather than repeating "countimp:" inside its own text.
.countimp_short_message <- function(msg, n = 160L) {
  s <- trimws(sub("^countimp: *", "", trimws(as.character(msg)[1L])))
  s <- sub("\n.*$", "", s)
  if (nchar(s) > n) paste0(substr(s, 1L, n - 1L), "\u2026") else s
}

## --- retry loop -------------------------------------------------------------
## The reader of .countimp_max_redraw, and the single place where a HARD draw
## failure is retried instead of aborting the whole imputation run (B56).
##
## Why here and not in the engines: mice.impute.*() must return a plain vector
## of draws, and the B40 tests pin the signature and return value of
## .countimp_1l_count(). Wrapping the engine keeps both intact -- the engine
## still does one draw, this function decides whether that draw stands.
##
## What counts as HARD is countimp_check()'s judgment, not ours: fit_failed,
## coef_nonfinite, mu_nonfinite, draw_nonfinite. A flagged draw (large but
## finite, extreme SEs, degenerate theta) is NOT retried -- retrying it would
## quietly bias the imputations towards the centre of the distribution, which
## is the opposite of what multiple imputation is for. Those draws are recorded
## and warned about, and they stand.
##
## zieh:  function() -> list(imp = <draws>, fit = <fit or NULL>, mu = <mu or NULL>)
##        It must draw afresh on every call, using the ambient RNG state, so
##        that a repeat is a genuinely new draw rather than the same one again.
.countimp_draw_retry <- function(zieh, y_obs, method = NA_character_,
                                 mu_ratio_max = 100,
                                 max_redraw = .countimp_max_redraw) {
  letzte <- NULL
  for (versuch in seq_len(max(1L, max_redraw))) {
    q <- .countimp_quietly(zieh())
    z <- q$value

    ## A hard error inside the engine is a hard failure like any other: record
    ## it through the normal path and retry, rather than aborting the run. Both
    ## branches fall through to the shared logging below -- an earlier version
    ## used next here and silently skipped it for exactly the case that needs
    ## the record most.
    if (inherits(z, "countimp_error")) {
      ## fit = NULL means "no fit information", which countimp_check() skips --
      ## the engine's failure would then be classified from the empty draw alone
      ## as no_draws, i.e. FLAGGED, and warned about even under action =
      ## "redraw". A try-error is the vocabulary countimp_check() already has
      ## for this, and it yields fit_failed -> HARD, which is what happened.
      letzte <- countimp_check(imp = numeric(0), y_obs = y_obs,
                               fit = structure(as.character(z), class = "try-error"),
                               mu = NULL, method = method, action = "redraw",
                               mu_ratio_max = mu_ratio_max,
                               warnings = q$critical)
      letzte$fehler <- .countimp_short_message(as.character(z))
    } else {
      letzte <- countimp_check(imp = z$imp, y_obs = y_obs, fit = z$fit,
                               mu = z$mu, method = method, action = "redraw",
                               mu_ratio_max = mu_ratio_max,
                               warnings = q$critical)
      letzte$fehler <- NULL
      if (!isTRUE(letzte$redraw)) {
        ## .countimp_quietly() muffles every warning so the engines can be
        ## retried without shouting about attempts nobody receives. But warnings
        ## of the ACCEPTED draw are the user's business -- countimp raises some
        ## deliberately (the underdispersion fallback in .countimp_rqpois(),
        ## B04, warns once per session). Wrapping the engines without this
        ## silenced them, so re-emit them here. Warnings of discarded attempts
        ## stay muffled: that draw is gone, and reporting it would describe a
        ## result the user never sees.
        for (wm in q$warnings) warning(wm, call. = FALSE)
        return(z$imp)
      }
    }

    if (versuch < max(1L, max_redraw))
      .countimp_note_event("draw_repeated",
        sprintf("%s: attempt %d of %d failed (%s)", method, versuch,
                max(1L, max_redraw),
                if (!is.null(letzte$fehler)) letzte$fehler
                else paste(letzte$problems, collapse = ",")))
  }

  ## Every attempt failed hard. Now -- and only now -- abort, with the reason
  ## and the fact that it was not a one-off.
  urs <- if (!is.null(letzte$fehler)) letzte$fehler else
    paste(letzte$problems, collapse = ", ")

  ## Underdispersion is the one cause the generic advice actively misleads
  ## about. A negative binomial contains the Poisson only as the limit
  ## theta -> Inf, so for variance < mean the likelihood keeps rising in theta
  ## and there is no finite estimate: the fit fails, and it fails on data whose
  ## predictors may be perfectly well behaved. Reported on 28 August 2026 from
  ## study 5, where `bnb` aborted on independent standard normal predictors
  ## with a message naming collinearity -- and a user with sum scores over K
  ## items has underdispersed data ALWAYS, so they would go and delete
  ## predictors.
  ##
  ## Said only when it applies, and only for the NB methods: it is an addition
  ## to the diagnosis, not a replacement, because underdispersion and a
  ## genuinely collinear design can occur together.
  vm <- NA_real_
  yo <- y_obs[is.finite(y_obs)]
  if (length(yo) >= 5L && mean(yo) > 0) vm <- stats::var(yo) / mean(yo)
  unter <- identical(method, "negbin") || identical(method, "cnb")
  unter <- unter && is.finite(vm) && vm < 1
  ## The inner message carries its own guess at a cause ("Collinear predictors
  ## are the usual cause; drop one of them"). Left standing it would appear
  ## immediately before the sentence saying that dropping predictors will not
  ## help -- so where the cause IS known, keep the technical finding and drop
  ## the guess with it.
  if (unter) urs <- sub("\\. .*$", "", urs)
  zusatz <- if (unter) sprintf(paste0(
    "\n  The observed values are UNDERDISPERSED (variance/mean = %.2f). A ",
    "negative binomial\n  has no solution there -- it contains the Poisson ",
    "only as the limit theta -> Inf, so\n  for variance < mean the ",
    "likelihood keeps rising in theta and no finite estimate\n  exists. This ",
    "is a property of the data, not of the predictors: dropping or\n  ",
    "rescaling them will not help. For underdispersed counts use \"cmp\" ",
    "(Conway-\n  Maxwell-Poisson); for bounded counts with a known upper ",
    "limit, \"bp\"."), vm) else ""

  stop(sprintf(paste0(
    "countimp: the imputation model for method \"%s\" failed on %d successive ",
    "draws (%s). The imputation was not completed.%s%s ",
    "See countimp_diagnostics() for what each attempt reported."),
    method, max(1L, max_redraw), urs, zusatz,
    if (nzchar(zusatz)) "" else paste0(
      " Extreme predictor values, a separated model or too few observed cases",
      " are the usual causes; rescale or drop collinear predictors, or choose",
      " a simpler count model.")), call. = FALSE)
}
