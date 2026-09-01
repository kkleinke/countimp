## predict2l.R -- predictions for units in clusters the model never saw
##
## The problem. In a two-level imputation the incomplete variable may be
## missing for every unit of some cluster j. That cluster then contributes
## nothing to the fit, and no random effect u_j is estimated for it -- yet
## those are exactly the units we must impute.
##
## The two backends countimp uses disagree on what to do:
##
##   glmmTMB  predict(..., allow.new.levels = TRUE) silently uses u_j = 0,
##            i.e. it imputes at the population average.
##   glmmPQL  predict(..., level = 1) returns NA for such units. The NA then
##            travels into .pmm.match(), where sort.int(d, partial = donors)
##            fails with "'partial' index 5 outside bounds" -- a message that
##            points at the donor pool and hides the real cause.
##
## Neither is right for multiple imputation. Setting u_j = 0 is a conditional
## mean: it makes imputed units of unseen clusters more alike than observed
## units are, shrinks between-cluster variance, and biases the intraclass
## correlation of the completed data downwards. Proper imputation draws the
## unknown effect from its own distribution,
##
##     u_j* ~ N(0, tau^2)      tau^2 = estimated between-cluster variance
##
## which is the predictive distribution given the model -- the same principle
## the rest of the package follows when drawing beta* from its posterior.
##
## .countimp_predict_2l() implements this for both backends and reports how
## many units were affected, so the user is not silently handed
## population-average imputations.

## Between-cluster SD of the random intercept.
##
## `grp` names the grouping factor to read; with several levels in the model
## (school and class, say) they have different variances and taking the first
## for both would draw the inner effect from the outer distribution. Measured
## on a three-level fit: 0.419 for schools against 0.201 for classes, i.e. the
## wrong one is off by a factor of two. NULL keeps the old behaviour of
## returning the first, which is what every single-level caller wants.
.countimp_tau <- function(fit, grp = NULL) {
  out <- NA_real_
  if (inherits(fit, "glmmTMB")) {
    vc <- try(glmmTMB::VarCorr(fit)$cond, silent = TRUE)
    if (!inherits(vc, "try-error") && length(vc)) {
      i <- if (is.null(grp)) 1L else match(as.character(grp)[1L], names(vc))
      ## An unmatched name is a programming error in countimp, not a user
      ## situation: the caller took the name from the fit's own formula.
      if (is.na(i))
        stop("countimp: the fitted model has no random-effect term for '",
             grp, "'. Terms: ", paste(names(vc), collapse = ", "), ".",
             call. = FALSE)
      out <- sqrt(as.numeric(vc[[i]][1L, 1L]))
    }
  } else if (inherits(fit, "lme")) {
    vc <- try(nlme::VarCorr(fit), silent = TRUE)
    if (!inherits(vc, "try-error"))
      out <- suppressWarnings(as.numeric(vc[1, "StdDev"]))
  }
  if (!is.finite(out)) NA_real_ else out
}

## Observed levels per grouping factor, in the shape .countimp_predict_2l()
## wants. dat[[grp]] with a VECTOR grp does not do this: [[ ]] on a character
## vector indexes recursively (dat$school$class) and fails. One helper so the
## five call sites cannot each get it wrong in their own way.
.countimp_obs_levels <- function(dat, grp) {
  stats::setNames(lapply(grp, function(g) unique(as.character(dat[[g]]))), grp)
}


## `grp` may name SEVERAL grouping levels since the three-level extension, and
## `obs_levels` is then a list with one entry per level. With more than one
## level a unit can be known on one and new on another -- a new class inside a
## known school -- and those are different situations: the school effect is
## estimated and must be kept, only the class effect has to be drawn.
##
## The glmmTMB path therefore works on the LINK scale and adds one draw per new
## level, instead of replacing the whole prediction by the population average.
## For a single level the two are identical, because predict(allow.new.levels =
## TRUE) sets an unseen effect to zero and the population average is exactly
## the prediction with that effect at zero -- so this generalises the old code
## rather than changing it. It also saves the second predict() call.
.countimp_predict_2l <- function(fit, newdata, type = "response",
                                 grp = NULL, obs_levels = NULL,
                                 draw_new = TRUE) {
  ## One level given the old way (a bare vector of levels) is wrapped, so both
  ## calling conventions stay valid.
  if (!is.null(grp) && !is.null(obs_levels) && !is.list(obs_levels)) {
    if (length(grp) > 1L)
      stop("countimp: `obs_levels` must be a list with one entry per grouping ",
           "level when `grp` names ", length(grp), " of them.", call. = FALSE)
    obs_levels <- stats::setNames(list(obs_levels), grp)
  }

  if (!inherits(fit, "glmmTMB")) {
    ## nlme/glmmPQL: one grouping level, and predict(level = 1) returns NA for
    ## unseen ones rather than a population average. Kept as it was.
    if (!is.null(grp) && length(grp) > 1L)
      stop("countimp: the nlme backend handles one grouping level, got ",
           length(grp), ".", call. = FALSE)
    p <- suppressWarnings(stats::predict(fit, newdata = newdata, type = type,
                                         level = 1, na.action = stats::na.pass))
    isnew <- if (!is.null(grp) && !is.null(obs_levels))
      !(as.character(newdata[[grp]]) %in% as.character(obs_levels[[grp]]))
      else !is.finite(p)
    if (!any(isnew)) return(p)
    p0 <- suppressWarnings(stats::predict(fit,
            newdata = newdata[isnew, , drop = FALSE], level = 0,
            na.action = stats::na.pass))
    tau <- .countimp_tau(fit)
    if (draw_new && is.finite(tau) && tau > 0 && !is.null(grp)) {
      lev  <- as.character(newdata[[grp]])[isnew]
      ulev <- unique(lev)
      ustar <- stats::setNames(stats::rnorm(length(ulev), 0, tau), ulev)
      p0 <- p0 + ustar[lev]
    }
    p[isnew] <- p0
    .countimp_note_event("new_cluster_levels",
      sprintf("%d unit(s) in %d cluster(s) absent from the fit; random effect drawn from N(0, %.3f^2)",
              sum(isnew), length(unique(as.character(newdata[[grp]])[isnew])),
              if (is.finite(tau)) tau else 0))
    return(p)
  }

  eta <- stats::predict(fit, newdata = newdata, type = "link",
                        allow.new.levels = TRUE)
  linkinv <- tryCatch(stats::family(fit)$linkinv, error = function(e) identity)
  fertig <- function(e)
    if (identical(type, "response")) linkinv(e) else e

  if (is.null(grp) || is.null(obs_levels)) return(fertig(eta))

  for (g in grp) {
    lev <- as.character(newdata[[g]])
    isnew <- !(lev %in% as.character(obs_levels[[g]]))
    if (!any(isnew)) next
    tau <- .countimp_tau(fit, g)
    if (draw_new && is.finite(tau) && tau > 0) {
      ulev  <- unique(lev[isnew])
      ustar <- stats::setNames(stats::rnorm(length(ulev), 0, tau), ulev)
      eta[isnew] <- eta[isnew] + ustar[lev[isnew]]
    }
    .countimp_note_event("new_cluster_levels",
      sprintf("%d unit(s) in %d cluster(s) of '%s' absent from the fit; random effect drawn from N(0, %.3f^2)",
              sum(isnew), length(unique(lev[isnew])), g,
              if (is.finite(tau)) tau else 0))
  }
  fertig(eta)
}


## Prediction at a DRAWN coefficient vector (glmmTMB) ------------------------
##
## Proper multiple imputation predicts at a draw b* ~ N(beta.hat, V), not at
## the point estimate: that is where between-imputation parameter variance
## comes from. glmmTMB has no "predict at these coefficients" argument of the
## kind lm() users expect, so countimp up to 2.6.0 overwrote an internal slot,
##
##     fitmis$fit$par[names(fitmis$fit$par) == "beta"] <- b.star
##
## and predicted from the mutated object. Verified to work with glmmTMB 1.1.14
## -- but it rests on an undocumented internal, and it leaves the object
## self-contradictory: predict.glmmTMB() reads fit$fit$par, while fixef(),
## ranef() and VarCorr() read fit$fit$parfull, which still held beta.hat.
##
## Two consequences, both silent:
##   1. any code path that reads fixef() after the injection sees the point
##      estimate, not the draw;
##   2. if a future glmmTMB has predict() read parfull instead, imputations
##      revert to the point estimate. Every imputation would then be drawn at
##      the same coefficients, between-imputation variance would collapse, and
##      pooled standard errors would be too small -- with nothing failing.
##
## So we write both slots (verified prediction-neutral relative to writing
## par alone, and it leaves fixef() consistent with the draw), and we check
## once per session that the injection actually reaches predict(). If a future
## glmmTMB ignores it, countimp stops with a named cause instead of quietly
## understating uncertainty.

.countimp_inject_beta <- function(fit, b.star) {
  if (!inherits(fit, "glmmTMB"))
    stop("countimp: .countimp_inject_beta() expects a glmmTMB fit, got ",
         paste(class(fit), collapse = "/"), ".", call. = FALSE)
  b.star <- as.numeric(b.star)
  i <- which(names(fit$fit$par) == "beta")
  if (!length(i))
    stop("countimp: this glmmTMB fit has no 'beta' entry in fit$fit$par, so the\n",
         "  drawn coefficients cannot be passed to predict(). This indicates a\n",
         "  change in glmmTMB's internal layout -- please report it.",
         call. = FALSE)
  if (length(i) != length(b.star))
    stop(sprintf(paste0("countimp: drew %d coefficient(s) for %d 'beta' slot(s).\n",
                        "  This is a bug in countimp -- please report it."),
                 length(b.star), length(i)), call. = FALSE)
  out <- fit
  out$fit$par[i] <- b.star
  j <- which(names(out$fit$parfull) == "beta")   # keep fixef()/ranef() consistent
  if (length(j) == length(b.star)) out$fit$parfull[j] <- b.star
  .countimp_check_beta_route(fit)
  out
}

## One-time check: does overwriting beta change what predict() returns?
## `predictor` exists so that the check itself can be tested: a guard that
## cannot be shown to fire is not a guard. Tests pass a function that ignores
## the mutated object -- the exact failure mode -- and assert that we stop.
.countimp_check_beta_route <- function(fit, predictor = NULL) {
  if (isTRUE(.countimp_state$beta_route_ok)) return(invisible(TRUE))
  if (is.null(predictor))
    predictor <- function(object, newdata)
      stats::predict(object, newdata = newdata, type = "link",
                     na.action = stats::na.pass, allow.new.levels = TRUE)
  ok <- tryCatch({
    nd <- stats::model.frame(fit)
    nd <- nd[seq_len(min(3L, nrow(nd))), , drop = FALSE]
    p0 <- predictor(fit, nd)
    shifted <- fit
    i <- which(names(shifted$fit$par) == "beta")
    shifted$fit$par[i] <- shifted$fit$par[i] + 1      # a full unit on the link scale
    j <- which(names(shifted$fit$parfull) == "beta")
    if (length(j) == length(i)) shifted$fit$parfull[j] <- shifted$fit$parfull[j] + 1
    p1 <- predictor(shifted, nd)
    any(is.finite(p0)) && !isTRUE(all.equal(as.numeric(p0), as.numeric(p1)))
  }, error = function(e) NA)
  if (isTRUE(ok)) {
    .countimp_state$beta_route_ok <- TRUE
    return(invisible(TRUE))
  }
  if (is.na(ok)) return(invisible(NA))   # could not probe; do not block the run
  stop("countimp: this version of glmmTMB (", as.character(utils::packageVersion("glmmTMB")),
       ") ignores the drawn\n",
       "  regression coefficients, so every imputation would be drawn at the same\n",
       "  point estimate. Between-imputation variance would be lost and pooled\n",
       "  standard errors too small, which is why countimp stops here rather than\n",
       "  returning imputations that look fine.\n",
       "  Please report this (the internal layout of glmmTMB fits has changed);\n",
       "  installing the previous glmmTMB version restores correct behaviour.",
       call. = FALSE)
}


## ---------------------------------------------------------------------------
## Joint draw of (beta, u, log theta) for two-level fits
##
## The problem this solves. Proper imputation of y_ij in a two-level count
## model needs a draw from the posterior predictive distribution
##
##     y_ij* ~ f(y | eta_ij*),   eta_ij* = x_ij' beta* + z_ij' u_j*
##
## with (beta*, u_j*) drawn from their JOINT posterior. countimp up to 2.6.0
## drew beta* from N(beta.hat, V.beta) and then called predict() on the fit
## carrying beta*. That is not the same thing, and the difference is large:
## predict.glmmTMB() re-derives the conditional modes of u given the parameter
## vector it is handed, so shifting beta moves eta by only
##
##     sigma^2 / (n_j tau^2 + sigma^2)
##
## of the intended amount -- the standard shrinkage factor. Measured on
## n = 300, 15 clusters of 20: an intercept shift of 1.0 moved eta by 0.0909,
## and 0.0909 is exactly that factor. The random effects absorb the rest.
##
## Consequence for the user: sd(eta) between imputations came out 4-8 times
## too small across configurations (n = 300-600, 10-30 units per cluster).
## Between-imputation variance was therefore far too small, pooled standard
## errors too narrow, and coverage of intervals from mice::pool() below
## nominal -- with nothing failing and no warning.
##
## What is drawn here. glmmTMB fits by Laplace approximation, so the joint
## posterior of the parameter vector AND the random effects is available as a
## Gaussian with precision matrix P = sdreport(getJointPrecision = TRUE):
##
##     (beta, u, log theta, ...) ~ N(point estimates, P^{-1})
##
## Drawing z ~ N(0, P^{-1}) and forming eta* = eta.hat + X z_beta + Z z_u
## reproduces the target exactly. Verified: sd(eta*) matches the documented
## se.fit from predict(..., se.fit = TRUE) to 2e-10 relative, for random
## intercepts, one random slope and two random slopes; the simulated draw
## reaches 97.7-101.4% of that target. Note that se.fit is the right
## reference precisely because it includes the conditional variance of u --
## adding Var(u_j) from ranef(condVar = TRUE) to x'Vx instead overstates the
## target by a factor of 2.2, because beta and u are posterior-correlated.
##
## log theta is drawn from the same matrix when the family has it, which makes
## the dispersion draw coherent with the mean draw. The posterior correlation
## between beta and log theta is small (|r| < 0.07 in the configurations
## measured), so this is not where the gain lies -- it costs nothing and
## removes the need for a separate .countimp_draw_theta() call at 2l sites.
##
## We draw by solving U z = w with P = U'U rather than inverting P: it is
## cheaper, and chol() failing is exactly the signal that P is not positive
## definite, which is the condition under which we must not draw at all.
## Cost of the whole route is comparable to one predict() call (0.014 s at
## n = 300, 0.11 s at n = 3000).

## Random-effect design matrix for newdata, in the column order of getME(fit, "b").
## glmmTMB orders b cluster-major: for each grouping factor, level 1's terms
## come first, then level 2's, and so on -- i.e. rowwise over the ranef matrix,
## not columnwise. Verified against getME(fit, "Z") (exact equality) for
## (1|g), (1+x|g) and (1+x+w|g). Rows of unobserved clusters are left at zero:
## their effect is supplied separately by .countimp_predict_2l().
.countimp_re_design <- function(fit, newdata, component = c("cond", "zi")) {
  component <- match.arg(component)
  rt <- fit$modelInfo$reTrms[[component]]
  if (is.null(rt) || !length(rt$cnms)) return(NULL)
  blocks <- vector("list", length(rt$cnms))
  for (k in seq_along(rt$cnms)) {
    gf  <- names(rt$cnms)[k]
    cn  <- rt$cnms[[k]]
    lev <- levels(rt$flist[[gf]])
    if (is.null(lev) || !length(lev)) return(NULL)
    if (is.null(newdata[[gf]])) return(NULL)          # grouping factor absent
    M <- vapply(cn, function(nm) {
      if (identical(nm, "(Intercept)")) rep(1, nrow(newdata))
      else if (is.null(newdata[[nm]])) rep(NA_real_, nrow(newdata))
      else as.numeric(newdata[[nm]])
    }, numeric(nrow(newdata)))
    M <- matrix(M, nrow = nrow(newdata))
    if (anyNA(M)) return(NULL)                        # a slope variable is missing
    j  <- match(as.character(newdata[[gf]]), lev)
    Zk <- matrix(0, nrow(newdata), length(lev) * length(cn))
    for (i in which(!is.na(j)))
      Zk[i, (j[i] - 1L) * length(cn) + seq_along(cn)] <- M[i, ]
    blocks[[k]] <- Zk
  }
  do.call(cbind, blocks)
}

## Draw eta and theta jointly. Returns a list(eta, mu, theta, joint) where
## `joint` reports whether the joint route was used, so callers and tests can
## tell the two regimes apart. On any obstacle we fall back to the old
## beta-only route: an imputation with understated between-imputation variance
## is still usable and is what previous versions produced, whereas stopping
## would make the package unusable the moment a backend changes.
## theta.min matches the default of .countimp_draw_theta(): the mean of the
## zero-truncated count distribution grows without bound as theta -> 0, which
## is the mechanism behind the imputation explosions recorded as B09a. The
## floor is a documented, deliberate bias.
.countimp_draw_2l <- function(fit, newdata, grp = NULL, obs_levels = NULL,
                              theta = NULL, theta.min = 0.25) {
  linkinv <- tryCatch(stats::family(fit)$linkinv, error = function(e) exp)
  if (!is.function(linkinv)) linkinv <- exp

  fallback <- function(why) {
    if (!isTRUE(.countimp_state$joint_noted)) {
      .countimp_state$joint_noted <- TRUE
      .countimp_note_event("joint_draw_unavailable", why)
    }
    beta <- glmmTMB::fixef(fit)$cond
    ## The last resort has its own way of failing, and it used to fail
    ## unreadably: chol() reports "the leading minor of order 1 is not
    ## positive", which names the arithmetic and not the situation. Measured on
    ## a three-level call (two columns coded -2) whose data carry no such
    ## structure -- the model is richer than the data support, the covariance
    ## of the fixed effects comes out singular, and that message is what the
    ## user saw.
    ##
    ## Same outcome as before -- the run stops, because drawing beta at its
    ## point estimate would understate the parameter uncertainty (B01) and
    ## silently so. What changes is that the message says what happened. The
    ## core words it the same way (glmm.hpp, glmm_slopes.hpp): the covariance
    ## is singular, and a variance component at zero is the usual cause.
    rv <- try(t(chol(stats::vcov(fit)$cond)), silent = TRUE)
    if (inherits(rv, "try-error")) {
      gf <- tryCatch(names(fit$modelInfo$reTrms$cond$cnms),
                     error = function(e) character(0))
      stop("countimp: the covariance of the fixed effects is singular, so no ",
           "draw can be taken.\n  The model is richer than these data support ",
           "-- a variance component at zero is the\n  usual cause",
           if (length(gf) > 1L)
             paste0(", and with ", length(gf), " grouping levels (",
                    paste(gf, collapse = ", "), ") it is the\n  likeliest one: ",
                    "code only one of them as -2, or use a method for a single\n",
                    "  level")
           else "",
           ".\n  countimp stops rather than drawing beta at its point ",
           "estimate, which would\n  understate the parameter uncertainty ",
           "without showing it (B01).", call. = FALSE)
    }
    b.star <- as.numeric(beta) + as.numeric(rv %*% stats::rnorm(ncol(rv)))
    eta <- .countimp_predict_2l(.countimp_inject_beta(fit, b.star), newdata,
                                type = "link", grp = grp, obs_levels = obs_levels)
    list(eta = as.numeric(eta), mu = as.numeric(linkinv(eta)),
         theta = if (is.null(theta)) .countimp_draw_theta(fit) else theta,
         joint = FALSE)
  }

  ## Both routes below are glmmTMB-specific (the fallback goes through
  ## .countimp_inject_beta). A different backend here is a programming error in
  ## countimp, not a user situation, so name it instead of failing later inside
  ## fixef() with a message that points nowhere.
  if (!inherits(fit, "glmmTMB"))
    stop("countimp: .countimp_draw_2l() expects a glmmTMB fit, got ",
         paste(class(fit), collapse = "/"), ".", call. = FALSE)

  ## The grouping factors must be present in newdata. Without them neither
  ## route works -- predict.glmmTMB() cannot even build the model frame -- so
  ## this is not a fallback situation but a programming error in countimp
  ## (every call site sets colnames(newdatamis) <- nam, which includes the
  ## grouping column). Name it here; otherwise it surfaces as
  ## "object 'g' not found" from deep inside model.frame().
  gfs <- names(fit$modelInfo$reTrms$cond$cnms)
  miss <- setdiff(gfs, names(newdata))
  if (length(miss))
    stop("countimp: newdata is missing the grouping factor(s) ",
         paste(sQuote(miss), collapse = ", "),
         " required by the fitted model.", call. = FALSE)

  ## suppressWarnings: sdreport() warns about non-positive-definite Hessians,
  ## which is precisely the condition we detect below and handle by falling
  ## back. Letting the warning through would emit it once per imputation.
  P <- suppressWarnings(tryCatch(
    TMB::sdreport(fit$obj, getJointPrecision = TRUE)$jointPrecision,
    error = function(e) NULL))
  if (is.null(P)) return(fallback("sdreport(getJointPrecision) failed"))
  nmz <- colnames(as.matrix(P))
  if (is.null(nmz) || !any(nmz == "beta"))
    return(fallback("joint precision has no 'beta' block"))

  ## A non-positive-definite P is an expected, handled condition (it happens
  ## when the fit did not converge, e.g. an overparameterised random-effect
  ## structure). CHOLMOD emits a warning on top of the error; suppressing it
  ## keeps one handled condition from producing dozens of warnings per
  ## imputation, which would bury the diagnostics the user should see.
  U <- suppressWarnings(tryCatch(Matrix::chol(P), error = function(e) NULL))
  if (is.null(U))
    return(fallback("joint precision is not positive definite (fit may not have converged)"))
  z <- tryCatch(as.numeric(Matrix::solve(U, stats::rnorm(nrow(P)))),
                error = function(e) NULL)
  if (is.null(z) || !all(is.finite(z)))
    return(fallback("solving the Cholesky system failed"))

  ## eta.hat, with a fresh effect for clusters the fit never saw
  eta <- as.numeric(.countimp_predict_2l(fit, newdata, type = "link",
                                         grp = grp, obs_levels = obs_levels))
  X <- tryCatch(stats::model.matrix(stats::delete.response(
                  stats::terms(stats::formula(fit, component = "cond", fixed.only = TRUE))),
                  newdata), error = function(e) NULL)
  ib <- which(nmz == "beta")
  if (is.null(X) || ncol(X) != length(ib))
    return(fallback("fixed-effect design for newdata could not be rebuilt"))
  eta <- eta + as.numeric(X %*% z[ib])

  Z  <- .countimp_re_design(fit, newdata)
  iu <- which(nmz == "b")
  if (!is.null(Z) && length(iu) == ncol(Z))
    eta <- eta + as.numeric(Z %*% z[iu])
  else if (length(iu))
    return(fallback("random-effect design for newdata could not be rebuilt"))

  ## Dispersion from the same draw; exp(betadisp) is theta in nbinom2.
  ## An explicitly supplied theta is used as given and NOT perturbed: the
  ## argument exists to override the fit's dispersion (e.g. a value already
  ## drawn elsewhere, or a fixed one), so drawing around it again would apply
  ## the uncertainty twice.
  if (!is.null(theta))
    return(list(eta = eta, mu = as.numeric(linkinv(eta)), theta = theta,
                joint = TRUE))
  th <- suppressWarnings(as.numeric(.countimp_theta(fit, pflicht = FALSE) %||% NA_real_))
  idisp <- which(nmz == "betadisp")
  if (length(idisp) == 1L && is.finite(th) && th > 0) {
    th <- exp(log(th) + z[idisp])
    if (th < theta.min) {
      .countimp_note_event("theta_draw_floored",
                       sprintf("joint draw theta = %.4g raised to %.4g", th, theta.min))
      th <- theta.min
    }
  }

  list(eta = eta, mu = as.numeric(linkinv(eta)), theta = th, joint = TRUE)
}

## Joint draw for a genuine zero-inflation fit -------------------------------
##
## For a single glmmTMB fit with a ziformula, one draw from the Laplace
## posterior must move BOTH linear predictors -- the conditional count part and
## the zero-inflation part -- because they are estimated jointly and are
## correlated. Drawing them from two separate fits (as the hurdle methods do,
## legitimately, because there the two parts ARE separate likelihoods) would
## discard that correlation.
##
## The joint precision carries the parts in named blocks:
##   beta (count fixed) | b (count random) | betazi | bzi | betadisp | theta...
## `bzi` appears only when the ziformula has a random part; when it does, the
## `b` block still refers to the count part alone, so the count-side indexing
## in .countimp_draw_2l() stays valid.
##
## Returns pi (zero-inflation probability), mu (count rate) and theta from ONE
## draw. On any failure it falls back to the two marginal Wald draws and says
## so through .countimp_note_event(), never silently.
.countimp_draw_2l_zi <- function(fit, newdata, grp = NULL, obs_levels = NULL,
                                 theta.min = 0.25) {
  if (!inherits(fit, "glmmTMB"))
    stop("countimp: .countimp_draw_2l_zi() expects a glmmTMB fit, got ",
         paste(class(fit), collapse = "/"), ".", call. = FALSE)

  eta.c <- as.numeric(.countimp_predict_2l(fit, newdata, type = "link",
                                           grp = grp, obs_levels = obs_levels))
  eta.z <- as.numeric(stats::predict(fit, newdata = newdata, type = "zlink",
                                     allow.new.levels = TRUE))
  th    <- suppressWarnings(as.numeric(.countimp_theta(fit, pflicht = FALSE) %||% NA_real_))

  fallback <- function(why) {
    if (!isTRUE(.countimp_state$joint_zi_noted)) {
      .countimp_state$joint_zi_noted <- TRUE
      .countimp_note_event("joint_zi_draw_unavailable", why)
    }
    ## marginal Wald draws, one per component
    pert <- function(comp, eta) {
      b <- glmmTMB::fixef(fit)[[comp]]
      V <- try(stats::vcov(fit)[[comp]], silent = TRUE)
      X <- try(stats::model.matrix(stats::delete.response(stats::terms(
             stats::formula(fit, component = if (comp == "cond") "cond" else "zi",
                            fixed.only = TRUE))), newdata), silent = TRUE)
      if (inherits(V, "try-error") || inherits(X, "try-error") ||
          is.null(V) || ncol(X) != length(b)) return(eta)
      rv <- try(t(chol(V)), silent = TRUE)
      if (inherits(rv, "try-error")) return(eta)
      d <- as.numeric(rv %*% stats::rnorm(ncol(rv)))
      eta + as.numeric(X %*% d)
    }
    list(mu = exp(pert("cond", eta.c)),
         pi = stats::plogis(pert("zi", eta.z)),
         theta = .countimp_draw_theta(fit), joint = FALSE)
  }

  P <- suppressWarnings(tryCatch(
    TMB::sdreport(fit$obj, getJointPrecision = TRUE)$jointPrecision,
    error = function(e) NULL))
  if (is.null(P)) return(fallback("sdreport(getJointPrecision) failed"))
  nmz <- colnames(as.matrix(P))
  if (is.null(nmz) || !any(nmz == "beta") || !any(nmz == "betazi"))
    return(fallback("joint precision lacks a 'beta' or 'betazi' block"))

  U <- suppressWarnings(tryCatch(Matrix::chol(P), error = function(e) NULL))
  if (is.null(U))
    return(fallback("joint precision is not positive definite (fit may not have converged)"))
  z <- tryCatch(as.numeric(Matrix::solve(U, stats::rnorm(nrow(P)))),
                error = function(e) NULL)
  if (is.null(z) || !all(is.finite(z)))
    return(fallback("solving the Cholesky system failed"))

  dm <- function(comp) tryCatch(stats::model.matrix(stats::delete.response(
          stats::terms(stats::formula(fit, component = comp, fixed.only = TRUE))),
          newdata), error = function(e) NULL)
  Xc <- dm("cond"); Xz <- dm("zi")
  ib <- which(nmz == "beta"); iz <- which(nmz == "betazi")
  if (is.null(Xc) || ncol(Xc) != length(ib))
    return(fallback("count-part design for newdata could not be rebuilt"))
  if (is.null(Xz) || ncol(Xz) != length(iz))
    return(fallback("zero-part design for newdata could not be rebuilt"))
  eta.c <- eta.c + as.numeric(Xc %*% z[ib])
  eta.z <- eta.z + as.numeric(Xz %*% z[iz])

  ## count-side random effects; `b` excludes the zi random part (see above)
  Z  <- .countimp_re_design(fit, newdata)
  iu <- which(nmz == "b")
  if (!is.null(Z) && length(iu) == ncol(Z))
    eta.c <- eta.c + as.numeric(Z %*% z[iu])
  else if (length(iu))
    return(fallback("random-effect design for newdata could not be rebuilt"))

  idisp <- which(nmz == "betadisp")
  if (length(idisp) == 1L && is.finite(th) && th > 0) {
    th <- exp(log(th) + z[idisp])
    if (th < theta.min) {
      .countimp_note_event("theta_draw_floored",
        sprintf("joint ZI draw theta = %.4g raised to %.4g", th, theta.min))
      th <- theta.min
    }
  }

  list(mu = exp(eta.c), pi = stats::plogis(eta.z), theta = th, joint = TRUE)
}
