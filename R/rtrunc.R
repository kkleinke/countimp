## rtrunc.R -- k-truncated Poisson and negative binomial draws
##
## Why this file exists. The hurdle methods (2l.hp, 2l.hnb) draw the positive
## part of the outcome from a zero-truncated count distribution. Until now that
## came from aster::rktp() and aster::rktnb(), which made aster a hard Depends:
## 3.1 MB of compiled code, R (>= 4.3.0), plus the package 'trust' -- for two
## random number generators. countimp used 2 of aster's 30 exports, and only
## ever with k = 0 and xpred = 1, i.e. plain zero truncation.
##
## The dependency is therefore replaced by base R. That removes a compiled
## package from the install chain, drops the R >= 4.3.0 floor that aster
## imposes, and eliminates a source of breakage countimp cannot control.
##
## Method. Naive rejection sampling ("draw until the value exceeds k") is the
## obvious approach and the wrong one: when mu is small, P(X = 0) approaches 1
## and the expected number of rejected draws per accepted one explodes -- at
## mu = 0.01 it is about 100. Hurdle models routinely produce small mu for the
## count part, so this matters in practice.
##
## Instead we invert the truncated cdf. For X ~ F truncated to X > k,
##
##     F_trunc(x) = (F(x) - F(k)) / (1 - F(k))
##
## so a draw is X = F^-1( F(k) + U * (1 - F(k)) ) with U ~ Unif(0, 1). Written
## that way it loses precision when F(k) is close to 1, which is exactly the
## small-mu case. Rewriting in terms of the survival function,
##
##     F(k) + U(1 - F(k)) = 1 - S(k)(1 - U),      S(k) = P(X > k)
##
## and passing that upper-tail probability to qpois()/qnbinom() with
## lower.tail = FALSE avoids the cancellation entirely. The result is exact
## (not approximate), needs one uniform per draw regardless of mu, and is
## vectorised over mu.
##
## xpred, in aster's parameterisation, is the number of i.i.d. truncated
## variates summed per observation. countimp only ever passes xpred = 1; the
## general case is implemented anyway so the functions are drop-in replacements.

.countimp_rktp <- function(n, k = 0, mu, xpred = 1) {
  mu <- rep_len(as.numeric(mu), n)
  if (any(!is.finite(mu)) || any(mu < 0))
    stop("mu must be finite and non-negative", call. = FALSE)
  xp <- rep_len(as.integer(xpred), n)
  if (any(xp < 1L)) stop("xpred must be a positive integer", call. = FALSE)

  draw1 <- function(mu.i) {
    s <- stats::ppois(k, mu.i, lower.tail = FALSE)          # P(X > k)
    ## S(k) == 0 to machine precision: the truncation region is numerically
    ## unreachable (mu far below k). Return the smallest admissible value
    ## rather than NaN -- the model is degenerate there in any case.
    out <- rep(k + 1L, length(mu.i))
    ok  <- s > 0
    if (any(ok))
      out[ok] <- stats::qpois(stats::runif(sum(ok)) * s[ok], mu.i[ok],
                              lower.tail = FALSE)
    out
  }

  if (all(xp == 1L)) return(as.numeric(draw1(mu)))
  vapply(seq_len(n), function(i) sum(draw1(rep(mu[i], xp[i]))), numeric(1))
}

.countimp_rktnb <- function(n, size, k = 0, mu, xpred = 1) {
  mu <- rep_len(as.numeric(mu), n)
  if (any(!is.finite(mu)) || any(mu < 0))
    stop("mu must be finite and non-negative", call. = FALSE)
  if (!is.finite(size) || size <= 0)
    stop("size (theta) must be finite and positive", call. = FALSE)
  xp <- rep_len(as.integer(xpred), n)
  if (any(xp < 1L)) stop("xpred must be a positive integer", call. = FALSE)

  draw1 <- function(mu.i) {
    s <- stats::pnbinom(k, size = size, mu = mu.i, lower.tail = FALSE)
    out <- rep(k + 1L, length(mu.i))
    ok  <- s > 0
    if (any(ok))
      out[ok] <- stats::qnbinom(stats::runif(sum(ok)) * s[ok], size = size,
                                mu = mu.i[ok], lower.tail = FALSE)
    out
  }

  if (all(xp == 1L)) return(as.numeric(draw1(mu)))
  vapply(seq_len(n), function(i) sum(draw1(rep(mu[i], xp[i]))), numeric(1))
}


## ---------------------------------------------------------------------------
## .countimp_rate() -- the RATE parameter, not the truncated mean
##
## glmmTMB's predict(type = "response") returns the mean of the *fitted*
## distribution. For a truncated family that mean is already conditional on
## Y > 0:
##
##     truncated_poisson:  E[Y | Y>0] = mu / (1 - exp(-mu))
##     truncated_nbinom2:  E[Y | Y>0] = mu / (1 - (k/(k+mu))^k)
##
## The zero-truncated samplers, however, expect the *untruncated* rate mu and
## apply the truncation themselves. Passing predict(type = "response") into
## them therefore truncates twice and inflates the imputations: at mu = 1 the
## drawn mean is 1.99 instead of 1.58 (+26 %), at mu = 0.2 it is 1.65 instead
## of 1.10 (+50 %). The error vanishes for large mu, where truncation stops
## mattering -- which is why it survives casual inspection on well-behaved
## data and bites hardest on the sparse, zero-heavy counts these methods
## exist for.
##
## The rate is recovered from the link scale, which predict() leaves alone.
## ---------------------------------------------------------------------------

.countimp_rate <- function(fit, newdata, grp = NULL, obs_levels = NULL) {
  eta <- if (exists(".countimp_predict_2l", mode = "function")) {
    .countimp_predict_2l(fit, newdata, type = "link",
                         grp = grp, obs_levels = obs_levels)
  } else {
    stats::predict(fit, newdata = newdata, type = "link",
                   na.action = stats::na.pass, allow.new.levels = TRUE)
  }
  linkinv <- tryCatch(stats::family(fit)$linkinv,
                      error = function(e) exp)
  if (!is.function(linkinv)) linkinv <- exp
  as.numeric(linkinv(eta))
}


## ---------------------------------------------------------------------------
## .countimp_stabilize() -- catch a non-identified theta
##
## In truncated_nbinom2 models theta is often poorly identified: the truncated
## distribution carries little information about overdispersion, because the
## zeros -- the most informative part -- are absent. When theta collapses
## towards 0 the profile likelihood in beta goes flat, SE(beta) jumps by orders
## of magnitude, and the subsequent draw b* ~ N(beta, V) produces linear
## predictors beyond 700, i.e. mu = Inf.
##
## Remedy: fix theta at a moment estimator computed from the observed positive
## values (map = betadisp) and re-estimate beta conditional on it. The
## intervention runs only when max SE(beta) > se.max, then costs about 0.03 s,
## and is kept only if it actually reduces SE(beta). Otherwise the original fit
## is returned unchanged.
##
## Returns a list with
##   fit          -- the fit (re-estimated where applicable)
##   theta.fixed  -- TRUE if theta was fixed
##   theta        -- the theta that was used
## ---------------------------------------------------------------------------
.countimp_mom_theta <- function(yp, lower = 0.1, upper = 50) {
  yp <- yp[is.finite(yp) & yp > 0]
  if (length(yp) < 5L) return(upper)
  m <- mean(yp); v <- stats::var(yp)
  if (!is.finite(v) || v <= m) return(upper)   # no detectable overdispersion
  min(max(m^2 / (v - m), lower), upper)
}

.countimp_max_se <- function(fit) {
  V <- try(as.matrix(stats::vcov(fit)$cond), silent = TRUE)
  if (inherits(V, "try-error") || !length(V)) return(NA_real_)
  s <- sqrt(diag(V))
  if (!any(is.finite(s))) return(Inf)
  max(s[is.finite(s)])
}

.countimp_stabilize <- function(fit, form, data, se.max = 5, theta.min = 0.25,
                                quiet = FALSE) {
  out <- list(fit = fit, theta.fixed = FALSE, theta = NA_real_)
  ## stats::sigma() and not glmmTMB::sigma(): sigma() is a generic in stats
  ## since R 3.3, and the qualified generic keeps working through S3 dispatch
  ## even if glmmTMB stops re-exporting it. Callers must pass a glmmTMB
  ## nbinom2 fit -- for an MASS::glm.nb object stats::sigma() returns the
  ## residual scale (1), not theta.
  th <- try(stats::sigma(fit), silent = TRUE)
  out$theta <- if (inherits(th, "try-error")) NA_real_ else as.numeric(th)

  ## Two independent symptoms of a non-identified theta:
  ##  (a) inflated SE(beta) -- the flat-likelihood route, blows up the beta draw
  ##  (b) theta below theta.min -- the heavy-tail route, blows up the draw itself
  ##      (E[Y|Y>0] grows without bound as theta -> 0 at fixed mu)
  ## Detected separately because either can occur without the other: in the
  ## validation runs (a) caught 2 of 4 explosions, (b) caught all 4.
  se0 <- .countimp_max_se(fit)
  bad.se <- !is.finite(se0) || se0 > se.max
  bad.th <- !is.finite(out$theta) || out$theta < theta.min
  if (!bad.se && !bad.th) return(out)

  tm  <- .countimp_mom_theta(data[["Y"]])
  fam <- try(stats::family(fit), silent = TRUE)
  if (inherits(fam, "try-error")) return(out)

  f2 <- try(suppressWarnings(glmmTMB::glmmTMB(
          formula = form, data = data, family = fam,
          start = list(betadisp = log(tm)),
          map   = list(betadisp = factor(NA)))), silent = TRUE)
  if (inherits(f2, "try-error")) return(out)

  ## Accept only if the refit is a genuine improvement.  With a pure theta
  ## alarm SE(beta) may already be small, so requiring a *smaller* SE would
  ## reject every such refit; there we require that it not get worse.
  se1 <- .countimp_max_se(f2)
  if (!is.finite(se1)) return(out)
  if (bad.se && se1 >= se0) return(out)
  if (!bad.se && se1 > 2 * se0) return(out)

  if (!quiet)
    warning(sprintf(paste0("countimp: theta was not identified (theta = %.3g, ",
            "max SE(beta) = %.3g; trigger: %s). Refitted with theta fixed at ",
            "the moment estimate %.3g; max SE(beta) is now %.3g. The ",
            "imputation model conditions on this value, so between-imputation ",
            "variability in theta is not reflected for this fit."),
            out$theta, se0,
            paste(c("SE(beta)", "theta")[c(bad.se, bad.th)], collapse = " + "),
            tm, se1), call. = FALSE)

  list(fit = f2, theta.fixed = TRUE, theta = tm)
}


## Quasi-Poisson draws --------------------------------------------------------
##
## A quasi-Poisson GLM estimates a mean function plus a dispersion phi with
## Var(Y) = phi * mu. There is no quasi-Poisson random number generator,
## because a quasi-likelihood is not a distribution -- so imputations must be
## drawn from a real distribution matched to the fitted mean and variance.
## countimp uses the negative binomial with
##
##     size = mu / (phi - 1),   since  Var = mu + mu^2/size = phi * mu.
##
## That identity holds only for phi > 1. The Bayesian variants applied it
## unconditionally, with two distinct consequences:
##
##   phi < 1  -> size < 0   -> rnbinom() returns NaN for EVERY imputed value
##   phi = 1  -> size = Inf -> rnbinom() is exactly Poisson. Verified against
##                             rpois (KS test p = 0.29), so this is correct
##                             behaviour, not a defect.
##
## Underdispersion is not an edge case. It arises whenever counts have an
## upper bound -- a six-item symptom count, days per week, doctor visits per
## month -- and phi is a sample quantity that can fall below 1 even when the
## truth is Poisson. Measured end to end: a binomial DGP (mean 3.48, variance
## 2.04) gave a fitted phi of 0.427, and 270 of 270 imputed values were NaN.
##
## No negative binomial can represent Var < mu, so within that family there is
## nothing to draw from. The fallback is the Poisson (Var = mu), which
## overstates the variance of the imputations and so errs towards intervals
## that are too wide rather than towards false precision. The bootstrap
## variants of this package already branched exactly this way; the Bayesian
## ones never got the branch.
##
## Non-finite mu is treated as an error, not silently passed on: rnbinom() and
## rpois() both answer a non-finite mean with NA, and an NA imputation written
## back into the data corrupts the analysis silently, which is worse than
## stopping.

.countimp_warn_underdisp <- function(dispersion) {
  if (isTRUE(.countimp_state$underdisp_warned)) return(invisible(FALSE))
  .countimp_state$underdisp_warned <- TRUE
  warning(sprintf(paste0(
    "countimp: the quasi-Poisson dispersion was estimated as %.3g, i.e. the ",
    "counts are underdispersed (variance below the mean). The negative ",
    "binomial used for quasi-Poisson imputation cannot represent that -- ",
    "size = mu/(phi - 1) would be negative and every imputed value would be ",
    "NaN. Imputations were drawn from a Poisson instead, which overstates ",
    "their variance rather than understating it. Consider whether a bounded ",
    "outcome would be better modelled directly (e.g. as a binomial ",
    "proportion). This warning is shown once per session."), dispersion),
    call. = FALSE)
  invisible(TRUE)
}

.countimp_rqpois <- function(mu, dispersion, quiet = FALSE) {
  mu <- as.numeric(mu)
  n  <- length(mu)
  if (n == 0L) return(numeric(0))
  if (!is.numeric(dispersion) || length(dispersion) != 1L || is.na(dispersion))
    stop("countimp: the quasi-Poisson dispersion estimate is missing or not ",
         "a single number. The imputation model could not be fitted.",
         call. = FALSE)
  bad <- !is.finite(mu) | mu < 0
  if (any(bad))
    stop(sprintf(paste0(
      "countimp: %d of %d fitted means are not finite and non-negative, so ",
      "the draw would return NaN. exp(x %%*%% beta) overflowed, which ",
      "usually means the imputation model is separated or a predictor is on ",
      "a very large scale. Rescale or drop collinear predictors."),
      sum(bad), n), call. = FALSE)
  if (dispersion > 1)
    return(stats::rnbinom(n = n, size = mu/(dispersion - 1), mu = mu))
  if (!quiet) .countimp_warn_underdisp(dispersion)
  stats::rpois(n, mu)
}


## Draws from zero-inflated and hurdle models ---------------------------------
##
## The single-level zi/hurdle methods used to impute like this:
##
##     pc <- predict(fit, newdata, type = "prob")
##     for (i in 1:nrow(pc))
##       pcvec[i] <- sample(as.numeric(names(pc[i, ])), 1, pc[i, ])
##
## i.e. they enumerated the probability mass function on a grid and sampled a
## grid point. The grid is the problem. pscl's predict(type = "prob") returns
## one column per value 0, 1, ..., max(y) of the TRAINING data, so the draw is
## confined to the range of the observed counts. Anything the model predicts
## beyond that range is unreachable, and the mass that belongs there is simply
## missing from the row (rows do not sum to 1 when extrapolating).
##
## Measured: with a training maximum of 48 and a case whose fitted count mean
## was 134.7, the support stopped at 48 and 72.9% of the probability mass was
## cut off. Imputations are therefore systematically too small precisely for
## the cases that need large values, which biases the imputation model towards
## the observed range -- the opposite of what multiple imputation is for.
##
## The replacement draws from the distribution itself. Both model families
## factor into the same two parts,
##
##     P(Y = 0) = p0,      P(Y = k | Y > 0) = zero-truncated count density,
##
## which was verified against predict(type = "prob") for all four
## combinations (zeroinfl/hurdle x poisson/negbin): the reconstructed p0 and
## the conditional positive distribution agree to machine precision
## (max deviation 1.4e-16). The two families differ only in how p0 is built:
##
##   zeroinfl: p0 = pi + (1 - pi) * f(0)      pi = P(structural zero)
##   hurdle  : p0 = 1 - z * (1 - f(0))        z  = predict(type = "zero")
##
## Two pscl pitfalls this code deliberately avoids:
##
##   * type = "response" is NOT the count mean. For zeroinfl it is
##     (1 - pi) * mu, for hurdle it is z * mu. Feeding it to a count
##     distribution builds the zeros in twice. Use type = "count".
##   * for hurdle, type = "zero" is a RATIO, not a probability:
##     z = P(Y > 0) / P(count part > 0). Verified numerically. Treating it as
##     P(zero) inverts the zero part.
##
## Robustness. The identities above are pscl semantics, and pscl could change
## them. Rather than trust them silently, .countimp_rzi() re-derives p0 and
## compares it against predict(type = "prob") for a single row on every call.
## A semantic change in pscl then produces an explicit error naming the
## deviation instead of silently wrong imputations. One extra one-row predict
## call is a small price for that.

.countimp_zi_theta <- function(fit) {
  ## Inf for a Poisson count part, else the count-part theta.
  d <- fit$dist
  d <- if (is.list(d)) d[["count"]] else d[[1L]]
  if (is.null(d) || !identical(as.character(d)[1L], "negbin")) return(Inf)
  ## Field knowledge lives in .countimp_theta() (see foreign-packages.R); the
  ## message below stays here because it is about the imputation, not the
  ## package. pflicht = FALSE so a missing theta reaches that message rather
  ## than the accessor's generic one.
  th <- unlist(.countimp_theta(fit, pflicht = FALSE), use.names = TRUE)
  if (length(th) == 0L)
    stop("countimp: the fitted model reports a negative binomial count part ",
         "but carries no theta. Cannot draw imputations.", call. = FALSE)
  th <- if (!is.null(names(th)) && "count" %in% names(th)) th[["count"]] else th[[1L]]
  th <- as.numeric(th)
  if (!is.finite(th) || th <= 0)
    stop(sprintf(paste0("countimp: the estimated theta of the count part is ",
         "%s, which is not a valid dispersion. The imputation model did not ",
         "converge."), format(th)), call. = FALSE)
  th
}

.countimp_zi_p0 <- function(kind, z, mu, theta) {
  f0 <- if (is.finite(theta)) stats::dnbinom(0, size = theta, mu = mu)
        else                  stats::dpois(0, mu)
  if (identical(kind, "zeroinfl")) z + (1 - z) * f0 else 1 - z * (1 - f0)
}

## Kept separate from .countimp_rzi() so the guard itself is testable without
## having to fake a pscl object: feeding a fit a different class makes pscl's
## own predict() switch families too, so both sides move together and the
## guard cannot be exercised through the fit. Called with plain numbers, it can.
.countimp_zi_check <- function(p0.own, p0.pscl, kind, tol = 1e-6) {
  if (!is.finite(p0.pscl) || !is.finite(p0.own)) return(invisible(TRUE))
  if (abs(p0.pscl - p0.own) <= tol) return(invisible(TRUE))
  pv <- tryCatch(as.character(utils::packageVersion("pscl")),
                 error = function(e) "unknown")
  stop(sprintf(paste0(
    "countimp: P(Y = 0) reconstructed from predict(type = 'count'/'zero') is ",
    "%.8f, but pscl's predict(type = 'prob') reports %.8f for the same case ",
    "(%s model). countimp assumes zeroinfl: p0 = pi + (1 - pi) f(0) and ",
    "hurdle: p0 = 1 - z (1 - f(0)). If pscl changed these semantics the ",
    "imputations would be wrong, so countimp stops rather than continue. ",
    "Please report this together with your pscl version (%s)."),
    p0.own, p0.pscl, kind, pv), call. = FALSE)
}

.countimp_rzi <- function(fit, newdata) {
  kind <- if (inherits(fit, "zeroinfl")) "zeroinfl"
          else if (inherits(fit, "hurdle")) "hurdle"
          else stop("countimp: .countimp_rzi() expects a zeroinfl or hurdle ",
                    "fit, got class ", paste(class(fit), collapse = "/"), ".",
                    call. = FALSE)
  n <- nrow(newdata)
  if (is.null(n) || n == 0L) return(numeric(0))

  mu <- as.numeric(stats::predict(fit, newdata = newdata, type = "count",
                                  na.action = stats::na.pass))
  z  <- as.numeric(stats::predict(fit, newdata = newdata, type = "zero",
                                  na.action = stats::na.pass))
  bad <- !is.finite(mu) | mu < 0 | !is.finite(z)
  if (any(bad))
    stop(sprintf(paste0(
      "countimp: %d of %d predicted values from the %s model are not finite ",
      "(count mean or zero part). exp(x %%*%% beta) overflowed, which usually ",
      "means the imputation model is separated or a predictor is on a very ",
      "large scale. Rescale or drop collinear predictors."), sum(bad), n, kind),
      call. = FALSE)

  theta <- .countimp_zi_theta(fit)
  p0 <- .countimp_zi_p0(kind, z, mu, theta)

  ## Guard against a change in pscl's predict() semantics: re-derive P(Y = 0)
  ## for one row and compare with pscl's own probability matrix.
  ref <- try(stats::predict(fit, newdata = newdata[1L, , drop = FALSE],
                            type = "prob", na.action = stats::na.pass),
             silent = TRUE)
  if (!inherits(ref, "try-error"))
    .countimp_zi_check(p0[1L], suppressWarnings(as.numeric(ref[1L, 1L])), kind)

  p0 <- pmin(pmax(p0, 0), 1)
  out <- numeric(n)
  pos <- stats::runif(n) >= p0
  if (any(pos))
    out[pos] <- if (is.finite(theta))
      .countimp_rktnb(sum(pos), size = theta, k = 0, mu = mu[pos])
    else
      .countimp_rktp(sum(pos), k = 0, mu = mu[pos])
  out
}


## Drawing the dispersion parameter (finding B09) ---------------------------
##
## Proper multiple imputation draws every unknown parameter, not just the
## regression coefficients. countimp up to 2.6.0 drew beta but passed the
## point estimate theta.hat to the count draw:
##
##     im <- rnegbin(length(p), theta = fit.sum$sigma, mu = p)
##
## Every imputation then conditions on the same dispersion. Rubin's
## between-imputation variance B captures no uncertainty about theta, so the
## total variance T = Wbar + (1 + 1/m) B is too small and pooled intervals for
## quantities that depend on the tail of the count distribution are too narrow.
## The effect is largest exactly where these models are used: small samples,
## weakly identified dispersion, high missingness.
##
## What is drawn. theta > 0, and all three backends report its uncertainty on
## the log scale, so the draw is made there:
##
##     log theta* ~ N(log theta.hat, se(log theta.hat))
##
## which keeps theta* positive and matches the parameterisation in which the
## asymptotic normal approximation actually holds.
##
## Where se(log theta) comes from, per backend:
##   glmmTMB      sqrt(diag(fit$sdr$cov.fixed))["betadisp"] -- betadisp is
##                log(theta), verified: exp(betadisp) == sigma(fit).
##   MASS::glm.nb fit$SE.theta / fit$theta (delta method, theta scale -> log).
##   pscl         fit$SE.logtheta (already on the log scale; note that
##                summary()$vcov deliberately excludes log(theta), which is
##                why the coefficient draw does not cover it -- see B17).
##
## When NOT to draw. If theta was fixed -- which .countimp_stabilize() does
## when theta is not identified -- there is no standard error to draw from,
## and the backends then report none: 'betadisp' disappears from cov.fixed.
## Poisson fits have no theta at all. In both cases se is NA and the point
## estimate is returned unchanged. The two mechanisms therefore compose
## correctly: a non-identified theta is fixed by the stabiliser and then not
## drawn, which is the honest outcome (and the stabiliser already warns that
## between-imputation variability in theta is not reflected for that fit).

.countimp_theta_se <- function(fit) {
  ## glmmTMB: log(theta) is the 'betadisp' entry of the joint precision
  if (inherits(fit, "glmmTMB")) {
    s <- try(sqrt(diag(fit$sdr$cov.fixed)), silent = TRUE)
    if (inherits(s, "try-error") || !length(s)) return(NA_real_)
    v <- s[names(s) == "betadisp"]
    if (!length(v) || !is.finite(v[1]) || v[1] <= 0) return(NA_real_)
    return(as.numeric(v[1]))
  }
  ## pscl::zeroinfl / pscl::hurdle report it directly on the log scale
  if (!is.null(fit$SE.logtheta)) {
    v <- as.numeric(fit$SE.logtheta)[1]
    return(if (is.finite(v) && v > 0) v else NA_real_)
  }
  ## MASS::glm.nb reports SE on the theta scale; delta method to log
  if (!is.null(fit$SE.theta) && !is.null(fit$theta)) {
    v <- as.numeric(fit$SE.theta)[1] / as.numeric(fit$theta)[1]
    return(if (is.finite(v) && v > 0) v else NA_real_)
  }
  NA_real_
}

.countimp_draw_theta <- function(fit, theta = NULL, theta.min = 0.25,
                                 se.max = 5, quiet = FALSE) {
  if (is.null(theta)) {
    theta <- try(as.numeric(stats::sigma(fit)), silent = TRUE)
    if (inherits(theta, "try-error")) theta <- NA_real_
  }
  ## Fail loudly on a non-scalar theta. Passing summary(fit) instead of
  ## summary(fit)$theta is an easy slip, and as.numeric() on a list would
  ## either error deep inside rnegbin() with an opaque message or -- worse --
  ## silently yield the wrong number.
  if (!is.numeric(theta) && !is.integer(theta))
    stop("countimp: .countimp_draw_theta() expects a numeric theta, got ",
         class(theta)[1], ". This is a bug in countimp -- please report it.",
         call. = FALSE)
  theta <- as.numeric(theta)[1]
  if (!is.finite(theta) || theta <= 0) return(theta)

  se <- .countimp_theta_se(fit)
  if (!is.finite(se)) {          # theta fixed, or backend reports no SE
    .countimp_note_theta_fixed(quiet)
    return(theta)
  }

  ## An UNIDENTIFIED theta is not drawn either, and it announces itself through
  ## se, not through theta.
  ##
  ## Fitting a negative binomial to data that carry no overdispersion -- the
  ## cautious choice, and therefore a frequent one -- sends theta towards
  ## infinity along a flat likelihood, and se(log theta) with it. The draw
  ## exp(rnorm(1, log(theta), se)) then spans orders of magnitude, and a draw
  ## near theta.min means Var = mu + 4 mu^2. Because every FCS cycle refits on
  ## the previous cycle's imputations, one unlucky draw feeds itself: measured
  ## on Poisson data with n = 800, one imputation of ten reached a standard
  ## deviation of 5.2 against 1.9 in the data, with a largest value of 44
  ## against 11 observed.
  ##
  ## The threshold is calibrated, not borrowed. .countimp_stabilize() uses 5
  ## for a different quantity (the largest SE among the coefficients), so the
  ## number carries over only by coincidence. Over 40 replications each of
  ## MASS::glm.nb, se(log theta) came out as
  ##
  ##   theta = 0.5            0.05 to 0.16    identified
  ##   theta = 2              0.06 to 0.29    identified
  ##   theta = 20, n >= 500   0.21 to 11      mostly identified
  ##   Poisson data           median 3.4 to 7.9, up to 32
  ##
  ## No identified case reaches 5, and at se = 5 the 95 % range of the draw
  ## already spans four orders of magnitude -- the draw then says less about
  ## the dispersion than the point estimate does.
  ##
  ## This composes with .countimp_stabilize(), which fixes theta in the
  ## multilevel paths: there se is NA afterwards and the branch above returns.
  ## The single-level path over MASS::glm.nb has no stabiliser, which is why
  ## the guard belongs here as well.
  if (se > se.max) {
    .countimp_note_theta_fixed(quiet)
    return(theta)
  }

  out <- exp(stats::rnorm(1L, log(theta), se))
  if (!is.finite(out) || out <= 0) return(theta)

  ## Guard the lower tail. As theta -> 0 at fixed mu the zero-truncated
  ## negative binomial mean grows without bound, which is the mechanism behind
  ## the imputation explosions recorded as B09a. A draw below theta.min is
  ## therefore raised to it. This is a deliberate, documented bias towards the
  ## bulk of the distribution: it trades a small distortion of the draw for
  ## imputations that stay on the scale of the data. It is rare when theta is
  ## identified -- and when it is not, .countimp_stabilize() has already fixed
  ## theta and this function does not draw at all.
  if (out < theta.min) {
    .countimp_note_event("theta_draw_floored",
                     sprintf("drawn theta = %.4g raised to %.4g", out, theta.min))
    out <- theta.min
  }
  out
}

.countimp_note_theta_fixed <- function(quiet = FALSE) {
  if (isTRUE(quiet) || isTRUE(.countimp_state$theta_fixed_noted))
    return(invisible(FALSE))
  .countimp_state$theta_fixed_noted <- TRUE
  .countimp_note_event("theta_not_drawn",
                   "no standard error for theta available; imputations condition on the point estimate")
  invisible(TRUE)
}
