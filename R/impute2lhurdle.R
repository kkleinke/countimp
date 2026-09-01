## impute2lhurdle.R -- one engine for the two-level hurdle methods
##
## Up to 3.0.0 the 16 exported names in mice.impute.2l.hnb.R were 16 copies of
## the same ~118-line body (pairwise similarity 0.98), differing in four
## places: the truncated count family (truncated_nbinom2 / truncated_poisson),
## the draw (Bayes / cluster bootstrap), and the two intercept flags -- which
## are already function arguments, so the .noint variants only fix their
## defaults. That is the structure that carried defect B30 in 16 copies of the
## zero-inflated family, and the reason for this rebuild.
##
## WHY A HURDLE IS FITTED IN TWO PIECES AND ZERO INFLATION IS NOT. The hurdle
## likelihood factorises: the gate is a binomial GLMM on 1{Y = 0} over ALL
## units, and the positive part is a zero-truncated count GLMM over the
## POSITIVE units only. The two factors share no parameters, so two separate
## fits are exact, not an approximation -- and each part may be drawn from its
## own posterior independently. Zero inflation does NOT factorise (a zero can
## come from either component), which is why .countimp_2l_zi() needs one joint
## fit and one joint draw; see impute2lzi.R.
##
## The consequence for the draw: a unit that passes the gate is imputed from a
## ZERO-TRUNCATED law and can never be zero. That is exactly what distinguishes
## these methods from 2l.zip / 2l.zinb, and it was the defect in B30 that the
## zero-inflated methods drew untruncated counts behind a gate.
##
## Validated against the pre-rebuild implementation by
## skripte/k01_referenz_hurdle.R (distribution of the imputations per method:
## mean, zero proportion, dispersion), because consolidation reorders the RNG
## stream and draws cannot be compared value by value.

## The engine. `family` is "nbinom2" (hnb) or "poisson" (hp); `draw` is
## "bayes" or "boot".
.countimp_2l_hurdle <- function(y, ry, x, type,
                                family = c("nbinom2", "poisson"),
                                draw = c("bayes", "boot"),
                                intercept.c = TRUE, intercept.z = TRUE,
                                wy = NULL, EV = FALSE) {
  family <- match.arg(family)
  draw   <- match.arg(draw)
  if (is.null(wy)) wy <- !ry

  X   <- data.frame(x[ry, , drop = FALSE])
  nam <- colnames(X)
  dec <- .countimp_decode_type(type, nam)
  dat <- data.frame(Y = y[ry], X)

  ## Resample whole clusters once, for BOTH parts: the gate and the truncated
  ## count model must see the same bootstrap sample, or the two parts describe
  ## different populations. Read the cluster levels after resampling -- cluster
  ## sampling leaves roughly a third of the clusters out of any one sample, and
  ## those have no u_j in either fit (see predict2l.R).
  ## nz is the response of the gate formula, so BOTH the resample fit and the
  ## full-data refit need it. It is row-wise, so computing it before the
  ## resample is equivalent and keeps the two in step -- computed after, as it
  ## was, dat.voll lacks the column and the gate refit fails into the fallback
  ## without anyone noticing. Measured: 1 of 2 refits silently degraded.
  dat$nz <- as.numeric(dat$Y == 0)
  dat.voll <- dat
  if (identical(draw, "boot")) dat <- .countimp_boot_clusters(dat, dec$group)
  ## from the FULL data, not the resample: a left-out cluster still has data
  ## and keeps its conditional posterior (see .countimp_boot_rate())
  obs.levels <- .countimp_obs_levels(dat.voll, dec$group)

  newdata <- data.frame(x[wy, , drop = FALSE])
  colnames(newdata) <- nam

  ## ---- part 1: the gate, on all units -------------------------------------
  ## Response is 1{Y = 0}, so the fitted probability is P(zero) and a unit is
  ## imputed positive with probability 1 - p. intercept.z and the type codes
  ## 5/6 act on this part alone.
  f.gate <- .countimp_2l_formula(dec, "zero", response = "nz",
                                 intercept = intercept.z)
  fit.gate <- glmmTMB::glmmTMB(formula = f.gate, data = dat,
                               family = .countimp_family("binomial"))
  p.zero <- if (identical(draw, "bayes")) {
    .countimp_draw_2l(fit.gate, newdata, grp = dec$group,
                      obs_levels = obs.levels)$mu
  } else {
    ## The bootstrap resample carries the parameter uncertainty for the gate.
    ## .countimp_rate() applies the fit's own inverse link (plogis here), so
    ## the binomial link is not hard-coded.
    .countimp_boot_rate(fit.gate, f.gate, dat.voll,
                        .countimp_family("binomial"), newdata,
                        grp = dec$group, obs_levels = obs.levels)
  }
  positive <- stats::rbinom(length(p.zero), size = 1L, prob = 1 - p.zero) == 1L

  ## ---- part 2: the truncated count model, on the positive units -----------
  ## intercept.c and the type codes 3/4 act on this part alone.
  pos.dat <- dat[dat$Y > 0, , drop = FALSE]
  f.cnt <- .countimp_2l_formula(dec, "count", response = "Y",
                               intercept = intercept.c)
  tfam <- if (identical(family, "poisson")) "truncated_poisson" else "truncated_nbinom2"
  fit.cnt <- glmmTMB::glmmTMB(formula = f.cnt, data = pos.dat,
                              family = .countimp_family(tfam, "log"))
  ## theta is often weakly identified in truncated models
  if (identical(family, "nbinom2"))
    fit.cnt <- .countimp_stabilize(fit.cnt, f.cnt, pos.dat)$fit

  ## NOTE on type = "link" / $mu here. For a truncated family glmmTMB's
  ## conditional mean is E[Y | Y > 0], which is NOT the rate parameter the
  ## truncated draw needs (measured: 2.198 vs 1.446). .countimp_draw_2l()
  ## returns the linear predictor and its inverse link, i.e. the rate -- which
  ## is the right quantity. Do not substitute type = "conditional".
  cnt.dr <- if (identical(draw, "bayes")) {
    .countimp_draw_2l(fit.cnt, newdata, grp = dec$group,
                      obs_levels = obs.levels)
  } else {
    ## The cluster bootstrap carries the uncertainty in the regression
    ## parameters, but NOT in theta -- a resample yields one point estimate per
    ## fit. .countimp_draw_theta() adds the posterior draw of theta and
    ## enforces the lower bound that guards against the truncated-NB
    ## imputation explosions recorded as B09a. Using sigma() directly would
    ## drop both. .countimp_rate() inverts the model's own link instead of
    ## hard-coding exp().
    ## the count part is fitted on the positive units, so its own "full
    ## data" is dat.voll restricted the same way
    list(mu = .countimp_boot_rate(fit.cnt, f.cnt,
                                  dat.voll[dat.voll$Y > 0, , drop = FALSE],
                                  .countimp_family(tfam, "log"), newdata,
                                  grp = dec$group, obs_levels = obs.levels),
         theta = .countimp_draw_theta(fit.cnt))
  }

  ## Zero-truncated draw: a unit past the gate is positive by construction.
  n <- length(cnt.dr$mu)
  cnt <- if (identical(family, "poisson")) {
    .countimp_rktp(n = n, k = 0, mu = cnt.dr$mu, xpred = 1)
  } else {
    .countimp_rktnb(n = n, size = cnt.dr$theta, k = 0, mu = cnt.dr$mu,
                    xpred = 1)
  }

  ## numeric, not integer: the truncated draws come back as double, and
  ## assigning them into an integer vector would coerce silently.
  imp <- numeric(length(positive))
  imp[positive] <- cnt[positive]

  if (EV) imp <- .countimp_ev_handle(imp, y, x, wy)
  imp
}
