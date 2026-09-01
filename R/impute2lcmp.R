## impute2lcmp.R -- two-level Conway-Maxwell-Poisson (underdispersion)
##
## The single-level case (R/compois.R, B82) is fitted by countimp itself. This
## one cannot be: the random effects need the Laplace approximation, which is
## why every 2l method in this package goes through glmmTMB. What remains
## countimp's own is everything that the register found fragile about that
## contact:
##
##   the DRAW   glmmTMB's simulate() ignores newdata and returns values for the
##              FITTING cases (measured in compois.R), which is unusable for
##              imputation. Draws come from .countimp_rcmp(), the exact
##              categorical sampler verified in test-B82.
##   the SHAPE  glmmTMB reports 1/nu through sigma(). Rather than believing
##              that, .countimp_cmp_nu_glmm() scores both readings against the
##              data with countimp's own likelihood and takes the better one.
##
## What is NOT reimplemented is the fit. Writing an adaptive quadrature over
## the random effects in R would mean one COM-Poisson normalising sum per
## quadrature node per case, and it would cover random intercepts only -- while
## the type codes of this package allow random slopes. A worse model that is
## ours is not better than a good model that is measured.
##
##
## THE DISPERSION FLOOR IS UPSIDE DOWN HERE, AND THAT MATTERS
##
## .countimp_draw_2l() floors the drawn dispersion at theta.min = 0.25. For
## nbinom2 that floor is a guard against imputation explosions: as theta -> 0
## the tails grow without bound (B09a). For compois the reported dispersion is
## sigma = 1/nu, so a SMALL sigma is strong UNDERdispersion -- the tails get
## tighter, not heavier, and the floor would silently cap nu at 4, i.e. throw
## away exactly the regime this family exists for. That is not a corner case:
## measured on three underdispersed designs (analyse/k25_2lcmp_boden.csv),
## 12 items at p ~ 0.83 gives sigma = 0.172 and at p ~ 0.90 sigma = 0.089 --
## both below the floor, which would have capped nu at 4 instead of 5.8 and
## 11.3. Only the mildest design (10 items at p ~ 0.6, nu = 2.7) stays above
## it.
##
## The floor is therefore switched off here (theta.min = 0), and nu is instead
## held inside the range the core can represent, [.countimp_cmp_nu_min,
## .countimp_cmp_nu_max] = [1e-3, 1e3], which is the same bound the single-level
## draw applies.


## The engine. `draw` is "bayes" or "boot".
##
## Shaped like .countimp_2l_count(): one model part, no zero mechanism, the
## type codes 3-6 unused. The two differences are the shape handling above and
## the sampler at the end.
.countimp_2l_cmp <- function(y, ry, x, type, draw = c("bayes", "boot"),
                             intercept = TRUE, wy = NULL, EV = FALSE) {
  draw <- match.arg(draw)
  if (is.null(wy)) wy <- !ry

  X   <- data.frame(x[ry, , drop = FALSE])
  nam <- colnames(X)
  dec <- .countimp_decode_type(type, nam)
  dat <- data.frame(Y = y[ry], X)

  dat.voll <- dat
  if (identical(draw, "boot")) dat <- .countimp_boot_clusters(dat, dec$group)
  ## from the FULL data, not the resample: a left-out cluster still has data
  ## and keeps its conditional posterior (see .countimp_boot_rate())
  obs.levels <- .countimp_obs_levels(dat.voll, dec$group)

  newdata <- data.frame(x[wy, , drop = FALSE])
  colnames(newdata) <- nam

  form <- .countimp_2l_formula(dec, "count", response = "Y",
                               intercept = intercept)
  fit <- glmmTMB::glmmTMB(formula = form, data = dat,
                          family = .countimp_family("compois", "log"))

  ## The shape, measured against this fit's own data and fitted means. Done
  ## once per fit and BEFORE the draw perturbs anything, so the reading is
  ## taken at the point estimate where the likelihood comparison is sharpest.
  mu.fit <- as.numeric(stats::predict(fit, type = "response"))
  ## The measurement returns the reading (`rezi`) alongside the value, because
  ## the drawn dispersion below has to be converted with the SAME reading. A
  ## caller that re-derives it from nu alone is wrong precisely where it
  ## matters: an unreadable sigma, or a sigma numerically equal to 1.
  sh <- .countimp_cmp_nu_glmm(fit, dat$Y, mu.fit)

  dr <- if (identical(draw, "bayes")) {
    ## theta.min = 0: see the note at the top of the file. The joint draw
    ## returns the perturbed sigma; exp(betadisp) IS sigma (measured), so the
    ## perturbation is on the right scale and in the right direction.
    .countimp_draw_2l(fit, newdata, grp = dec$group, obs_levels = obs.levels,
                      theta.min = 0)
  } else {
    ## The cluster bootstrap carries the parameter uncertainty; the shape comes
    ## from the resampled fit as a point estimate, exactly as the NB bootstrap
    ## variant takes theta from its own fit.
    list(mu = .countimp_boot_rate(fit, form, dat.voll,
                                  .countimp_family("compois", "log"), newdata,
                                  grp = dec$group, obs_levels = obs.levels),
         theta = .countimp_theta(fit, pflicht = FALSE))
  }

  nu.star <- suppressWarnings(as.numeric(dr$theta %||% NA_real_))
  nu.star <- if (!is.finite(nu.star) || nu.star <= 0) sh$nu else
    if (sh$rezi) 1 / nu.star else nu.star
  nu.star <- min(max(nu.star, .countimp_cmp_nu_min), .countimp_cmp_nu_max)

  mu <- as.numeric(dr$mu)
  imp <- if (any(!is.finite(mu)) || max(mu) > .countimp_cmp_mu_max)
    rep(NA_real_, length(mu)) else .countimp_rcmp(mu, nu.star)
  imp <- as.numeric(imp)

  if (EV) imp <- .countimp_ev_handle(imp, y, x, wy)
  imp
}
