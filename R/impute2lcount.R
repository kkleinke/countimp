## impute2lcount.R -- one engine for the plain two-level count methods
##
## Up to 3.0.0 the 12 exported names in mice.impute.2l.poisson.R were 12 copies
## of one ~68-line body, differing in the count family (poisson / nbinom2), the
## intercept flag and the draw. Four of the twelve were never distinct at all:
## 2l.nb2* is documented as "identical to 2l.nb*; kept for backward
## compatibility", and a body diff confirmed it -- the only differences were a
## variable name, a brace placement and a comment. Those four are now real
## assignments, so they cannot drift apart from their originals.
##
## These methods have ONE model part, which is what separates this engine from
## impute2lzi.R (joint fit, joint draw) and impute2lhurdle.R (two factorised
## fits). There is no zero mechanism: a zero is an ordinary count that the
## Poisson or NB law produces on its own.

## The engine. `family` is "poisson" or "nbinom2"; `draw` is "bayes" or "boot".
.countimp_2l_count <- function(y, ry, x, type,
                               family = c("poisson", "nbinom2"),
                               draw = c("bayes", "boot"),
                               intercept = TRUE, wy = NULL, EV = FALSE) {
  family <- match.arg(family)
  draw   <- match.arg(draw)
  if (is.null(wy)) wy <- !ry

  X   <- data.frame(x[ry, , drop = FALSE])
  nam <- colnames(X)
  dec <- .countimp_decode_type(type, nam)
  dat <- data.frame(Y = y[ry], X)

  ## Resample whole clusters -- but read the observed levels from the FULL
  ## data, not from the resample. A left-out cluster is not an unseen one: its
  ## rows are in `dat.voll`, so it must keep its conditional posterior rather
  ## than fall back to the prior (see .countimp_boot_rate()).
  dat.voll <- dat
  if (identical(draw, "boot")) dat <- .countimp_boot_clusters(dat, dec$group)
  obs.levels <- .countimp_obs_levels(dat.voll, dec$group)

  newdata <- data.frame(x[wy, , drop = FALSE])
  colnames(newdata) <- nam

  ## Single model part: the type codes 3-6 that split predictors between a
  ## count and a zero model have no meaning here, so "count" collects codes
  ## 1-4 and the zero-only codes are simply unused.
  form <- .countimp_2l_formula(dec, "count", response = "Y",
                               intercept = intercept)
  fit <- glmmTMB::glmmTMB(formula = form, data = dat,
                          family = .countimp_family(family, "log"))
  if (identical(family, "nbinom2"))
    fit <- .countimp_stabilize(fit, form, dat)$fit

  dr <- if (identical(draw, "bayes")) {
    ## Joint draw of (beta, u, log theta) from the Laplace posterior. Drawing
    ## beta alone and predicting from the mutated fit lets predict() re-optimise
    ## u and absorb most of the draw -- see .countimp_draw_2l().
    .countimp_draw_2l(fit, newdata, grp = dec$group, obs_levels = obs.levels)
  } else {
    ## The cluster bootstrap carries the uncertainty in the regression
    ## parameters, but a resample yields one point estimate of theta per fit --
    ## .countimp_draw_theta() adds theta's own draw and enforces the lower
    ## bound that guards against the NB imputation explosions of B09a.
    ## .countimp_boot_rate() inverts the fit's own link instead of hard-coding
    ## exp(), and draws u_j for every cluster that has data -- including the
    ## ones the resample left out.
    list(mu = .countimp_boot_rate(fit, form, dat.voll,
                                  .countimp_family(family, "log"), newdata,
                                  grp = dec$group, obs_levels = obs.levels),
         theta = if (identical(family, "nbinom2")) .countimp_draw_theta(fit) else NA_real_)
  }

  imp <- if (identical(family, "poisson")) {
    stats::rpois(length(dr$mu), lambda = dr$mu)
  } else {
    MASS::rnegbin(length(dr$mu), theta = dr$theta, mu = dr$mu)
  }
  imp <- as.numeric(imp)

  if (EV) imp <- .countimp_ev_handle(imp, y, x, wy)
  imp
}
