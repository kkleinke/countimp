## Shared core for the single-level count methods --------------------------
##
## Ten exported names (poisson, pois, quasipoisson, qpois, nb, each with a
## .boot variant) reduce to three distributions x two draw routes. The four
## short names were byte-identical copies of the long ones, maintained twice.
##
## The three distributions differ in exactly three respects, and nothing else:
##
##   family        fitted with              variance of beta      count draw
##   ------------- ------------------------ --------------------- ------------
##   poisson       glm.fit(poisson)         cov.unscaled          rpois
##   quasipoisson  glm.fit(quasipoisson)    phi * cov.unscaled    .countimp_rqpois
##   nb            MASS::glm.nb             cov.unscaled          MASS::rnegbin
##
## The middle column is where the old code was wrong for quasi-Poisson (B40).
## summary.glm() reports cov.unscaled (the inverse information at phi = 1) and
## cov.scaled = phi * cov.unscaled separately. The sampling variance of beta
## under a quasi-likelihood with estimated dispersion is the SCALED one; the
## old code drew from the unscaled matrix for all three families, which for
## quasi-Poisson understates parameter uncertainty by a factor sqrt(phi).
## Poisson and glm.nb both fix phi = 1, so for those two the matrices coincide
## and the old code was right by coincidence, not by construction.
##
## Foreign-package contact is confined to .countimp_1l_fit(): MASS::glm.nb for
## the negative binomial, and summary.glm() -- which stats exports and which
## works on the plain list glm.fit() returns even though that list is not a
## glm object. The wrapper isolates that assumption in one place and falls back
## to computing the unscaled covariance from the QR decomposition itself if a
## future stats release tightens the method.


## Draw beta from its sampling distribution, or return the ML estimate.
##
## `scale` is the dispersion multiplier: 1 for families with fixed dispersion,
## phi for quasi-likelihood. It multiplies the VARIANCE, hence sqrt(scale) on
## the Cholesky factor.
.countimp_1l_draw_beta <- function(beta, cov_unscaled, scale = 1, bayes = TRUE) {
  if (!isTRUE(bayes)) return(beta)
  if (!is.finite(scale) || scale <= 0)
    stop("countimp: the dispersion estimate is ", format(scale),
         ", which cannot be used to scale the coefficient covariance. ",
         "This usually means the count model did not converge.", call. = FALSE)
  ch <- try(chol(cov_unscaled * scale), silent = TRUE)
  if (inherits(ch, "try-error"))
    stop("countimp: the coefficient covariance matrix is not positive ",
         "definite, so no coefficients can be drawn. Collinear predictors ",
         "are the usual cause; drop one of them or use the .boot variant.",
         call. = FALSE)
  as.vector(beta + t(ch) %*% stats::rnorm(ncol(ch)))
}


## Fit one single-level count model and report what the draw needs.
##
## Returns beta, the unscaled coefficient covariance, the dispersion scale for
## the beta draw, and theta for the negative binomial. This is the only place
## in the single-level count family that touches MASS or summary.glm().
.countimp_1l_fit <- function(x, y, dist) {
  if (identical(dist, "negbin")) {
    dat <- data.frame(y = y, x[, -1L, drop = FALSE])
    nam <- c("y", paste0("V", seq_len(ncol(x) - 1L)))
    names(dat) <- nam
    form <- stats::as.formula(paste("y ~", paste(nam[-1L], collapse = " + ")))
    fit <- try(MASS::glm.nb(form, data = dat), silent = TRUE)
    if (inherits(fit, "try-error"))
      stop("countimp: MASS::glm.nb() could not fit the negative binomial ",
           "model (", sub("^Error[^:]*: *", "", as.character(fit)),
           "). Try method \"poisson\" or \"quasipoisson\" for this variable.",
           call. = FALSE)
    fs <- summary(fit)
    return(list(beta = stats::coef(fit), cov = fs$cov.unscaled,
                scale = 1, theta = .countimp_draw_theta(fit, fs[["theta"]])))
  }

  fam <- switch(dist,
    poisson      = stats::poisson(link = "log"),
    quasipoisson = stats::quasipoisson(link = "log"),
    stop("countimp: unknown count distribution \"", dist, "\".", call. = FALSE))
  fit <- stats::glm.fit(x, y, family = fam)
  fs  <- .countimp_glm_summary(fit, ncol(x))
  ## Poisson fixes phi = 1; quasi-Poisson estimates it and it scales Var(beta).
  scale <- if (identical(dist, "quasipoisson")) fs$dispersion else 1
  list(beta = fit$coefficients, cov = fs$cov.unscaled,
       scale = scale, theta = NA_real_,
       dispersion = fs$dispersion)
}


## summary.glm() on the bare list glm.fit() returns.
##
## stats exports summary.glm and it happens to work on that list, but the list
## is not a glm object (no $terms, class "list"), so this is an undocumented
## convenience. Keep the assumption in one place and compute the two quantities
## we actually need from the QR factor if it ever stops holding.
.countimp_glm_summary <- function(fit, p) {
  fs <- try(stats::summary.glm(fit), silent = TRUE)
  if (!inherits(fs, "try-error") &&
      is.numeric(fs$dispersion) && !is.null(fs$cov.unscaled))
    return(fs)
  ## fallback: (X'WX)^-1 from the QR of the weighted design, and Pearson phi
  cu <- chol2inv(fit$qr$qr[seq_len(p), seq_len(p), drop = FALSE])
  df <- max(1L, length(fit$residuals) - fit$rank)
  phi <- if (fit$family$family %in% c("poisson", "binomial")) 1
         else sum((fit$weights * fit$residuals^2)[fit$weights > 0]) / df
  .countimp_note_event("glm_summary_fallback")
  list(cov.unscaled = cu, dispersion = phi)
}


## The single-level count engine.
##
## dist:  "poisson" | "quasipoisson" | "negbin"
## bayes: TRUE  -> draw beta from its sampling distribution (proper imputation)
##        FALSE -> bootstrap the observed cases, then use the ML estimate
.countimp_1l_count <- function(y, ry, x, wy = NULL, dist, bayes,
                               EV = FALSE, ...) {
  if (is.null(wy)) wy <- !ry
  x <- cbind(1, as.matrix(x))

  ## One draw. Kept as a closure so .countimp_draw_retry() can repeat it: every
  ## call refits (bootstrap) or redraws beta (Bayes) from the ambient RNG state,
  ## so a repeat is a genuinely new draw. Returning fit and mu alongside the
  ## draws is what lets countimp_check() see the fit-level problems (B56).
  ein_zug <- function() {
    if (isTRUE(bayes)) {
      f <- .countimp_1l_fit(x[ry, , drop = FALSE], y[ry], dist)
      beta.star <- .countimp_1l_draw_beta(f$beta, f$cov, f$scale, bayes = TRUE)
    } else {
      obs <- which(ry)
      sel <- sample(obs, length(obs), replace = TRUE)
      f <- .countimp_1l_fit(x[sel, , drop = FALSE], y[sel], dist)
      beta.star <- f$beta
    }

    ## No stop() on non-finite mu any more: countimp_check() classifies it as
    ## mu_nonfinite and the retry loop repeats the draw. Only when every
    ## attempt fails does the run abort -- with a message naming the count.
    mu <- exp(drop(x[wy, , drop = FALSE] %*% beta.star))

    im <- if (any(!is.finite(mu))) rep(NA_real_, length(mu)) else switch(dist,
      poisson      = stats::rpois(length(mu), mu),
      quasipoisson = .countimp_rqpois(mu, f$dispersion),
      negbin       = MASS::rnegbin(n = length(mu), mu = mu, theta = f$theta))

    if (isTRUE(EV) && all(is.finite(im))) im <- .countimp_ev_screen(im, y, ry, x, wy)
    list(imp = im, fit = f, mu = mu)
  }

  .countimp_draw_retry(ein_zug, y_obs = y[ry], method = dist)
}


## Extreme-value screening, shared by all six variants.
##
## Flag outliers and NaN among the draws, blank them, and redraw exactly those
## positions from the fit on the retained values. Previously copied verbatim
## into each of the ten methods.
.countimp_ev_screen <- function(im, y, ry, x, wy) {
  outl <- getOutliers(im, rho = c(0.3, 0.3), FLim = c(0.15, 0.85))
  idx  <- c(outl$iLeft, outl$iRight, which(is.nan(im)))
  if (length(idx) == 0L) return(im)
  im[idx] <- NA
  y[wy] <- im
  wy.ev <- logical(length(y))
  wy.ev[which(wy)[idx]] <- TRUE
  im[idx] <- .countimp_ev_refill(y, !is.na(y), x, wy = wy.ev)
  im
}
