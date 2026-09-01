## ztrunc.R -- zero-truncated Poisson and negative binomial, single level
##
## Counts that structurally cannot be zero: hospital length of stay measured in
## days for admitted patients, number of purchases among buyers, litter size
## among mothers who gave birth. Fitting an ordinary Poisson to those data
## understates mu, because the model spends probability mass on a zero that the
## sampling scheme excluded. B04 (the underdispersion note) and the manual's
## hurdle chapter both describe the zero-truncated count model, but until now
## countimp offered it only as the count half of a hurdle -- never on its own.
##
## Why fit it here instead of calling glmmTMB::truncated_poisson()?
##
## The first reason is dependency isolation, the project's stated goal: a
## single-level zero-truncated imputation must not require glmmTMB, TMB and a
## working C++ toolchain. The multilevel hurdle methods still use glmmTMB,
## because random effects genuinely need it.
##
## The second is speed, but only because of the analytic gradients below. Both
## sides measured on 1200 cases with the compiled TMB template already loaded
## (a cold first call costs glmmTMB about 0.9 s of template loading, which
## would flatter this code by a factor of ~45 and mean nothing). Five data sets
## per distribution, best of 15 timings each -- a single timing is worthless
## here: one ztp run came out at 0.7x, i.e. apparently slower, purely as noise.
##
##   ztp   median 3.8x   (range 3.0 - 4.3)
##   ztnb  median 5.1x   (range 4.0 - 6.0)
##
## Reproduced by analyse/k20_ztrunc_genauigkeit.csv.
##
## With numeric gradients the advantage disappears -- BFGS then needs p+1 extra
## likelihood evaluations per step, and for ztnb it wandered off in theta on
## samples where the truncated data carry little information about it. The
## gradients are what makes the pure-R fit competitive, not R being fast.
##
## Accuracy against glmmTMB on the same data is at optimiser tolerance: beta to
## 3e-5, standard errors to 2e-6, theta to 3e-5, log-likelihood to 1e-6
## (test-B69).


## log(1 - exp(-a)) for a >= 0.
##
## The naive expression loses all precision twice over: for small a,
## 1 - exp(-a) rounds to zero (exp(-1e-8) is 1 in double), and for large a,
## exp(-a) underflows. Maechler's (2012) split at log 2 keeps both ends exact:
## expm1 below, log1p above.
.countimp_log1mexp <- function(a) {
  ifelse(a <= log(2), log(-expm1(-a)), log1p(-exp(-a)))
}


## Negative log-likelihood, zero-truncated Poisson.
##
## ll_i = dpois(y_i; mu_i, log) - log(1 - exp(-mu_i))
.countimp_ztp_nll <- function(beta, x, y) {
  mu <- exp(drop(x %*% beta))
  if (any(!is.finite(mu)) || any(mu <= 0)) return(1e10)
  -sum(stats::dpois(y, mu, log = TRUE) - .countimp_log1mexp(mu))
}


## Gradient, zero-truncated Poisson.
##
## d ll_i / d beta_j = (y_i - mu_i - mu_i / (exp(mu_i) - 1)) x_ij
## expm1(mu) rather than exp(mu) - 1: for mu below about 1e-8 the latter is
## catastrophically inaccurate, and mu that small occurs whenever a predictor
## takes an extreme value during a Bayesian draw.
.countimp_ztp_gr <- function(beta, x, y) {
  mu <- exp(drop(x %*% beta))
  if (any(!is.finite(mu)) || any(mu <= 0)) return(rep(0, length(beta)))
  -drop(crossprod(x, y - mu - mu / expm1(mu)))
}


## Negative log-likelihood, zero-truncated NB2.
##
## P(Y = 0) = (theta / (theta + mu))^theta = exp(a) with a <= 0, so
## log(1 - P(Y = 0)) = .countimp_log1mexp(-a). theta is optimised on the log
## scale, which keeps it positive without a constraint.
.countimp_ztnb_nll <- function(par, x, y) {
  p  <- ncol(x)
  mu <- exp(drop(x %*% par[seq_len(p)]))
  th <- exp(par[p + 1L])
  if (any(!is.finite(mu)) || any(mu <= 0) || !is.finite(th) || th <= 0)
    return(1e10)
  a <- th * (log(th) - log(th + mu))
  -sum(stats::dnbinom(y, size = th, mu = mu, log = TRUE) -
       .countimp_log1mexp(-a))
}


## Gradient, zero-truncated NB2, with respect to (beta, log theta).
##
## The sign on the truncation term is the part that is easy to get wrong:
## ll_i carries -log(1 - P0), so d/dth of that is +P0 * (da/dth) / (1 - P0).
## Getting it backwards leaves the beta components exact and the theta
## component off by a factor of -90 -- the finite-difference check in test-B69
## is there because that error is invisible in the fitted coefficients.
.countimp_ztnb_gr <- function(par, x, y) {
  p    <- ncol(x)
  mu   <- exp(drop(x %*% par[seq_len(p)]))
  th   <- exp(par[p + 1L])
  if (any(!is.finite(mu)) || any(mu <= 0) || !is.finite(th) || th <= 0)
    return(rep(0, length(par)))
  s    <- th + mu
  P0   <- exp(th * (log(th) - log(s)))
  om   <- 1 - P0
  ## d ll_i / d mu, then chain rule d mu / d beta = mu * x
  dmu  <- y / mu - (y + th) / s - th * P0 / (s * om)
  gb   <- crossprod(x, dmu * mu)
  ## d ll_i / d theta, then chain rule d theta / d log theta = theta
  da   <- log(th) - log(s) + 1 - th / s
  dth  <- digamma(y + th) - digamma(th) + log(th) - log(s) + 1 - th / s -
          y / s + P0 * da / om
  c(-drop(gb), -sum(dth) * th)
}


## Fit a zero-truncated count model by maximum likelihood.
##
## Starting values come from an ordinary Poisson GLM on the truncated data.
## That overstates mu -- the fit has to explain the missing zeros somehow --
## but it is in the basin of attraction, and it costs one IRLS run.
##
## Returns beta, its covariance from the numerically differentiated Hessian,
## theta (NA for ztp), the log-likelihood and a convergence flag. Same shape as
## .countimp_1l_fit() so the draw engine treats both alike.
## Refuse zero-truncated methods on data that contain zeros.
##
## Called from the engine before the retry loop, not from the fit: a zero is a
## statement about the data, and repeating the draw cannot make it go away. The
## message names the methods that do model zeros, because "wrong method" is the
## actual finding and the user needs the alternative, not a diagnosis.
.countimp_zt_check_zeros <- function(y, dist) {
  n0 <- sum(y < 1, na.rm = TRUE)
  if (n0 == 0L) return(invisible(TRUE))
  stop("countimp: method \"", dist, "\" is for zero-truncated counts, but ",
       n0, " observed value(s) are zero or negative.\n",
       "  Zeros cannot occur under this model. Use \"hp\"/\"hnb\" (hurdle) or ",
       "\"zip\"/\"zinb\" (zero-inflated)\n  when zeros are present, or ",
       "\"poisson\"/\"nb\" when they are ordinary counts.", call. = FALSE)
}


.countimp_zt_fit <- function(x, y, dist) {
  start <- stats::glm.fit(x, y, family = stats::poisson())$coefficients
  if (any(!is.finite(start))) start <- c(log(mean(y)), rep(0, ncol(x) - 1L))

  if (identical(dist, "ztp")) {
    o <- try(stats::optim(start, .countimp_ztp_nll, .countimp_ztp_gr,
                          x = x, y = y, method = "BFGS", hessian = TRUE),
             silent = TRUE)
    if (inherits(o, "try-error") || o$value >= 1e10)
      stop("countimp: could not fit the zero-truncated Poisson model.",
           call. = FALSE)
    cov <- .countimp_zt_cov(o$hessian, ncol(x))
    return(list(beta = o$par, cov = cov, scale = 1, theta = NA_real_,
                ll = -o$value, conv = o$convergence, nobs = length(y),
                npar = ncol(x)))
  }

  o <- try(stats::optim(c(start, 0), .countimp_ztnb_nll, .countimp_ztnb_gr,
                        x = x, y = y, method = "BFGS", hessian = TRUE),
           silent = TRUE)
  if (inherits(o, "try-error") || o$value >= 1e10)
    stop("countimp: could not fit the zero-truncated negative binomial ",
         "model. Try method \"ztp\" for this variable.", call. = FALSE)
  p <- ncol(x)
  list(beta = o$par[seq_len(p)], cov = .countimp_zt_cov(o$hessian, p),
       scale = 1, theta = exp(o$par[p + 1L]), ll = -o$value,
       conv = o$convergence, nobs = length(y), npar = p + 1L)
}


## Coefficient covariance from the Hessian of the negative log-likelihood.
##
## solve() fails on a singular Hessian, which happens with collinear predictors
## or a theta that ran to its boundary. Fall back to the generalised inverse so
## the draw can proceed with whatever the data identified; the retry loop in
## .countimp_draw_retry() catches a genuinely unusable draw downstream.
.countimp_zt_cov <- function(H, p) {
  cov <- try(solve(H), silent = TRUE)
  if (inherits(cov, "try-error")) {
    ev  <- eigen(H, symmetric = TRUE)
    pos <- ev$values > max(ev$values) * 1e-10
    cov <- ev$vectors[, pos, drop = FALSE] %*%
           diag(1 / ev$values[pos], nrow = sum(pos)) %*%
           t(ev$vectors[, pos, drop = FALSE])
    .countimp_note_event("zt_hessian_singular")
  }
  cov[seq_len(p), seq_len(p), drop = FALSE]
}


## The single-level zero-truncated draw engine.
##
## Same structure as .countimp_1l_count(): one closure per draw, handed to
## .countimp_draw_retry() so a hard failure is repeated rather than fatal.
##
## dist:  "ztp" | "ztnb"
## bayes: TRUE  -> draw beta from N(betahat, cov) before drawing counts
##        FALSE -> bootstrap the observed cases, use the resulting estimate
.countimp_1l_ztrunc <- function(y, ry, x, wy = NULL, dist, bayes,
                                EV = FALSE, ...) {
  if (is.null(wy)) wy <- !ry
  x <- cbind(1, as.matrix(x))

  ## Observed zeros are a data problem, not a failed draw, so the check belongs
  ## here and not inside .countimp_zt_fit(): raised from inside the closure,
  ## .countimp_draw_retry() repeated it three times and then reported it as a
  ## convergence failure ("extreme predictor values, a separated model ...")
  ## with the actual cause buried in a parenthesis. Wrong diagnosis for a user
  ## who simply picked the wrong method (B73).
  .countimp_zt_check_zeros(y[ry], dist)

  ein_zug <- function() {
    if (isTRUE(bayes)) {
      f <- .countimp_zt_fit(x[ry, , drop = FALSE], y[ry], dist)
      beta.star <- .countimp_1l_draw_beta(f$beta, f$cov, f$scale, bayes = TRUE)
    } else {
      obs <- which(ry)
      sel <- sample(obs, length(obs), replace = TRUE)
      f <- .countimp_zt_fit(x[sel, , drop = FALSE], y[sel], dist)
      beta.star <- f$beta
    }

    mu <- exp(drop(x[wy, , drop = FALSE] %*% beta.star))

    ## The draw itself is zero-truncated, so imputations are >= 1 by
    ## construction -- no rejection loop, no post-hoc filtering. This is the
    ## whole point of the method: an ordinary Poisson draw would put mass on
    ## zero for a variable that cannot be zero.
    im <- if (any(!is.finite(mu))) rep(NA_real_, length(mu)) else switch(dist,
      ztp  = .countimp_rktp(length(mu), k = 0, mu = mu),
      ztnb = .countimp_rktnb(length(mu), size = f$theta, k = 0, mu = mu))

    if (isTRUE(EV) && all(is.finite(im))) im <- .countimp_ev_screen(im, y, ry, x, wy)
    list(imp = im, fit = f, mu = mu)
  }

  .countimp_draw_retry(ein_zug, y_obs = y[ry], method = dist)
}


## logLik method for the zero-truncated fit.
##
## Registered so that stats::AIC(), .countimp_npar() and .countimp_fit_usable()
## work on this object without a special case each. Adding four `inherits()`
## branches in R/fitdiag.R would have been the direct route and the wrong one:
## the generic exists precisely so that a new model class costs one method, not
## one branch per consumer.
#' @exportS3Method stats::logLik countimp_zt_fit
logLik.countimp_zt_fit <- function(object, ...) {
  structure(object$ll, df = object$npar, nobs = object$nobs,
            class = "logLik")
}
