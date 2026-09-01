## bounds.R -- interval-truncated counts: structural lower and upper bounds
##
## Why this file exists. Questionnaire scales have a maximum. "Number of
## symptoms out of 8", "days absent in a 5-day week", "correct answers out of
## 10" cannot exceed their ceiling, and a count with a floor above zero (see
## ztrunc.R) can have a ceiling as well. Until now countimp imputed such
## variables from an unbounded Poisson or NB and produced values the scale does
## not contain.
##
## Two things have to change, and only one of them is obvious.
##
## The obvious one is the draw: imputations must fall inside [lo, hi]. The
## tempting fix is to draw unbounded and clamp -- pmin(pmax(y, lo), hi). That is
## wrong, and the direction of the error is not what one would guess. Measured
## on a Poisson with bounds [1, 7], variance of the clamped draws divided by
## variance of the correct truncated draws:
##
##   mu =  1   0.74      mass outside 0.368   clamping INFLATES the floor
##   mu =  2   1.00      mass outside 0.136   about right by coincidence
##   mu =  4   1.16      mass outside 0.069   piles mass on both edges
##   mu =  8   0.73      mass outside 0.547
##   mu = 12   0.14      mass outside 0.911
##   mu = 20   0.00      mass outside 0.999   every draw becomes hi
##
## So clamping does not simply understate variance -- it overstates it in the
## interior of the range and destroys it once mu leaves the interval. There is no
## correction factor to apply; the draw has to come from the truncated
## distribution. (Reproduced by analyse/k21_bounds_ziehung.csv.)
##
## The less obvious one is the fit. If mass of the fitted distribution lies
## beyond the bound, an unbounded fit reads the piled-up ceiling observations as
## evidence for a flatter relationship. Measured on 3000 cases, true slope 0.45:
##
##   ceiling hi   mass beyond hi   unbounded slope   truncated slope
##          20          0.0000       0.420             0.420
##          10          0.0016       0.437             0.446
##           6          0.0261       0.379             0.435
##           4          0.1075       0.314             0.456
##           3          0.2005       0.241             0.410
##
## At hi = 4 the unbounded fit loses 30% of the slope. Where the bound does not
## bind (hi = 20) both agree, so the truncated fit costs nothing there. This is
## why bounds change the likelihood and not just the draw.
##
## Method. With a finite ceiling the support is a finite grid, so likelihood,
## gradient and draw all use the same summation over lo:hi -- no cdf
## differencing (which cancels catastrophically when mu is far outside the
## interval; ppois(7, 800) underflows to exactly 0 and the naive inverse-cdf
## draw then returns lo, the opposite of where the conditional mass is) and no
## rejection sampling (which would reject essentially every draw in that same
## regime). Everything runs on the log scale with a log-sum-exp shift, so a
## conditional distribution concentrated at one end is handled exactly:
## mu = 1e5 with hi = 7 returns 7 with probability 1, as it must.
##
## With hi = Inf this reduces to the pure lower truncation of ztrunc.R, which is
## the better route because it needs no grid; .countimp_bd_fit() therefore
## delegates to .countimp_zt_fit() when lo == 1 and hi == Inf.


## Normalise and validate a bounds specification.
##
## Accepted: c(lo, hi), a length-2 numeric with lo <= hi. Inf is allowed as the
## ceiling, -Inf is not a meaningful floor for a count and becomes 0.
## Non-integer bounds are an error rather than a rounding, because rounding
## silently changes the model the user asked for.
.countimp_bd_parse <- function(bounds, y = NULL, ry = NULL) {
  if (is.null(bounds)) return(NULL)
  if (!is.numeric(bounds) || length(bounds) != 2L)
    stop("'bounds' must be a numeric vector of length 2, c(lower, upper).",
         call. = FALSE)
  lo <- bounds[1L]; hi <- bounds[2L]
  if (is.na(lo) || is.na(hi))
    stop("'bounds' must not contain NA.", call. = FALSE)
  if (is.infinite(lo) && lo < 0) lo <- 0
  if (lo < 0)
    stop("the lower bound of a count cannot be negative (got ", lo, ").",
         call. = FALSE)
  if (lo >= hi)
    stop("'bounds' requires lower < upper (got ", lo, " and ", hi, ").",
         call. = FALSE)
  for (b in c(lo, hi))
    if (is.finite(b) && abs(b - round(b)) > 1e-8)
      stop("'bounds' must be whole numbers; got ", b,
           ". Rounding is not applied because it would change the model.",
           call. = FALSE)
  lo <- as.numeric(round(lo)); hi <- if (is.finite(hi)) as.numeric(round(hi)) else Inf

  ## Observed values outside the declared bounds mean the bounds are wrong, in
  ## the same way that an observed zero contradicts a zero-truncated model
  ## (B73). Naming the offending values is what makes this fixable.
  if (!is.null(y)) {
    yo <- if (is.null(ry)) y[!is.na(y)] else y[ry]
    bad <- yo[yo < lo | yo > hi]
    if (length(bad) > 0L)
      stop(length(bad), " observed value(s) lie outside bounds [", lo, ", ",
           format(hi), "]: ", paste(utils::head(sort(unique(bad)), 5L),
                                    collapse = ", "),
           if (length(unique(bad)) > 5L) ", ..." else "",
           ". Either the bounds are wrong or the variable is not bounded.",
           call. = FALSE)
  }
  list(lo = lo, hi = hi)
}


## Log KERNEL of the pmf on the grid lo:hi, one row per case.
##
## Returns an n x length(lo:hi) matrix, up to a row-wise additive constant. That
## constant is deliberately dropped: everything downstream uses only ratios
## within a row (the truncated likelihood is f(y)/sum_k f(k), the weights are
## f(k)/sum_j f(j)), and the dropped factor cancels exactly in both.
##
## Dropping it is what makes the extreme range work. Called with mu = 1e304 and
## hi = 5, stats::dpois() returns the same value -1.01e304 for every k, because
## the k-dependent part k*log(mu) - lgamma(k+1) is 24 orders of magnitude below
## the ulp of the common -mu term. The weights then come out uniform, so the
## truncated mean is the midpoint of the interval instead of hi, and the score
## points the wrong way -- which is how an optimiser ended up reporting a
## gradient of 1e186 on a flat objective (B78). In kernel form the same case is
## exact: k*700 - lgamma(k+1) is representable, and the weights concentrate at
## hi as they must.
##
## Parameterised by eta = log(mu) rather than mu, so no logarithm of a possibly
## underflowed mu is ever taken.
.countimp_bd_lgrid <- function(eta, lo, hi, dist, theta = NULL) {
  k <- lo:hi
  lgk <- lgamma(k + 1)
  if (dist == "negbin") {
    ## log(mu/(theta+mu)) written as eta - log(theta+mu): both terms are
    ## representable for every eta in [-745, 700].
    lp <- eta - log(theta + exp(eta))
    outer(lp, k, function(a, kk) a * kk) +
      matrix(lgamma(k + theta) - lgk, nrow = length(eta),
             ncol = length(k), byrow = TRUE)
  } else {
    outer(eta, k, function(e, kk) e * kk) -
      matrix(lgk, nrow = length(eta), ncol = length(k), byrow = TRUE)
  }
}


## Row-wise log-sum-exp. The shift is what keeps a concentrated conditional
## distribution representable: without it, exp(-800) underflows to 0 and the
## normalising constant becomes 0/0.
.countimp_bd_lse <- function(lg) {
  mx <- apply(lg, 1L, max)
  mx + log(rowSums(exp(lg - mx)))
}


## Negative log-likelihood of the interval-truncated count model.
##
## par = beta for Poisson, c(beta, log theta) for NB.
##
## Written as a log-softmax over the grid rather than as
## log f(y) - log sum_k f(k) computed from two separate quantities. Because the
## observed y is inside [lo, hi] (guaranteed by .countimp_bd_parse), it *is* one
## of the grid cells, so the likelihood is simply that cell's share of the row.
## The two-quantity form has both terms around -1e304 when mu is huge and their
## difference is then formed by cancellation; the log-softmax never leaves a
## moderate range.
.countimp_bd_nll <- function(par, x, y, lo, hi, dist) {
  p <- ncol(x)
  beta <- par[seq_len(p)]
  theta <- if (dist == "negbin") exp(par[p + 1L]) else NULL
  eta <- pmin(pmax(drop(x %*% beta), -700), 700)
  lg <- .countimp_bd_lgrid(eta, lo, hi, dist, theta)
  idx <- cbind(seq_along(y), as.integer(y - lo + 1L))
  val <- -sum(lg[idx] - .countimp_bd_lse(lg))
  if (!is.finite(val)) 1e10 else val
}


## Analytic gradient of the interval-truncated negative log-likelihood.
##
## The score of a truncated count model has a form that is both exact and
## numerically harmless, and it is worth writing out because the obvious
## implementation is neither.
##
## Obvious route: d/dbeta = -x' ((dmu_obs - dmu_norm) * mu), where dmu_* are
## derivatives of the log density and of the log normaliser with respect to mu.
## Correct, and it overflows: when eta hits its cap, mu is 1e304, the bracket is
## a difference of two nearly equal small numbers, and their product with mu came
## out as 1e186 -- an optimiser thrown off a flat objective by a gradient that
## disagreed with it (B78).
##
## Exact route: the mu factor cancels analytically before it is ever formed.
## With w_k the truncated probabilities and Ew = sum_k w_k * k,
##
##   Poisson:  (dmu_obs - dmu_norm) * mu = (y - mu)   - (Ew - mu) = y - Ew
##   NB2:      (dmu_obs - dmu_norm) * mu = (y - Ew) * theta / (theta + mu)
##
## i.e. the score is the residual against the *truncated* mean, which is the
## textbook exponential-family form. Nothing large is ever multiplied by
## anything large, and the Poisson case needs no division at all.
##
## For theta the terms are all logarithmic in magnitude, so the direct form is
## used there. The sign on the truncation term is the part that is easy to get
## wrong (B69), so test-B75 checks each component against finite differences
## separately.
.countimp_bd_gr <- function(par, x, y, lo, hi, dist) {
  p <- ncol(x)
  beta <- par[seq_len(p)]
  theta <- if (dist == "negbin") exp(par[p + 1L]) else NULL
  eta <- pmin(pmax(drop(x %*% beta), -700), 700)
  mu <- exp(eta)
  k <- lo:hi

  lg <- .countimp_bd_lgrid(eta, lo, hi, dist, theta)
  w <- exp(lg - .countimp_bd_lse(lg))              # truncated weights, rows sum 1
  Ew <- drop(w %*% k)                              # truncated mean, per case

  resid <- y - Ew
  if (dist == "negbin") resid <- resid * theta / (theta + mu)
  g_beta <- -drop(t(x) %*% resid)

  if (dist != "negbin") return(g_beta)

  ## d/dtheta of log f, then of the log normaliser; reported for log theta
  dth_cell <- outer(mu, k, function(m, kk)
    digamma(kk + theta) - digamma(theta) + log(theta / (theta + m)) + 1 -
      (kk + theta) / (theta + m))
  dth_obs <- digamma(y + theta) - digamma(theta) + log(theta / (theta + mu)) +
    1 - (y + theta) / (theta + mu)
  g_theta <- -sum((dth_obs - rowSums(w * dth_cell)) * theta)
  c(g_beta, g_theta)
}


## Fit the interval-truncated count model by maximum likelihood.
##
## Delegates to .countimp_zt_fit() for the pure lower-truncation case, which
## needs no grid. Starting values come from the unbounded glm, which is the
## right starting point precisely because it is biased in a known direction.
.countimp_bd_fit <- function(x, y, dist, lo, hi, max_grid = 5000L) {
  if (lo == 1 && !is.finite(hi))
    return(.countimp_zt_fit(x, y, if (dist == "negbin") "ztnb" else "ztp"))
  if (!is.finite(hi))
    stop("an infinite upper bound is only supported with lower bound 0 or 1; ",
         "use method \"", if (dist == "negbin") "ztnb" else "ztp",
         "\" for a lower bound of 1.", call. = FALSE)
  if (hi - lo + 1L > max_grid)
    stop("the bounds span ", format(hi - lo + 1),
         " values, above the grid limit of ", max_grid,
         ". A bound that wide almost never binds; impute without 'bounds'.",
         call. = FALSE)

  fam <- if (dist == "negbin") stats::poisson() else
    if (dist == "quasipoisson") stats::quasipoisson() else stats::poisson()
  st <- suppressWarnings(stats::glm.fit(x, y, family = fam))$coefficients
  st[!is.finite(st)] <- 0
  if (dist == "negbin") st <- c(st, log(.countimp_mom_theta(y)))

  op <- stats::optim(st, .countimp_bd_nll, .countimp_bd_gr,
                     x = x, y = y, lo = lo, hi = hi, dist = dist,
                     method = "BFGS", hessian = TRUE,
                     control = list(reltol = 1e-10, maxit = 500L))
  p <- ncol(x)
  theta <- if (dist == "negbin") exp(op$par[p + 1L]) else NULL
  ## Field names follow the .countimp_1l_fit()/.countimp_zt_fit() contract
  ## exactly -- beta, cov, scale, theta, ll, conv, nobs, npar -- and the object
  ## carries NO class, like the other two. Both are deliberate: logLik.countimp_zt_fit()
  ## reads $ll/$npar/$nobs, and countimp_check() recognises an internal fit by
  ## its fields. Inventing $loglik/$converged here, or attaching a class, broke
  ## both silently (B77).
  list(beta = op$par[seq_len(p)],
       cov = .countimp_zt_cov(op$hessian, length(op$par))[seq_len(p),
                                                         seq_len(p),
                                                         drop = FALSE],
       scale = 1, theta = if (is.null(theta)) NA_real_ else theta,
       ll = -op$value, conv = op$convergence, nobs = length(y),
       npar = length(op$par), dist = dist, lo = lo, hi = hi)
}


## Draw from the interval-truncated distribution.
##
## Exact categorical draw on the log-scale grid. Not inverse cdf: the cdf
## difference underflows to zero when mu is far outside [lo, hi] and the draw
## then collapses onto the wrong bound.
.countimp_rint <- function(mu, lo, hi, dist = "poisson", theta = NULL) {
  n <- length(mu)
  if (n == 0L) return(numeric(0))
  if (!is.finite(hi))
    return(if (dist == "negbin")
      .countimp_rktnb(n, size = theta, k = lo - 1, mu = mu)
      else .countimp_rktp(n, k = lo - 1, mu = mu))
  eta <- pmin(pmax(log(pmax(mu, .Machine$double.xmin)), -700), 700)
  lg <- .countimp_bd_lgrid(eta, lo, hi, dist, theta)
  w <- exp(lg - .countimp_bd_lse(lg))
  k <- lo:hi
  ## One uniform per case against the row cumulative -- vectorised, so cost is
  ## the grid width and not a per-case R loop.
  cw <- t(apply(w, 1L, cumsum))
  if (length(k) == 1L) return(rep(k, n))
  k[1L + rowSums(cw < stats::runif(n))]
}


## Single-level interval-truncated imputation, Poisson and NB.
##
## Same shape as .countimp_1l_ztrunc(): bounds are validated once, before the
## retry loop, because a value outside the bounds is a data problem and no
## number of redraws will fix it (B73).
.countimp_1l_bounded <- function(y, ry, x, wy = NULL, dist, bayes,
                                 bounds, EV = FALSE, ...) {
  bd <- .countimp_bd_parse(bounds, y = y, ry = ry)
  if (is.null(bd)) stop("'bounds' is required for this method.", call. = FALSE)
  if (is.null(wy)) wy <- !ry
  x <- cbind(1, as.matrix(x))

  ein_zug <- function() {
    if (isTRUE(bayes)) {
      f <- .countimp_bd_fit(x[ry, , drop = FALSE], y[ry], dist, bd$lo, bd$hi)
      beta.star <- .countimp_1l_draw_beta(f$beta, f$cov, f$scale, bayes = TRUE)
    } else {
      obs <- which(ry)
      sel <- sample(obs, length(obs), replace = TRUE)
      f <- .countimp_bd_fit(x[sel, , drop = FALSE], y[sel], dist, bd$lo, bd$hi)
      beta.star <- f$beta
    }
    mu <- exp(drop(x[wy, , drop = FALSE] %*% beta.star))
    th <- if (is.null(f$theta) || is.na(f$theta)) NULL else f$theta
    im <- if (any(!is.finite(mu))) rep(NA_real_, length(mu)) else
      .countimp_rint(mu, bd$lo, bd$hi, dist, th)
    if (isTRUE(EV) && all(is.finite(im))) im <- .countimp_ev_screen(im, y, ry, x, wy)
    ## EV screening refills from an unbounded fit, so re-impose the support.
    ## This is a redraw of the offending positions, not a clamp.
    if (isTRUE(EV) && all(is.finite(im))) {
      out <- which(im < bd$lo | im > bd$hi)
      if (length(out) > 0L)
        im[out] <- .countimp_rint(mu[out], bd$lo, bd$hi, dist, th)
    }
    list(imp = im, fit = f, mu = mu)
  }

  .countimp_draw_retry(ein_zug, y_obs = y[ry], method = dist)
}
