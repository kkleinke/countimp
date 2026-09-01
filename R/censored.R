## censored.R -- right-censored counts: the observed value is a lower bound
##
## Why this file exists. Top-coded counts are everywhere in survey data: "10 or
## more visits", "5+ children", an item scale that stops recording above its
## last category, a diary that was closed after 30 events. The recorded number
## is then not the count -- it is the statement "the count was at least this
## large". Until now countimp had no way to say that, and both halves of the
## job went wrong:
##
##   1. the FIT read the piled-up ceiling values as if they were the actual
##      counts, which pulls mu down and flattens the slope, and
##   2. the DRAW produced imputations on a scale the data do not have.
##
## This is a different statement from mice.impute.bp()'s bounds. A bound says
## the count CANNOT exceed hi (a five-day week has no sixth day); censoring says
## the count can, but was not recorded. The likelihoods differ accordingly: the
## bounded model renormalises over [lo, hi], the censored model keeps the whole
## support and gives a censored case the tail probability P(Y >= c). Handing
## censored data to "bp" understates the mean, handing bounded data to "cp"
## invents values the scale does not contain.
##
## Two things the file does NOT do, deliberately:
##   * left censoring ("2 or fewer") is not implemented. It is the mirror image
##     and would be a small addition, but nothing in the manual or the register
##     asks for it, and an untested mirror is worse than an honest gap.
##   * interval censoring is not implemented for the same reason.
##
## Method. The likelihood needs P(Y >= c), which stats gives exactly and on the
## log scale (ppois/pnbinom with lower.tail = FALSE, log.p = TRUE), so no grid
## is needed for the value. The gradient with respect to the linear predictor is
## a hazard-type ratio f(c-1)/S(c), formed as exp(log f - log S) so that a
## censored case whose tail probability has underflowed to zero does not become
## Inf/Inf. Only the NB dispersion needs a grid, and only a FINITE one: because
##
##     d/dtheta log S(c) = -[sum_{j < c} f(j) dlog f(j)/dtheta] / S(c)
##
## the sum runs over the observable part 0..c-1, not over the censored tail.
## That is what makes the analytic theta gradient affordable at all.


## Normalise and validate a censoring specification.
##
## Accepted:
##   censor = 10        one ceiling for every case (top coding)
##   censor = vec       per case, length n; Inf marks a case that is not
##                      censored at all
##
## A case counts as censored when y >= c, which is the only rule the data
## support: top coding records "10 or more" as 10, so an observed 10 is
## indistinguishable from a censored one. Where an explicit indicator exists,
## the per-case form expresses it exactly -- set c_i = y_i for the censored
## cases and Inf for the rest.
##
## Observed values ABOVE the ceiling are an error, not a warning: they mean the
## ceiling is wrong, in the same way that an observed zero contradicts a
## zero-truncated model (B73).
.countimp_cs_parse <- function(censor, n, y = NULL, ry = NULL) {
  if (is.null(censor))
    stop("'censor' is required for this method: give the censoring limit, ",
         "e.g. censor = 10 for a variable top-coded at 10.", call. = FALSE)
  if (!is.numeric(censor))
    stop("'censor' must be numeric: a single limit, or one limit per case.",
         call. = FALSE)
  if (anyNA(censor))
    stop("'censor' must not contain NA. Use Inf for a case that is not ",
         "censored.", call. = FALSE)
  if (length(censor) != 1L && length(censor) != n)
    stop("'censor' has length ", length(censor), "; it must have length 1 ",
         "(one limit for every case) or ", n, " (one per case).", call. = FALSE)
  cc <- rep_len(as.numeric(censor), n)
  if (any(cc < 1))
    stop("a censoring limit below 1 leaves no observable counts (got ",
         min(cc), ").", call. = FALSE)
  fin <- is.finite(cc)
  if (any(abs(cc[fin] - round(cc[fin])) > 1e-8))
    stop("'censor' must be whole numbers; got ",
         format(cc[fin][which(abs(cc[fin] - round(cc[fin])) > 1e-8)[1L]]),
         ". Rounding is not applied because it would change the model.",
         call. = FALSE)
  cc[fin] <- round(cc[fin])

  if (!is.null(y)) {
    obs <- if (is.null(ry)) !is.na(y) else ry
    bad <- which(obs & y > cc)
    if (length(bad))
      stop(length(bad), " observed value(s) exceed their censoring limit, ",
           "e.g. y = ", y[bad[1L]], " with censor = ", format(cc[bad[1L]]),
           ".\n  A censored observation is recorded AT the limit, never above ",
           "it, so either the\n  limit is wrong or the variable is not ",
           "censored.", call. = FALSE)
  }
  cc
}


## Which observed cases are censored. One place, so the fit, the draw and the
## scale decision cannot disagree about it.
.countimp_cs_is_cen <- function(y, cc) is.finite(cc) & y >= cc


## Negative log-likelihood of the right-censored count model.
##
## par = beta for Poisson, c(beta, log theta) for NB.
##   uncensored:  log f(y_i)
##   censored:    log P(Y >= c_i) = log S(c_i - 1)
##
## eta is capped at +-700 exactly as in bounds.R: exp(710) is Inf and every
## quantity downstream then becomes NaN, which the optimiser cannot walk back
## out of.
.countimp_cs_nll <- function(par, x, y, cc, dist) {
  p <- ncol(x)
  eta <- pmin(pmax(drop(x %*% par[seq_len(p)]), -700), 700)
  mu <- exp(eta)
  th <- if (dist == "negbin") exp(par[p + 1L]) else NULL
  if (dist == "negbin" && (!is.finite(th) || th <= 0)) return(1e10)
  cen <- .countimp_cs_is_cen(y, cc)
  ll <- numeric(length(y))
  if (any(!cen))
    ll[!cen] <- if (dist == "negbin")
      stats::dnbinom(y[!cen], size = th, mu = mu[!cen], log = TRUE)
    else stats::dpois(y[!cen], mu[!cen], log = TRUE)
  if (any(cen))
    ll[cen] <- if (dist == "negbin")
      stats::pnbinom(cc[cen] - 1, size = th, mu = mu[cen],
                     lower.tail = FALSE, log.p = TRUE)
    else stats::ppois(cc[cen] - 1, mu[cen], lower.tail = FALSE, log.p = TRUE)
  val <- -sum(ll)
  if (is.finite(val)) val else 1e10
}


## Analytic gradient of the right-censored negative log-likelihood.
##
## Uncensored cases give the ordinary score. For a censored case the derivative
## of log S(c-1) with respect to the linear predictor is
##
##   Poisson:  mu f(c-1) / S(c-1)
##   NB2:      mu f(c-1) (theta + c - 1) / ((theta + mu) S(c-1))
##
## using d/dmu P(Y <= k) = -f(k) (Poisson) and -f(k)(theta+k)/(theta+mu) (NB2);
## the NB2 form follows from d I_p(theta, k+1)/dp and reduces to the Poisson one
## as theta -> Inf, which is the cheap check that the (theta+k)/(theta+mu)
## factor is the right way up.
##
## Every one of those is formed as exp(log f - log S + ...): S underflows to 0
## for a censored case far in the tail of the current fit, and the ratio is then
## 0/0 = NaN in the direct form while the log form stays finite (it is a hazard
## rate, which is well behaved even where both parts have underflowed).
##
## For theta the sum runs over 0..c-1 as described at the top of the file. The
## row mask is what allows per-case limits: the grid is rectangular at
## max(c) - 1 and every case ignores the cells beyond its own limit.
.countimp_cs_gr <- function(par, x, y, cc, dist) {
  p <- ncol(x)
  eta <- pmin(pmax(drop(x %*% par[seq_len(p)]), -700), 700)
  mu <- exp(eta)
  th <- if (dist == "negbin") exp(par[p + 1L]) else NULL
  if (dist == "negbin" && (!is.finite(th) || th <= 0)) return(rep(0, p + 1L))
  cen <- .countimp_cs_is_cen(y, cc)
  d_eta <- numeric(length(y))
  d_th <- numeric(length(y))

  if (any(!cen)) {
    m <- mu[!cen]; yy <- y[!cen]
    if (dist == "negbin") {
      d_eta[!cen] <- (yy - m) * th / (th + m)
      d_th[!cen] <- digamma(yy + th) - digamma(th) + log(th / (th + m)) + 1 -
        (yy + th) / (th + m)
    } else {
      d_eta[!cen] <- yy - m
    }
  }

  if (any(cen)) {
    m <- mu[cen]; k <- cc[cen] - 1
    if (dist == "negbin") {
      lS <- stats::pnbinom(k, size = th, mu = m, lower.tail = FALSE,
                           log.p = TRUE)
      lf <- stats::dnbinom(k, size = th, mu = m, log = TRUE)
      d_eta[cen] <- exp(lf - lS + log(m) + log(th + k) - log(th + m))
      J <- 0:max(k)
      lfJ <- outer(m, J, function(mm, jj)
        stats::dnbinom(jj, size = th, mu = mm, log = TRUE))
      G <- outer(m, J, function(mm, jj)
        digamma(jj + th) - digamma(th) + log(th / (th + mm)) + 1 -
          (jj + th) / (th + mm))
      d_th[cen] <- -rowSums(exp(lfJ - lS) * G * outer(k, J, ">="))
    } else {
      lS <- stats::ppois(k, m, lower.tail = FALSE, log.p = TRUE)
      d_eta[cen] <- exp(stats::dpois(k, m, log = TRUE) - lS + log(m))
    }
  }

  d_eta[!is.finite(d_eta)] <- 0
  g <- -drop(crossprod(x, d_eta))
  if (dist != "negbin") return(g)
  d_th[!is.finite(d_th)] <- 0
  c(g, -sum(d_th) * th)
}


## Fit the right-censored count model by maximum likelihood.
##
## Starting values come from the ordinary glm ON THE CENSORED DATA. That fit is
## biased downwards by construction -- it is the bias this model exists to
## remove -- but it is in the basin of attraction and costs one IRLS run, which
## is the same reasoning as in ztrunc.R.
##
## Field names follow the .countimp_1l_fit()/.countimp_zt_fit() contract
## exactly: beta, cov, scale, theta, ll, conv, nobs, npar, and NO class (B77).
.countimp_cs_fit <- function(x, y, cc, dist, max_grid = 5000L) {
  cen <- .countimp_cs_is_cen(y, cc)
  if (all(cen))
    stop("every observed value is at its censoring limit, so the data carry ",
         "no information\n  about the count distribution above it. A censored ",
         "model cannot be fitted here.", call. = FALSE)
  ## The theta gradient sums over 0..max(c)-1. A limit that far out cannot bind
  ## in any meaningful way, so the honest answer is to say the grid is the
  ## obstacle and name the two ways around it.
  if (dist == "negbin" && any(cen) && max(cc[cen]) > max_grid)
    stop("the censoring limit ", format(max(cc[cen])), " exceeds the grid ",
         "limit of ", max_grid, " used for the\n  dispersion gradient. Use ",
         "method \"cp\" (no dispersion parameter), or impute without ",
         "'censor'\n  if a limit that high does not bind.", call. = FALSE)

  st <- suppressWarnings(stats::glm.fit(x, y,
                                        family = stats::poisson()))$coefficients
  st[!is.finite(st)] <- 0
  if (dist == "negbin") st <- c(st, log(.countimp_mom_theta(y)))

  ## `...` of optim() reaches BOTH fn and gr, so the two must accept the same
  ## arguments. They did not while max_grid was a gradient argument, and the
  ## try() below reported that programming error as "could not fit the model" --
  ## a message that sends the user to look at their data (B43 is the same shape:
  ## a missing package appearing as a data problem). The optimiser's own message
  ## is therefore carried through now, and max_grid is checked above, where it
  ## belongs: it is a precondition of the fit, not an option of the gradient.
  op <- try(stats::optim(st, .countimp_cs_nll, .countimp_cs_gr,
                         x = x, y = y, cc = cc, dist = dist,
                         method = "BFGS", hessian = TRUE,
                         control = list(reltol = 1e-10, maxit = 500L)),
            silent = TRUE)
  if (inherits(op, "try-error") || op$value >= 1e10)
    stop("countimp: could not fit the right-censored ",
         if (dist == "negbin") "negative binomial" else "Poisson", " model",
         if (inherits(op, "try-error"))
           paste0(" (", sub("^Error[^:]*: *", "",
                            trimws(as.character(op))), ")")
         else "",
         ".",
         if (dist == "negbin") " Try method \"cp\" for this variable." else "",
         call. = FALSE)
  p <- ncol(x)
  list(beta = op$par[seq_len(p)],
       cov = .countimp_zt_cov(op$hessian, length(op$par))[seq_len(p),
                                                          seq_len(p),
                                                          drop = FALSE],
       scale = 1,
       theta = if (dist == "negbin") exp(op$par[p + 1L]) else NA_real_,
       ll = -op$value, conv = op$convergence, nobs = length(y),
       npar = length(op$par), dist = dist, ncens = sum(cen))
}


## Draw counts from the LATENT distribution, i.e. the one the censoring hides.
.countimp_cs_rlatent <- function(mu, dist, theta = NULL) {
  if (identical(dist, "negbin"))
    as.numeric(MASS::rnegbin(length(mu), mu = mu, theta = theta))
  else as.numeric(stats::rpois(length(mu), mu))
}


## Draw from the latent distribution conditioned on Y >= k + 1.
##
## This is what a censored OBSERVATION is worth once it is imputed: the case is
## known to lie at or above its limit, so the imputation must come from the left
## truncation of the fitted distribution at that limit, not from the whole of it.
##
## .countimp_rktp()/.countimp_rktnb() take a scalar k, so per-case limits are
## handled by grouping. With top coding there is one group and the loop runs
## once.
.countimp_cs_rleft <- function(mu, k, dist, theta = NULL) {
  out <- numeric(length(mu))
  for (kk in unique(k)) {
    i <- which(k == kk)
    out[i] <- if (identical(dist, "negbin"))
      .countimp_rktnb(length(i), size = theta, k = kk, mu = mu[i])
    else .countimp_rktp(length(i), k = kk, mu = mu[i])
  }
  out
}


## The single-level right-censored draw engine.
##
## Three kinds of cell can be asked for, and they are genuinely different
## questions -- this is where the method earns its keep:
##
##   missing (wy & !ry)                draw from the fitted distribution
##   censored observation (wy & ry     draw from that distribution CONDITIONED
##     & y >= c)                       on Y >= c: the value is known to be at
##                                     least c, and that information is kept
##   observed, not censored, but       draw from the fitted distribution; the
##     marked in `where`               user asked for a synthetic value
##
## WHICH SCALE the missing cells are drawn on is decided by whether censored
## observations are being imputed, because the completed data set must be on one
## scale:
##
##   no censored cell in `where` -> the completed data stay top-coded, so a
##     missing cell is drawn and then capped at its own limit. Imputing above
##     the ceiling would put values into a column that cannot contain them.
##   censored cells in `where` -> the completed data are on the latent scale,
##     so missing cells are drawn uncapped. Capping them here while resolving
##     the censored ones would mix two scales in one column.
##
## `latent = TRUE/FALSE` overrides the rule; the automatic choice is recorded
## through .countimp_note_event(), never taken silently.
##
## Note on the FCS cycle: a resolved censored cell is returned in imp, but the
## sampler writes only genuinely missing cells back into its working copy of the
## data (see zisampler.R). The censoring information therefore survives every
## iteration, which is what the fit needs -- and other variables that use this
## one as a predictor see the recorded, censored value. That is the conservative
## choice and it is stated in the method documentation.
.countimp_1l_censored <- function(y, ry, x, wy = NULL, dist, bayes,
                                  censor = NULL, latent = NULL, EV = FALSE, ...) {
  if (is.null(wy)) wy <- !ry
  cc <- .countimp_cs_parse(censor, n = length(y), y = y, ry = ry)
  x <- cbind(1, as.matrix(x))

  cen_obs <- ry & .countimp_cs_is_cen(y, cc)
  if (!any(cen_obs))
    .countimp_note_event("censor_never_binds",
                         "no observed value reaches its censoring limit; the fit equals the uncensored one")

  ## Which cells of the imputation target are resolved censored observations.
  ## Indexed within which(wy), because that is the order the return value is in.
  idx <- which(wy)
  res <- cen_obs[idx]
  lat <- if (is.null(latent)) any(res) else isTRUE(latent)
  if (is.null(latent))
    .countimp_note_event("censor_scale",
                         sprintf("%s scale chosen automatically (%d censored observation(s) in `where`)",
                                 if (lat) "latent" else "censored", sum(res)))

  ## The other direction, and the one that stays silent otherwise.
  ##
  ## `latent = TRUE` asks for the latent scale, but the scale a completed
  ## column carries is decided by `where`: only cells listed there are
  ## replaced. With no censored cell in `where`, the missing cells are drawn
  ## uncapped while the censored observations stay at their limit -- one column,
  ## two scales, and every summary of it is a mixture of the two.
  ##
  ## This is not hypothetical. A simulation study called it exactly this way and
  ## the error went unnoticed for weeks: coverage against the latent truth was
  ## 6.6 % instead of 91.7 %, and nothing in the run said why. The automatic
  ## choice was recorded, as documented -- but the choice that can go wrong is
  ## the one passed by argument, and that one was recorded nowhere.
  ##
  ## The combination is not forbidden: someone may deliberately want the
  ## missing cells latent and the censoring left in place. So it is a warning
  ## and a log entry, not an error -- but a warning, because a log entry alone
  ## demonstrably was not enough to make it visible.
  if (isTRUE(latent) && !any(res)) {
    n_cen <- sum(cen_obs)
    .countimp_note_event("censor_scale",
                         sprintf(paste("latent scale requested, but no censored cell is marked in",
                                       "`where`: the %d missing cell(s) are drawn uncapped while",
                                       "%d censored observation(s) stay at their limit, so the",
                                       "completed column carries two scales"),
                                 sum(wy), n_cen))
    if (!isTRUE(.countimp_state$mixed_scale_warned)) {
      warning("latent = TRUE was requested, but `where` marks no censored cell. ",
              "The imputed\n  cells are then drawn on the latent scale while the ",
              n_cen, " censored observation(s) stay at\n  their limit -- the ",
              "completed column carries two scales, and summaries of it mix\n  ",
              "them. To resolve the censoring, add the censored cells to ",
              "`where`; to keep it,\n  use latent = FALSE. See ?mice.impute.cp, ",
              "section \"Which scale the imputations are on\".\n  This warning is ",
              "shown once per session.", call. = FALSE)
      .countimp_state$mixed_scale_warned <- TRUE
    }
  }

  if (isFALSE(latent) && any(res))
    stop("latent = FALSE was requested, but ", sum(res), " censored ",
         "observation(s) are marked for\n  imputation in `where`. Resolving a ",
         "censored value means going above its limit, so\n  the two cannot be ",
         "combined: drop those cells from `where`, or use latent = TRUE.",
         call. = FALSE)

  ein_zug <- function() {
    if (isTRUE(bayes)) {
      f <- .countimp_cs_fit(x[ry, , drop = FALSE], y[ry], cc[ry], dist)
      ## beta is drawn, theta is not -- the same choice the zero-truncated and
      ## bounded families make, and for the same reason: it is register entry
      ## B36, whose effect was measured and found not to be systematic. Listed
      ## here so that the omission reads as a decision rather than an oversight.
      beta.star <- .countimp_1l_draw_beta(f$beta, f$cov, f$scale, bayes = TRUE)
    } else {
      obs <- which(ry)
      sel <- sample(obs, length(obs), replace = TRUE)
      f <- .countimp_cs_fit(x[sel, , drop = FALSE], y[sel], cc[sel], dist)
      beta.star <- f$beta
    }
    mu <- exp(drop(x[idx, , drop = FALSE] %*% beta.star))
    th <- if (is.na(f$theta)) NULL else f$theta

    im <- rep(NA_real_, length(mu))
    if (all(is.finite(mu))) {
      if (any(!res))
        im[!res] <- .countimp_cs_rlatent(mu[!res], dist, th)
      ## A censored cell is drawn from the left truncation at its own limit.
      if (any(res))
        im[res] <- .countimp_cs_rleft(mu[res], cc[idx][res] - 1, dist, th)
      ## On the censored scale the drawn latent value is what the recording
      ## would have produced: min(y*, c). This is NOT the clamping that
      ## bounds.R warns about -- there the support ends at hi and clamping
      ## misstates a truncated distribution, here the censoring genuinely maps
      ## the tail onto the limit and min() IS the observation process.
      if (!lat) im[!res] <- pmin(im[!res], cc[idx][!res])
    }

    if (isTRUE(EV) && all(is.finite(im))) {
      im <- .countimp_ev_screen(im, y, ry, x, wy)
      ## EV screening refills from an unbounded pmm draw, which knows nothing
      ## about the censoring; re-impose both properties it can break.
      if (all(is.finite(im))) {
        if (!lat) im <- pmin(im, cc[idx])
        low <- which(res & im < cc[idx])
        if (length(low))
          im[low] <- .countimp_cs_rleft(mu[low], cc[idx][low] - 1, dist, th)
      }
    }
    list(imp = im, fit = f, mu = mu)
  }

  .countimp_draw_retry(ein_zug, y_obs = y[ry],
                       method = if (dist == "negbin") "cnb" else "cp")
}
