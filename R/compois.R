## Conway-Maxwell-Poisson counts (underdispersion).
##
## Rate form:  P(Y = y) = lam^y / (y!)^nu / Z(lam, nu),  Z = sum_j lam^j/(j!)^nu
## nu > 1 underdisperses, nu = 1 is Poisson, nu < 1 overdisperses.
##
## The regression is in the MEAN parametrisation (Huang 2017): the linear
## predictor gives E[Y] = mu and lam is solved per case from (mu, nu). That is
## the parametrisation the mice.impute.* interface needs -- imputation draws for
## new cases from a fitted mean, and a rate-form intercept has no direct reading.
##
## Why this is written out here instead of calling glmmTMB(family = compois):
##   1. glmmTMB's simulate() IGNORES newdata. Measured: a fit on n = 3000 with
##      newdata of 500 rows returns 3000 values, i.e. draws for the FITTING
##      cases. Imputation needs draws for the cases being imputed, so the
##      generator has to be ours regardless.
##   2. glmmTMB reports 1/nu through sigma(), not nu. Measured against known
##      truth: lam = 10, nu = 2 gives sigma = 0.4924, i.e. 1/sigma = 2.031.
##      Depending on an undocumented reciprocal is exactly the fragility this
##      package is removing -- if an update flips it, countimp would draw with
##      inverted dispersion and the values would still look plausible.
## Both were verified against glmmTMB 1.1.14; the fit below reproduces its
## coefficients, nu and log-likelihood to four decimals on the same data.


## Grid window for the moment sums.
##
## Asymptotics (Shmueli et al. 2005): mode ~ lam^(1/nu), E[Y] ~ lam^(1/nu) -
## (nu-1)/(2 nu), Var ~ mu/nu. So the spread is sqrt(mu/nu), which for nu > 1 is
## NARROWER than Poisson: the window length follows the spread, not mu.
##
## lam^(1/nu) is computed as exp(log(lam)/nu) and capped: a grid running from 0
## would need lam^(1/nu) cells, which overflows to Inf for large mu and made
## 0:ymax abort with "result would be too long a vector".
.countimp_cmp_win <- function(lam, nu) {
  lmd <- log(lam) / nu
  md  <- exp(pmin(lmd, 700))
  list(md = md, sd = sqrt(pmax(md, 1e-8) / nu))
}


## Moments of the COM-Poisson on a windowed grid.
##
## One pass returns everything the log-likelihood and its gradient need:
## logZ, E[Y], Var(Y), E[log Y!] and Cov(Y, log Y!). The window is [md - kappa
## sd, md + kappa sd] per case, widened while either edge still carries mass.
##
## capped = TRUE means the parameters are NOT representable on a grid of at most
## max_w cells. That is reported rather than swallowed: past mu of roughly 1e15
## the mode is not distinguishable from mode + 1 and the moments would be noise.
.countimp_cmp_mom <- function(lam, nu, kappa = 12, min_w = 30L, max_w = 2e5,
                              max_cells = 2e6) {
  n <- length(lam)
  ## Non-representable input is reported here, not passed on: lam = Inf would
  ## make the window NaN and the caller would fail later at an unrelated line
  ## ("missing value where TRUE/FALSE needed" -- measured).
  if (!is.finite(nu) || nu <= 0 || any(!is.finite(lam)) || any(lam <= 0))
    return(list(logZ = rep(NA_real_, n), mu = rep(NA_real_, n),
                va = rep(NA_real_, n), elf = rep(NA_real_, n),
                cov = rep(NA_real_, n), W = NA_integer_, capped = TRUE))

  loglam <- log(lam)
  w0 <- .countimp_cmp_win(lam, nu)
  ## Width spans BOTH sides. A width of kappa*sd covered only the left half and
  ## the widening loop then fixed it after the fact, on a needlessly large grid.
  ##
  ## The width is checked BEFORE anything is allocated. At nu = 0.0023 with
  ## mu = 257 the mode is lam^(1/nu) = 1.5e71 and the required width 1.9e38 --
  ## an earlier version clamped W to max_w and built the matrix anyway, so a
  ## fit on overdispersed data allocated hundreds of MB per likelihood call and
  ## did not finish. Reporting capped here lets the caller return 1e300 and the
  ## optimiser walk back out of that region.
  need <- 2 * kappa * max(w0$sd)
  if (!is.finite(need) || need > max_w)
    return(list(logZ = rep(NA_real_, n), mu = rep(NA_real_, n),
                va = rep(NA_real_, n), elf = rep(NA_real_, n),
                cov = rep(NA_real_, n), W = NA_integer_, capped = TRUE))

  ## Grouped by REQUIRED width, then blocked by a cell budget.
  ##
  ## The grid is rectangular, so one block costs (rows x widest row in it). A
  ## single block over all cases therefore pays the widest row for every row:
  ## measured at nu = 0.144, n = 400, the widths ranged over 73 cells but the
  ## block was 282 wide -- 113 000 cells for 30 000 of content, and the whole
  ## COM-Poisson fit inside countimp_fit_diag() took 4.2 of its 4.6 seconds.
  ## Sorting by width and cutting at doublings keeps the waste inside a block
  ## below a factor of two.
  ##
  ## The cell budget then bounds peak memory by max_cells independently of n,
  ## so a wide grid costs time rather than the session.
  wi  <- pmax(min_w, ceiling(2 * kappa * w0$sd))
  ord <- order(wi)
  grp <- cumsum(c(TRUE, diff(log2(wi[ord])) > 1))
  out <- list(logZ = numeric(n), mu = numeric(n), va = numeric(n),
              elf = numeric(n), cov = numeric(n))
  capped <- FALSE
  Wmax <- 0L
  for (g in unique(grp)) {
    gi <- ord[grp == g]
    W  <- max(min_w, max(wi[gi]))
    bs <- max(1L, min(length(gi), floor(max_cells / (W + 1L))))
    for (st in seq(1L, length(gi), by = bs)) {
      id <- gi[st:min(st + bs - 1L, length(gi))]
      b  <- .countimp_cmp_mom_blk(loglam[id], nu, w0$md[id], w0$sd[id],
                                  kappa, W, max_w)
      for (nm in names(out)) out[[nm]][id] <- b[[nm]]
      capped <- capped || b$capped
      Wmax <- max(Wmax, b$W)
    }
  }
  c(out, list(W = Wmax, capped = capped))
}


## One block of cases. Split out of .countimp_cmp_mom() so the widening loop
## works on a bounded matrix.
.countimp_cmp_mom_blk <- function(loglam, nu, md, sd, kappa, W, max_w) {
  n <- length(loglam)
  capped <- FALSE
  lo <- pmax(0, floor(md - kappa * sd))
  repeat {
    W <- min(W, max_w)
    J  <- matrix(rep.int(0:W, rep.int(n, W + 1L)), n, W + 1L) + lo
    ## lf is kept: it is needed again below for E[log y!] and cov, and lgamma
    ## over the whole grid is the single most expensive operation here.
    ## Computing it twice was 26% of the fit (Rprof, attributed to
    ## is.data.frame -- rowSums forces its lazily evaluated argument there).
    lf <- lgamma(J + 1)
    LP <- J * loglam - nu * lf
    ## Row maxima via max.col(). do.call(pmax, as.data.frame(LP)) turns a wide
    ## matrix into a list with one element PER COLUMN -- at W = 2e5 that alone
    ## hung the fit.
    rm_ <- LP[cbind(seq_len(n), max.col(LP, ties.method = "first"))]
    ok_l <- lo == 0 | LP[, 1L] - rm_ < -55
    ok_r <- LP[, W + 1L] - rm_ < -55
    if (all(ok_l) && all(ok_r)) break
    if (W >= max_w) { capped <- TRUE; break }
    if (!all(ok_r)) W <- ceiling(W * 1.7) + 20L
    if (!all(ok_l)) lo <- pmax(0, lo - ceiling(0.7 * W) - 10)
  }
  P  <- exp(LP - rm_)
  sw <- rowSums(P)
  P  <- P / sw
  ## PJ reused for mu, va and cov: three full-grid products instead of five.
  PJ  <- P * J
  mu  <- rowSums(PJ)
  elf <- rowSums(P * lf)
  list(logZ = rm_ + log(sw), mu = mu, va = rowSums(PJ * J) - mu^2,
       elf = elf, cov = rowSums(PJ * lf) - mu * elf,
       W = W, capped = capped)
}


## Solve lam from (mu, nu) by Newton on t = log lam.
##
## d mu / d t = Var(Y) > 0, so the step is (mu - mu(t)) / Var and the map is
## monotone -- no bracketing needed. Steps are clipped to +-5 in t, and t itself
## to +-700, because the Huang (2017) start value nu * log(mu + (nu-1)/(2 nu))
## has nu in the exponent: at nu = 1e150 exp(t) would be Inf.
.countimp_cmp_lam <- function(mu, nu, start = NULL, tol = 1e-11, maxit = 80L) {
  if (is.null(start)) {
    ## Huang (2017): lam ~ (mu + (nu-1)/(2 nu))^nu. The bracket turns NEGATIVE
    ## for nu < 1/(1 + 2 mu) -- at mu = 6, nu = 0.05 it is 6 - 9.5 -- so that
    ## branch needs its own start. As nu -> 0 the distribution is geometric with
    ## mean lam/(1-lam), giving lam = mu/(1+mu), which is exact in the limit and
    ## close throughout the small-nu range (0.857 against 0.877 at nu = 0.01).
    br <- mu + (nu - 1) / (2 * nu)
    t <- ifelse(br > 0, nu * log(pmax(br, 1e-300)), log(mu / (1 + mu)))
  } else t <- start
  t <- pmax(pmin(t, 700), -700)
  ## The step is limited in the LOG MODE, not in log lam: the mode is
  ## lam^(1/nu), so a step of 5 in t moves it by exp(5/nu) -- at nu = 0.02 that
  ## is exp(250). Measured before this bound: the solve for mu = 6 at nu = 0.05
  ## returned lam = 59, whose actual mean is 3e5. The limit is also capped at 5
  ## in t so that large nu does not take unbounded steps in lam.
  lim <- 5 * min(1, nu)
  m <- NULL
  for (i in seq_len(maxit)) {
    m <- .countimp_cmp_mom(exp(t), nu)
    if (isTRUE(m$capped)) break
    d <- (mu - m$mu) / pmax(m$va, 1e-12)
    d[!is.finite(d)] <- 0
    t <- pmax(pmin(t + pmax(pmin(d, lim), -lim), 700), -700)
    if (max(abs(m$mu - mu)) < tol * max(1, max(mu))) break
  }
  list(lam = exp(t), t = t, iter = i, mom = m)
}


## Beyond these the distribution is degenerate: nu -> Inf is a point mass,
## nu -> 0 has no finite moments. Without a bound BFGS jumped from the start
## value nu = 1 (gradient there was -366) straight to zeta = 366 in one step.
.countimp_cmp_nu_min <- 1e-3
.countimp_cmp_nu_max <- 1e3
.countimp_cmp_mu_max <- 1e7


## Negative log-likelihood and analytic gradient in (beta, zeta = log nu).
##
## Derivation -- every quantity comes from ONE grid pass:
##   l_i = y_i log lam_i - nu log(y_i!) - log Z(lam_i, nu),  mu_i = exp(x_i'b)
##   d l_i / d log lam_i      = y_i - mu_i
##   d log lam_i / d log mu_i = mu_i / V_i        (because d mu / d log lam = V)
##   => d l / d beta = sum_i (y_i - mu_i) (mu_i / V_i) x_i
##      At nu = 1, mu/V = 1 and this is exactly the Poisson score. That is the
##      independent check: the two agree to 2e-14.
##   For nu, lam has to move because mu is held fixed:
##   d log lam_i / d nu   = cov_i / V_i,  cov_i = Cov(Y_i, log Y_i!)
##   d log Z / d nu |_lam = -E[log Y!]  = -elf_i
##   => d l_i / d nu   = (y_i - mu_i) cov_i / V_i - log(y_i!) + elf_i
##      d l   / d zeta = nu * sum_i d l_i / d nu
##
## Returns 1e300 / a zero gradient outside the representable region so that
## optim() walks back instead of marching into parameters whose moments would
## be nonsense.
.countimp_cmp_nll <- function(par, x, y, gr = FALSE, memo = NULL) {
  p <- ncol(x)
  raus <- function() if (gr) rep(0, p + 1L) else 1e300
  if (!all(is.finite(par))) return(raus())
  nu <- exp(par[p + 1L])
  mu <- exp(drop(x %*% par[seq_len(p)]))
  if (nu < .countimp_cmp_nu_min || nu > .countimp_cmp_nu_max ||
      any(!is.finite(mu)) || any(mu <= 0) || max(mu) > .countimp_cmp_mu_max)
    return(raus())
  ## optim() calls fn and gr separately, and the mean solve dominates the cost
  ## (it iterates the moment sum; 8 iterations measured at nu = 0.14). Without
  ## this cache a par visited by both is solved twice.
  ##
  ## Measured, and the honest result is that this buys little: on the
  ## countimp_fit_diag() example (n = 400, nu = 0.14) optim made 41 fn calls
  ## and 14 gr calls in 4.71 s. All 14 gr calls hit the cache, but they are a
  ## quarter of the total, so the saving is small -- the win in that example
  ## came from not computing lgamma twice (see .countimp_cmp_mom_blk). Kept
  ## because it is correct and free; do not expect it to matter.
  ##
  ## The cache holds ONE point, which is all optim needs; passing it in as an
  ## argument keeps it per-fit rather than package-level state.
  if (!is.null(memo)) {
    if (!is.null(memo$par) && identical(memo$par, par)) {
      m <- memo$m
    } else {
      m <- .countimp_cmp_lam(mu, nu)
      memo$par <- par
      memo$m <- m
    }
  } else m <- .countimp_cmp_lam(mu, nu)
  if (isTRUE(m$mom$capped) || any(!is.finite(m$mom$logZ))) return(raus())
  lf <- lgamma(y + 1)
  if (!gr) {
    ll <- y * log(m$lam) - nu * lf - m$mom$logZ
    return(if (all(is.finite(ll))) -sum(ll) else 1e300)
  }
  w <- (y - m$mom$mu) / pmax(m$mom$va, 1e-300)
  g <- -c(drop(crossprod(x, w * m$mom$mu)),
          nu * sum(w * m$mom$cov - lf + m$mom$elf))
  if (all(is.finite(g))) g else rep(0, p + 1L)
}


## Fit the COM-Poisson regression.
##
## Field names follow the .countimp_1l_fit()/.countimp_zt_fit() contract
## exactly -- beta, cov, scale, theta, ll, conv, nobs, npar -- and the object
## carries NO class, like the other two. Both are deliberate (see B77).
## theta holds nu, which is what countimp_fit_diag() and the draw read.
##
## par_full/cov_full are extra: they carry zeta as well, so the Bayesian draw
## can take (beta, zeta) JOINTLY from the asymptotic posterior instead of
## holding the dispersion at its point estimate.
##
## The start value for zeta comes from the MEASURED dispersion (nu ~ 1/phi,
## because V/M ~ 1/nu), not from nu = 1: that start produced a gradient of -366
## and a wild first step.
.countimp_cmp_fit <- function(x, y) {
  g0  <- suppressWarnings(stats::glm.fit(x, y, family = stats::poisson()))
  mu0 <- stats::fitted(g0)
  phi <- sum((y - mu0)^2 / pmax(mu0, 1e-8)) / max(1L, length(y) - ncol(x))
  nu0 <- min(max(1 / max(phi, 1e-6), .countimp_cmp_nu_min * 10),
             .countimp_cmp_nu_max / 10)
  st  <- c(g0$coefficients, log(nu0))
  st[!is.finite(st)] <- 0

  ## One environment per fit, shared by fn and gr so the mean solve at a given
  ## par happens once instead of twice. NOT package-level state: a second fit
  ## gets a fresh one, and nothing survives the call.
  memo <- new.env(parent = emptyenv())
  op <- stats::optim(st, .countimp_cmp_nll,
                     function(par, ...) .countimp_cmp_nll(par, ..., gr = TRUE),
                     x = x, y = y, memo = memo, method = "BFGS", hessian = TRUE,
                     control = list(reltol = 1e-10, maxit = 500L))
  p  <- ncol(x)
  cv <- .countimp_zt_cov(op$hessian, p + 1L)
  list(beta = op$par[seq_len(p)],
       cov = cv[seq_len(p), seq_len(p), drop = FALSE],
       scale = 1, theta = exp(op$par[p + 1L]),
       ll = -op$value, conv = op$convergence, nobs = length(y),
       npar = length(op$par), par_full = op$par, cov_full = cv)
}


## Draw from the COM-Poisson with given means.
##
## Exact categorical draw on the log-scale window grid, the same construction
## the bounded family uses: an inverse-cdf draw underflows when the mass sits
## far from where the cdf is evaluated. Cases are sorted by mu and taken in
## blocks so each block gets a tight window rather than one grid sized for
## max(mu).
.countimp_rcmp <- function(mu, nu, block = 500L) {
  n <- length(mu)
  if (n == 0L) return(numeric(0))
  out <- numeric(n)
  o <- order(mu)
  for (s in seq(1L, n, by = block)) {
    idx <- o[s:min(s + block - 1L, n)]
    lam <- .countimp_cmp_lam(mu[idx], nu)$lam
    if (any(!is.finite(lam))) { out[idx] <- NA_real_; next }
    w0 <- .countimp_cmp_win(lam, nu)
    ## Same pre-check as in .countimp_cmp_mom(): a width that cannot be built is
    ## reported as NA rather than clamped and allocated.
    if (!is.finite(2 * 12 * max(w0$sd)) || 2 * 12 * max(w0$sd) > 2e5) {
      out[idx] <- NA_real_; next
    }
    W  <- max(30L, ceiling(24 * max(w0$sd)))
    lo <- pmax(0, floor(w0$md - 12 * w0$sd))
    repeat {
      W <- min(W, 2e5)
      J  <- matrix(rep.int(0:W, rep.int(length(idx), W + 1L)),
                   length(idx), W + 1L) + lo
      LP <- J * log(lam) - nu * lgamma(J + 1)
      rm_ <- LP[cbind(seq_len(nrow(LP)), max.col(LP, ties.method = "first"))]
      ok_l <- lo == 0 | LP[, 1L] - rm_ < -55
      ok_r <- LP[, W + 1L] - rm_ < -55
      if ((all(ok_l) && all(ok_r)) || W >= 2e5) break
      if (!all(ok_r)) W <- ceiling(W * 1.7) + 20L
      if (!all(ok_l)) lo <- pmax(0, lo - ceiling(0.7 * W) - 10)
    }
    CS <- .countimp_rowcumsum(exp(LP - rm_))
    u  <- stats::runif(length(idx)) * CS[, ncol(CS)]
    out[idx] <- J[cbind(seq_along(idx),
                        max.col(CS >= u, ties.method = "first"))]
  }
  out
}


## Row-wise cumulative sums. t(apply(W, 1, cumsum)) transposes twice and
## silently drops to a vector when the matrix has one row.
.countimp_rowcumsum <- function(M) {
  if (nrow(M) == 1L) return(matrix(cumsum(M[1L, ]), nrow = 1L))
  matrix(t(apply(M, 1L, cumsum)), nrow(M), ncol(M))
}


## The single-level COM-Poisson draw engine.
##
## Same shape as .countimp_1l_count() and .countimp_1l_bounded(): one closure
## per draw handed to .countimp_draw_retry(), so a hard failure is repeated
## rather than fatal.
.countimp_1l_cmp <- function(y, ry, x, wy = NULL, bayes = TRUE, EV = FALSE, ...) {
  if (is.null(wy)) wy <- !ry
  x <- cbind(1, as.matrix(x))

  ein_zug <- function() {
    if (isTRUE(bayes)) {
      f <- .countimp_cmp_fit(x[ry, , drop = FALSE], y[ry])
      ## Joint draw of (beta, zeta): the dispersion is a parameter with
      ## posterior uncertainty like any other, and the Hessian already gives
      ## its covariance with beta. Holding nu at its point estimate would make
      ## the imputations too similar across m.
      ps <- .countimp_1l_draw_beta(f$par_full, f$cov_full, f$scale, bayes = TRUE)
      beta.star <- ps[seq_along(f$beta)]
      nu.star <- exp(ps[length(ps)])
    } else {
      obs <- which(ry)
      sel <- sample(obs, length(obs), replace = TRUE)
      f <- .countimp_cmp_fit(x[sel, , drop = FALSE], y[sel])
      beta.star <- f$beta
      nu.star <- f$theta
    }
    nu.star <- min(max(nu.star, .countimp_cmp_nu_min), .countimp_cmp_nu_max)
    mu <- exp(drop(x[wy, , drop = FALSE] %*% beta.star))
    im <- if (any(!is.finite(mu)) || max(mu) > .countimp_cmp_mu_max)
      rep(NA_real_, length(mu)) else .countimp_rcmp(mu, nu.star)
    if (isTRUE(EV) && all(is.finite(im))) im <- .countimp_ev_screen(im, y, ry, x, wy)
    list(imp = im, fit = f, mu = mu)
  }

  .countimp_draw_retry(ein_zug, y_obs = y[ry], method = "cmp")
}


## ---------------------------------------------------------------------------
## Contact with glmmTMB's compois family (used by the two-level methods)
##
## The single-level methods above deliberately avoid glmmTMB. The two-level ones
## cannot: random effects need the Laplace machinery, which is why every 2l
## method in this package goes through glmmTMB. What CAN be kept is the reason
## the single-level core exists -- not depending on an undocumented convention
## for what the reported dispersion MEANS.
##
## Measured against glmmTMB 1.1.14, 500 cases in 25 clusters, data drawn from
## this file's own generator at nu = 2.5 (analyse/k25_2lcmp_glmmtmb.csv):
##
##   fixef            1.2984, 0.3899   (true 1.3, 0.4)
##   sqrt(VarCorr)    0.3457           (true 0.35)
##   sigma(fit)       0.4025           -> 1/sigma = 2.4845, i.e. sigma = 1/nu
##   fit$fit$par      names beta, betadisp, theta; exp(betadisp) = 0.4025,
##                    so the INTERNAL dispersion parameter is log(sigma) and
##                    the joint draw in .countimp_draw_2l() perturbs the right
##                    quantity on the right scale
##   predict(response) equals exp(predict(link)) to 0 (max absolute difference
##                    over 5 rows), i.e. the fit is in the mean parametrisation,
##                    like the core above
##   the fit itself   2.3 s at n = 500, which is what makes these methods
##                    several times more expensive than 2l.poisson
##
## So sigma() = 1/nu today. The point of the function below is that countimp
## does not have to believe that tomorrow.


## Conditional log-likelihood of the COM-Poisson at GIVEN means.
##
## Used to compare two candidate shapes against the same data; the mean solve is
## the same one the fit uses, so this shares the core rather than duplicating it.
.countimp_cmp_ll_mu <- function(y, mu, nu) {
  if (!is.finite(nu) || nu <= 0) return(-Inf)
  ok <- is.finite(mu) & mu > 0
  if (!any(ok)) return(-Inf)
  m <- .countimp_cmp_lam(mu[ok], nu)
  if (isTRUE(m$mom$capped) || any(!is.finite(m$mom$logZ))) return(-Inf)
  ll <- y[ok] * log(m$lam) - nu * lgamma(y[ok] + 1) - m$mom$logZ
  if (all(is.finite(ll))) sum(ll) else -Inf
}


## The shape nu of a glmmTMB compois fit -- MEASURED, not assumed.
##
## glmmTMB reports 1/nu through sigma(). A future release could flip that
## without anything looking wrong: the imputations would still be plausible
## counts, drawn at inverted dispersion. The register calls that class of defect
## the dangerous one, so the orientation is decided from the data instead:
##
##   both candidates, s and 1/s, are scored by the conditional log-likelihood of
##   the observed counts at the fitted means, and the better one wins.
##
## The comparison is legitimate because mu is held fixed at the fit's own values
## (which already contain the random effects) and only the shape varies -- so
## this is a one-parameter likelihood comparison inside countimp's own,
## independently verified core (B82), not a second opinion from the same source.
##
## Two guards:
##   * |log s| < 0.05 means the two candidates differ by under 5%, so the choice
##     cannot matter; the documented reading is taken without a note.
##   * a disagreement with the documented reading is REPORTED through the
##     diagnostics log. It is the signal that glmmTMB changed, and it should
##     never pass unseen even though countimp handles it correctly.
##
## Measured on the fit above (analyse/k25_2lcmp_lesart.csv):
## ll(nu = 1/sigma = 2.4845) = -814.6 against ll(nu = sigma = 0.4025) =
## -1011.6, a gap of 197 -- the two candidates are not close, so the rule is
## not deciding on noise.
## Returns list(nu, rezi, sigma): the shape, and WHICH READING produced it.
## The reading is returned rather than recomputed by the caller, because a
## caller that infers it from nu alone gets it wrong in exactly the cases where
## it matters -- an unreadable sigma, or sigma numerically equal to 1.
.countimp_cmp_nu_glmm <- function(fit, y, mu) {
  s <- suppressWarnings(as.numeric(.countimp_theta(fit, pflicht = FALSE) %||% NA_real_))
  if (!is.finite(s) || s <= 0) {
    .countimp_note_event("compois_sigma_unreadable",
                         "no usable dispersion from the fit; nu set to 1 (Poisson)")
    return(list(nu = 1, rezi = TRUE, sigma = NA_real_))
  }
  if (abs(log(s)) < 0.05) return(list(nu = 1 / s, rezi = TRUE, sigma = s))
  kand <- c(1 / s, s)
  ll <- vapply(kand, function(nu) .countimp_cmp_ll_mu(y, mu, nu), numeric(1))
  if (all(!is.finite(ll))) return(list(nu = 1 / s, rezi = TRUE, sigma = s))
  best <- which.max(ll)
  if (best != 1L)
    .countimp_note_event("compois_sigma_orientation",
      sprintf(paste0("glmmTMB %s reports a dispersion that fits as nu = sigma, ",
                     "not nu = 1/sigma (ll %.1f vs %.1f); countimp used the ",
                     "measured orientation"),
              .countimp_pkg_version("glmmTMB"), ll[2], ll[1]))
  list(nu = kand[best], rezi = best == 1L, sigma = s)
}
