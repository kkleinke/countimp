## ===========================================================================
## Base imputation methods for non-count variables (since 3.0.0)
##
## WHY THIS FILE EXISTS
##
## countimp implements imputation models for count data. A realistic data set,
## however, holds continuous, binary and categorical variables next to the
## counts, and those need imputation models too. Until 3.0.0 countimp obtained
## them from mice, which made mice a de-facto requirement for any mixed data
## set even though it is declared in Suggests.
##
## The methods below close that gap. They implement published algorithms --
## the same ones mice implements -- and they depend only on packages that ship
## with R (stats, MASS, nnet). countimp is therefore usable end to end without
## mice, while mice's methods remain available through the pass-through path
## (see .countimp_find_method).
##
## NAMING. These are exported as countimp.impute.<method>, NOT
## mice.impute.<method>. Using mice's prefix would mask mice's own functions
## whenever countimp is attached after mice, which would change results in
## code that does not mention countimp at all. The resolver knows both prefixes
## (see .countimp_find_method) and users keep writing method = "pmm".
##
## RELATION TO mice'S IMPLEMENTATIONS. norm, norm.nob, norm.predict, pmm,
## logreg, mean and sample follow the same algorithms and should be
## statistically interchangeable. polyreg and polr deliberately differ: they
## draw the regression coefficients from their posterior, which mice's versions
## do not do (they impute from the ML fit). See the note at those functions.
## Set options(countimp.methods = "mice") to prefer mice's implementations
## throughout.
##
## REFERENCES
##   Rubin, D. B. (1987) Multiple Imputation for Nonresponse in Surveys. Wiley.
##     Chapter 5.4 (the normal-model draw used by .countimp_norm_draw).
##   Little, R. J. A. (1988) Missing-data adjustments in large surveys.
##     Journal of Business and Economic Statistics 6, 287-296. (PMM)
##   White, I. R., Daniel, R., Royston, P. (2010) Avoiding bias due to perfect
##     prediction in multiple imputation of incomplete categorical variables.
##     Computational Statistics and Data Analysis 54, 2267-2275. (augmentation)
##   van Buuren, S. (2018) Flexible Imputation of Missing Data. 2nd ed. CRC.
##   Kleinke, K., Reinecke, J., Salfran, D., Spiess, M. (2020) Applied Multiple
##     Imputation. Springer.
## ===========================================================================


## ---------------------------------------------------------------------------
## Bayesian draw from the normal linear model
##
## Draws (beta, sigma) from their joint posterior under the standard
## non-informative prior (Rubin 1987, 5.4): sigma^2 from the residual sum of
## squares scaled by a chi-square draw, then beta given sigma^2.
##
## Numerically this goes through the QR decomposition rather than forming
## crossprod(x) and inverting it: with collinear predictors -- which happens
## routinely inside FCS, where every other variable enters as a predictor --
## the normal equations lose about twice as many digits as the QR route.
## A ridge penalty is applied only if the QR route reports rank deficiency,
## and it is reported rather than applied silently.
## ---------------------------------------------------------------------------

.countimp_norm_draw <- function(y, ry, x, ridge = 1e-05, ...) {
  xobs <- x[ry, , drop = FALSE]
  yobs <- y[ry]
  p    <- ncol(x)
  df   <- max(length(yobs) - p, 1L)

  fit  <- stats::lm.fit(x = xobs, y = yobs)
  cf   <- fit$coefficients
  R    <- qr.R(fit$qr)

  ## Rank deficiency shows up as NA coefficients (lm.fit pivots them out) or
  ## as a singular R. Both are handled the same way: penalise, and say so.
  singular <- anyNA(cf) || nrow(R) < p
  if (!singular) {
    V <- try(chol2inv(R), silent = TRUE)          # (x'x)^-1 without forming x'x
    singular <- inherits(V, "try-error")
  }
  if (singular) {
    xtx <- crossprod(xobs)
    V   <- solve(xtx + diag(diag(xtx) * ridge, nrow = p))
    cf[is.na(cf)] <- 0
    .countimp_note("collinear predictors: ridge penalty applied",
                     where = ".countimp_norm_draw")
  }

  sigma.star <- sqrt(sum(fit$residuals^2) / stats::rchisq(1L, df))
  ## V is symmetric by construction; enforce it so chol() cannot fail on
  ## accumulated floating-point asymmetry.
  V <- (V + t(V)) / 2
  L <- try(t(chol(V)), silent = TRUE)
  if (inherits(L, "try-error")) {
    ev <- eigen(V, symmetric = TRUE)
    L  <- ev$vectors %*% diag(sqrt(pmax(ev$values, 0)), nrow = p)
  }
  beta.star <- as.vector(cf) + as.vector(L %*% stats::rnorm(p)) * sigma.star

  list(coef = as.vector(cf), beta = beta.star, sigma = sigma.star)
}


## ---------------------------------------------------------------------------
## Data augmentation against perfect prediction (White, Daniel, Royston 2010)
##
## A categorical outcome that is perfectly predicted by a covariate drives the
## corresponding coefficient to infinity, and the imputations inherit that: all
## cases in the affected cell get the same category, with no uncertainty. The
## remedy is to append a small set of pseudo-observations covering every
## combination of category and covariate direction, down-weighted so that they
## contribute the equivalent of (p + 1) observations in total.
##
## The pseudo-rows are placed at mean +/- sd/2 of each predictor, clipped to the
## observed range so that augmentation never extrapolates beyond the data.
## ---------------------------------------------------------------------------

.countimp_augment <- function(y, ry, x, wy, maxcat = 50L) {
  x    <- as.matrix(x)
  p    <- ncol(x)
  cats <- sort(unique(unclass(y)))
  k    <- length(cats)
  if (k > maxcat)
    stop("Variable has ", k, " categories; the maximum for augmentation is ",
         maxcat, ".", call. = FALSE)

  ## Nothing to protect against without predictors, or with a single case to
  ## impute -- return unchanged rather than adding noise.
  if (p == 0L || sum(!ry) <= 1L)
    return(list(y = y, ry = ry, x = x, wy = wy, w = rep(1, length(y))))

  mu <- colMeans(x, na.rm = TRUE)
  sd <- apply(x, 2L, stats::sd, na.rm = TRUE)
  sd[!is.finite(sd)] <- 0
  lo <- apply(x, 2L, min, na.rm = TRUE)
  hi <- apply(x, 2L, max, na.rm = TRUE)

  nr <- 2L * p * k
  ## One block of 2*k rows per predictor: within a block that predictor is
  ## shifted by +/- sd/2 in alternation, the others sit at their mean.
  shift <- matrix(0, nr, p)
  for (j in seq_len(p)) {
    rows <- ((j - 1L) * 2L * k + 1L):(j * 2L * k)
    shift[rows, j] <- rep(c(0.5, -0.5), length.out = length(rows))
  }
  xa <- matrix(mu, nr, p, byrow = TRUE) + shift * matrix(sd, nr, p, byrow = TRUE)
  xa <- pmin(pmax(xa, matrix(lo, nr, p, byrow = TRUE)),
             matrix(hi, nr, p, byrow = TRUE))
  dimnames(xa) <- list(paste0("AUG", seq_len(nr)), colnames(x))

  ## Every category appears twice per predictor, so no category can be absent
  ## from the augmented data.
  ycat <- rep(rep(cats, each = 2L), times = p)
  ya <- if (is.factor(y)) {
    lv <- levels(y)
    if (is.ordered(y)) ordered(lv[c(unclass(y), ycat)], levels = lv)
    else               factor(lv[c(unclass(y), ycat)],  levels = lv)
  } else c(y, ycat)

  list(y  = ya,
       ry = c(ry, rep(TRUE,  nr)),
       x  = rbind(x, xa),
       wy = c(wy, rep(FALSE, nr)),
       w  = c(rep(1, length(y)), rep((p + 1L) / nr, nr)))
}


## ---------------------------------------------------------------------------
## Donor selection for predictive mean matching
##
## Returns, for each element of yhatmis, the index of one randomly chosen donor
## among its `donors` nearest neighbours in yhatobs (Little 1988).
##
## Implementation note: the k nearest neighbours of a point in a sorted vector
## always form a contiguous window of length k, and among all such windows the
## nearest-neighbour window is the one whose farther endpoint is closest to the
## target. Both facts together turn the search into a binary search plus k + 1
## vectorised comparisons, which avoids a per-case loop over the data -- FCS
## calls this function once per variable per iteration per imputation, so its
## cost is multiplied by maxit * m.
## ---------------------------------------------------------------------------

.countimp_match <- function(yhatobs, yhatmis, donors = 5L) {
  yhatobs <- as.vector(yhatobs)
  yhatmis <- as.vector(yhatmis)
  n1 <- length(yhatobs)
  n2 <- length(yhatmis)
  if (n1 == 0L) stop("No donors available: variable has no observed cases.",
                     call. = FALSE)
  if (n2 == 0L) return(integer(0))
  k <- max(1L, min(as.integer(donors), n1))
  if (k == n1)
    return(sample.int(n1, n2, replace = TRUE))

  ord    <- order(yhatobs)
  sorted <- yhatobs[ord]
  pos    <- findInterval(yhatmis, sorted)         # 0 .. n1

  ## Candidate window starts, clamped to a valid range; keep the best.
  best   <- pmin(pmax(pos - k + 1L, 1L), n1 - k + 1L)
  bestd  <- pmax(abs(sorted[best] - yhatmis), abs(sorted[best + k - 1L] - yhatmis))
  for (j in seq_len(k)) {
    s <- pmin(pmax(pos - k + 1L + j, 1L), n1 - k + 1L)
    d <- pmax(abs(sorted[s] - yhatmis), abs(sorted[s + k - 1L] - yhatmis))
    take <- d < bestd
    best[take]  <- s[take]
    bestd[take] <- d[take]
  }

  ord[best + sample.int(k, n2, replace = TRUE) - 1L]
}


## ---------------------------------------------------------------------------
## Continuous variables
## ---------------------------------------------------------------------------

#' Imputation of Continuous Variables under the Normal Linear Model
#'
#' Base imputation methods for continuous variables, so that \pkg{countimp}
#' can impute mixed data sets without \pkg{mice}.
#'
#' \describe{
#'   \item{\code{countimp.impute.norm}}{Bayesian linear regression: draws
#'     \eqn{(\beta, \sigma)} from their posterior and adds a residual draw.
#'     This is the method that propagates all sources of uncertainty.}
#'   \item{\code{countimp.impute.norm.nob}}{Linear regression ignoring the
#'     uncertainty about \eqn{(\beta, \sigma)}; keeps the residual draw. Too
#'     narrow for inference, useful as a reference.}
#'   \item{\code{countimp.impute.norm.predict}}{Fitted values without a
#'     residual draw. Deterministic; understates variability and attenuates
#'     associations. Provided for comparison, not for inference.}
#'   \item{\code{countimp.impute.pmm}}{Predictive mean matching: imputes an
#'     observed value from a donor with a similar predicted value. Keeps
#'     imputations inside the observed range and does not assume normality.}
#'   \item{\code{countimp.impute.mean}}{Unconditional mean. Distorts the
#'     distribution and biases nearly every estimate other than the mean
#'     itself; included for teaching and for reproducing analyses that used it.}
#'   \item{\code{countimp.impute.sample}}{Random draw from the observed values,
#'     ignoring all predictors.}
#' }
#'
#' @param y vector to be imputed.
#' @param ry logical vector, \code{TRUE} where \code{y} is observed.
#' @param x numeric matrix of predictors, without an intercept column.
#' @param wy logical vector, \code{TRUE} where an imputation is wanted.
#'   Defaults to \code{!ry}.
#' @param ridge ridge penalty, applied only if the predictors are rank
#'   deficient. Default \code{1e-05}.
#' @param donors number of candidate donors for \code{pmm}. Default \code{5}.
#' @param matchtype matching type for \code{pmm}: \code{0} matches predicted
#'   means from the point estimate, \code{1} (the default) matches observed
#'   predictions from the point estimate against missing predictions from the
#'   drawn coefficients, \code{2} uses drawn coefficients on both sides.
#' @param ... further arguments, ignored.
#' @return A vector of length \code{sum(wy)}.
#' @references Rubin, D. B. (1987) \emph{Multiple Imputation for Nonresponse in
#'   Surveys}. New York: Wiley.
#'
#'   Little, R. J. A. (1988) Missing-data adjustments in large surveys.
#'   \emph{Journal of Business and Economic Statistics}, 6, 287-296.
#'
#'   Kleinke, K., Reinecke, J., Salfran, D., Spiess, M. (2020)
#'   \emph{Applied Multiple Imputation}. Cham: Springer.
#' @name countimp.impute.norm
#' @export
countimp.impute.norm <- function(y, ry, x, wy = NULL, ridge = 1e-05, ...) {
  if (is.null(wy)) wy <- !ry
  x <- cbind(1, as.matrix(x))
  parm <- .countimp_norm_draw(y, ry, x, ridge = ridge)
  as.vector(x[wy, , drop = FALSE] %*% parm$beta) +
    stats::rnorm(sum(wy)) * parm$sigma
}

#' @rdname countimp.impute.norm
#' @export
countimp.impute.norm.nob <- function(y, ry, x, wy = NULL, ridge = 1e-05, ...) {
  if (is.null(wy)) wy <- !ry
  x <- cbind(1, as.matrix(x))
  parm <- .countimp_norm_draw(y, ry, x, ridge = ridge)
  ## Point estimate instead of the drawn coefficients: no parameter
  ## uncertainty. The residual draw uses the estimated sigma.
  sigma <- sqrt(sum((y[ry] - x[ry, , drop = FALSE] %*% parm$coef)^2) /
                  max(sum(ry) - ncol(x), 1L))
  as.vector(x[wy, , drop = FALSE] %*% parm$coef) + stats::rnorm(sum(wy)) * sigma
}

#' @rdname countimp.impute.norm
#' @export
countimp.impute.norm.predict <- function(y, ry, x, wy = NULL, ridge = 1e-05, ...) {
  if (is.null(wy)) wy <- !ry
  x <- cbind(1, as.matrix(x))
  parm <- .countimp_norm_draw(y, ry, x, ridge = ridge)
  as.vector(x[wy, , drop = FALSE] %*% parm$coef)
}

#' @rdname countimp.impute.norm
#' @export
countimp.impute.pmm <- function(y, ry, x, wy = NULL, donors = 5L,
                                matchtype = 1L, ridge = 1e-05, ...) {
  if (is.null(wy)) wy <- !ry
  x <- cbind(1, as.matrix(x))
  ## A factor y is matched on its integer codes. That treats the levels as
  ## equally spaced, which is defensible for ordered factors and crude for
  ## unordered ones -- use polyreg for the latter.
  ynum <- if (is.factor(y)) as.integer(y) else y
  parm <- .countimp_norm_draw(ynum, ry, x, ridge = ridge)

  cf <- switch(as.character(matchtype),
               "0" = list(parm$coef, parm$coef),
               "1" = list(parm$coef, parm$beta),
               "2" = list(parm$beta, parm$beta),
               stop("`matchtype` must be 0, 1 or 2.", call. = FALSE))
  yhatobs <- x[ry, , drop = FALSE] %*% cf[[1L]]
  yhatmis <- x[wy, , drop = FALSE] %*% cf[[2L]]
  y[ry][.countimp_match(yhatobs, yhatmis, donors)]
}

#' @rdname countimp.impute.norm
#' @export
countimp.impute.mean <- function(y, ry, x = NULL, wy = NULL, ...) {
  if (is.null(wy)) wy <- !ry
  rep(mean(y[ry]), sum(wy))
}

#' @rdname countimp.impute.norm
#' @export
countimp.impute.sample <- function(y, ry, x = NULL, wy = NULL, ...) {
  .countimp_impute_sample(y, ry, wy = wy)
}


## ---------------------------------------------------------------------------
## Categorical variables
## ---------------------------------------------------------------------------

#' Imputation of Categorical Variables
#'
#' Base imputation methods for binary, nominal and ordered variables, so that
#' \pkg{countimp} can impute mixed data sets without \pkg{mice}.
#'
#' \describe{
#'   \item{\code{countimp.impute.logreg}}{Binary variables. Logistic
#'     regression with coefficients drawn from their asymptotic posterior;
#'     imputations are Bernoulli draws from the resulting probabilities.}
#'   \item{\code{countimp.impute.polyreg}}{Unordered categorical variables.
#'     Multinomial logistic regression via \code{\link[nnet]{multinom}}.}
#'   \item{\code{countimp.impute.polr}}{Ordered categorical variables.
#'     Proportional-odds model via \code{\link[MASS]{polr}}, falling back to
#'     \code{multinom} if the proportional-odds fit fails.}
#' }
#'
#' All three protect against perfect prediction by data augmentation (White,
#' Daniel and Royston 2010): a small set of down-weighted pseudo-observations
#' covering every category keeps coefficients finite, so that a perfectly
#' predicted cell does not produce imputations without uncertainty.
#'
#' @section Parameter uncertainty in polyreg and polr:
#' These two functions draw the regression coefficients from their asymptotic
#' posterior \eqn{N(\hat\theta, H^{-1})} before computing the category
#' probabilities. Imputing from the fitted probabilities of the ML solution
#' instead -- which is what several implementations do -- omits the uncertainty
#' about \eqn{\theta} and produces between-imputation variance that is too
#' small, hence confidence intervals that are too narrow. Set
#' \code{bayes = FALSE} to obtain the ML behaviour, for instance to reproduce
#' results from such an implementation.
#'
#' @param y vector to be imputed; a factor for \code{polyreg} and \code{polr}.
#' @param ry logical vector, \code{TRUE} where \code{y} is observed.
#' @param x numeric matrix of predictors, without an intercept column.
#' @param wy logical vector, \code{TRUE} where an imputation is wanted.
#'   Defaults to \code{!ry}.
#' @param bayes draw the coefficients from their posterior (default
#'   \code{TRUE}) or impute from the ML fit (\code{FALSE}).
#' @param maxit maximum number of iterations passed to \code{multinom}.
#' @param ... further arguments, ignored.
#' @return A vector of length \code{sum(wy)}, a factor with the levels of
#'   \code{y} where \code{y} is a factor.
#' @references White, I. R., Daniel, R., Royston, P. (2010) Avoiding bias due
#'   to perfect prediction in multiple imputation of incomplete categorical
#'   variables. \emph{Computational Statistics and Data Analysis}, 54,
#'   2267-2275.
#'
#'   Kleinke, K., Reinecke, J., Salfran, D., Spiess, M. (2020)
#'   \emph{Applied Multiple Imputation}. Cham: Springer.
#' @name countimp.impute.logreg
#' @export
countimp.impute.logreg <- function(y, ry, x, wy = NULL, bayes = TRUE, ...) {
  if (is.null(wy)) wy <- !ry
  yf <- if (is.factor(y)) y else factor(y)
  if (nlevels(yf) > 2L)
    stop("`logreg` needs a variable with two categories; this one has ",
         nlevels(yf), ". Use `polyreg` or `polr`.", call. = FALSE)

  aug <- .countimp_augment(yf, ry, x, wy)
  X   <- cbind(1, as.matrix(aug$x))
  y01 <- as.integer(unclass(aug$y)) - 1L

  ## quasibinomial keeps glm.fit from warning about the fractional weights
  ## that augmentation introduces; it does not change the point estimate.
  fit <- stats::glm.fit(x = X[aug$ry, , drop = FALSE], y = y01[aug$ry],
                        weights = aug$w[aug$ry],
                        family = stats::quasibinomial(link = "logit"))
  beta <- fit$coefficients
  beta[is.na(beta)] <- 0

  if (bayes) {
    ## fit$qr is the QR of the IRLS-weighted design matrix, so R'R = X'WX and
    ## chol2inv(R) is the unscaled coefficient covariance -- no need to refit
    ## or to form and invert X'WX.
    R <- try(qr.R(fit$qr), silent = TRUE)
    V <- if (inherits(R, "try-error") || nrow(R) < ncol(X)) NULL else
           try(chol2inv(R), silent = TRUE)
    if (is.null(V) || inherits(V, "try-error")) {
      .countimp_note("singular logistic fit: imputing from the ML estimate",
                       where = "countimp.impute.logreg")
    } else {
      V <- (V + t(V)) / 2
      L <- try(t(chol(V)), silent = TRUE)
      if (inherits(L, "try-error")) {
        ev <- eigen(V, symmetric = TRUE)
        L  <- ev$vectors %*% diag(sqrt(pmax(ev$values, 0)), nrow = ncol(V))
      }
      beta <- beta + as.vector(L %*% stats::rnorm(length(beta)))
    }
  }

  p   <- stats::plogis(as.vector(X[aug$wy, , drop = FALSE] %*% beta))
  hit <- stats::runif(length(p)) <= p
  ## Numeric in, numeric out -- see .countimp_draw_category() for why.
  ## A 0/1 variable stored as numeric is the common case here, and returning a
  ## factor would put the level INDEX (1/2) into the data.
  if (!is.factor(y)) {
    z <- suppressWarnings(as.numeric(levels(yf)[1L + hit]))
    if (!any(is.na(z))) return(z)
  }
  factor(levels(yf)[1L + hit], levels = levels(yf))
}

#' @rdname countimp.impute.logreg
#' @export
countimp.impute.polyreg <- function(y, ry, x, wy = NULL, bayes = TRUE,
                                    maxit = 100L, ...) {
  if (is.null(wy)) wy <- !ry
  if (!requireNamespace("nnet", quietly = TRUE))
    stop("Method 'polyreg' needs package 'nnet' (shipped with R).",
         call. = FALSE)
  yf  <- as.factor(y)
  aug <- .countimp_augment(yf, ry, x, wy)
  fy  <- as.factor(aug$y)

  ## A category that covers every observed case leaves nothing to model.
  ## Same type contract as the draw below: numeric in, numeric out.
  tab <- table(aug$y[aug$ry])
  if (any(tab == sum(aug$ry))) {
    ein <- levels(fy)[which(tab == sum(aug$ry))]
    if (!is.factor(y)) {
      z <- suppressWarnings(as.numeric(ein))
      if (!is.na(z)) return(rep(z, sum(wy)))
    }
    return(factor(rep(ein, sum(wy)), levels = levels(yf)))
  }

  xy  <- data.frame(y = fy, as.data.frame(as.matrix(aug$x)))
  fit <- nnet::multinom(stats::formula(xy), data = xy[aug$ry, , drop = FALSE],
                        weights = aug$w[aug$ry], maxit = maxit, trace = FALSE,
                        Hess = TRUE, MaxNWts = 1500L)
  post <- .countimp_multinom_probs(fit, xy[aug$wy, , drop = FALSE],
                                   nlevels(fy), bayes)
  .countimp_draw_category(post, levels(fy), levels(yf),
                          want_numeric = !is.factor(y))
}

#' @rdname countimp.impute.logreg
#' @export
countimp.impute.polr <- function(y, ry, x, wy = NULL, bayes = TRUE,
                                 maxit = 100L, ...) {
  if (is.null(wy)) wy <- !ry
  yf  <- if (is.ordered(y)) y else factor(y, ordered = TRUE)
  aug <- .countimp_augment(yf, ry, x, wy)
  fy  <- aug$y
  xy  <- data.frame(y = fy, as.data.frame(as.matrix(aug$x)))

  fit <- try(suppressWarnings(
    MASS::polr(stats::formula(xy), data = xy[aug$ry, , drop = FALSE],
               weights = aug$w[aug$ry], Hess = TRUE)), silent = TRUE)

  if (inherits(fit, "try-error")) {
    ## The proportional-odds assumption can make the fit fail outright. A
    ## multinomial model is the weaker but always-fittable fallback; it is
    ## recorded so the user can see it happened.
    .countimp_note("proportional-odds fit failed: fell back to multinom",
                     where = "countimp.impute.polr")
    return(countimp.impute.polyreg(y, ry, x, wy = wy, bayes = bayes,
                                   maxit = maxit))
  }

  post <- .countimp_polr_probs(fit, xy[aug$wy, , drop = FALSE],
                              nlevels(fy), bayes)
  .countimp_draw_category(post, levels(fy), levels(yf),
                          want_numeric = !is.factor(y))
}


## --- helpers for the categorical methods -----------------------------------

## Category probabilities from a multinom fit, optionally after drawing the
## coefficients from N(theta.hat, H^-1). multinom stores the coefficients as a
## (nlev-1) x npar matrix; the Hessian is over the same parameters in row-major
## order, which is why the drawn vector is folded back with byrow = TRUE.
.countimp_multinom_probs <- function(fit, newdata, nlev, bayes) {
  if (bayes) {
    H <- fit$Hessian
    V <- if (is.null(H)) NULL else try(chol2inv(chol((H + t(H)) / 2)),
                                       silent = TRUE)
    if (is.null(V) || inherits(V, "try-error")) {
      .countimp_note("singular multinomial Hessian: imputing from the ML fit",
                       where = "polyreg")
    } else {
      cf  <- stats::coef(fit)
      mat <- if (is.matrix(cf)) cf else matrix(cf, nrow = 1L)
      L   <- t(chol(V))
      dr  <- as.vector(t(mat)) + as.vector(L %*% stats::rnorm(nrow(V)))
      X  <- stats::model.matrix(stats::formula(fit), data = newdata)
      B  <- matrix(dr, nrow = nrow(mat), byrow = TRUE)
      eta <- cbind(0, X %*% t(B))
      ex  <- exp(eta - apply(eta, 1L, max))
      return(ex / rowSums(ex))
    }
  }
  post <- stats::predict(fit, newdata, type = "probs")
  .countimp_as_prob_matrix(post, nrow(newdata), nlev)
}

## Category probabilities from a polr fit. The parameter vector is
## (beta, zeta); drawing it jointly from N(theta.hat, H^-1) respects the
## correlation between slopes and cut-points, which a separate draw would not.
.countimp_polr_probs <- function(fit, newdata, nlev, bayes) {
  if (bayes) {
    H <- fit$Hessian
    V <- if (is.null(H)) NULL else try(chol2inv(chol((H + t(H)) / 2)),
                                       silent = TRUE)
    if (is.null(V) || inherits(V, "try-error")) {
      .countimp_note("singular proportional-odds Hessian: imputing from the ML fit",
                       where = "polr")
    } else {
      nb <- length(fit$coefficients)
      th <- c(fit$coefficients, fit$zeta) +
        as.vector(t(chol(V)) %*% stats::rnorm(nrow(V)))
      X  <- stats::model.matrix(stats::formula(fit), data = newdata)
      X  <- X[, setdiff(colnames(X), "(Intercept)"), drop = FALSE]
      eta  <- as.vector(X %*% th[seq_len(nb)])
      zeta <- th[-seq_len(nb)]
      ## Cut-points must stay ordered; a draw can cross them.
      zeta <- sort(zeta)
      cum  <- vapply(zeta, function(z) stats::plogis(z - eta), numeric(length(eta)))
      cum  <- cbind(matrix(cum, nrow = length(eta)), 1)
      return(cbind(cum[, 1L, drop = FALSE],
                   cum[, -1L, drop = FALSE] - cum[, -ncol(cum), drop = FALSE]))
    }
  }
  post <- stats::predict(fit, newdata, type = "probs")
  .countimp_as_prob_matrix(post, nrow(newdata), nlev)
}

.countimp_as_prob_matrix <- function(post, n, nlev) {
  if (is.vector(post)) {
    post <- if (n == 1L) matrix(post, nrow = 1L) else matrix(c(1 - post, post), ncol = 2L)
  }
  post
}

## One draw per row from the categorical distribution in that row.
## `want_numeric`: what the CALLER handed in. A numeric column has to come
## back numeric -- a factor is inserted by the imputation loop through
## as.numeric(), and that yields the LEVEL INDEX rather than the value.
##
## Measured before this correction, y = 3 + Poisson(3), method = "polr":
## observed 3 to 12, imputed 1 to 9 -- every value shifted by the distance
## from the lower bound to 1, without a warning from the package. For y from 0
## it barely showed (a shift of one); for y from 10, the range 10 to 18 came
## back as 2 to 12. mice::mice() stops in the same situation with "'x' must be
## numeric or complex"; countimp ran on silently, which is the worse of the
## two outcomes.
.countimp_draw_category <- function(post, lev_fit, lev_out, want_numeric = FALSE) {
  post <- pmax(as.matrix(post), 0)
  post <- post / rowSums(post)
  u    <- stats::runif(nrow(post))
  idx  <- 1L + rowSums(t(apply(post, 1L, cumsum)) < u)
  idx  <- pmin(idx, ncol(post))
  if (want_numeric) {
    z <- suppressWarnings(as.numeric(lev_fit[idx]))
    ## Levels that will not convert to numbers would be a silent loss of
    ## data; the factor form is then the more honest result.
    if (!any(is.na(z))) return(z)
    .countimp_note("levels are not numeric: returning a factor",
                   where = ".countimp_draw_category")
  }
  factor(lev_fit[idx], levels = lev_out)
}
