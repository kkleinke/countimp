#' Single-level count imputation
#'
#' Imputes a count variable with missing values under a Poisson, quasi-Poisson
#' or negative binomial model. Each distribution comes in two variants: a
#' Bayesian one, which draws the regression coefficients from their sampling
#' distribution before drawing the counts, and a bootstrap one, which
#' bootstraps the observed cases and uses the resulting estimate. Both are
#' proper imputation in the sense of Rubin; the Bayesian variant is the
#' default choice, the bootstrap variant is useful when the coefficient
#' covariance is unstable.
#'
#' @param y Incomplete numeric count vector.
#' @param ry Logical vector: \code{TRUE} where \code{y} is observed.
#' @param x Numeric matrix or data frame of predictors, complete.
#' @param wy Logical vector: positions to impute. Defaults to \code{!ry}.
#' @param EV Logical. If \code{TRUE}, imputations flagged as extreme are
#'   redrawn by predictive mean matching. Requires package \pkg{mice}.
#' @param ... Ignored, present for compatibility with the method contract.
#'
#' @return A numeric vector of length \code{sum(wy)} with the imputed counts.
#'
#' @details
#' For the quasi-Poisson variants the estimated dispersion \eqn{\phi} enters
#' twice: it scales the covariance of the drawn coefficients, and it governs
#' the spread of the drawn counts. Earlier versions of this package scaled
#' only the counts and drew the coefficients from the unscaled covariance,
#' which understated parameter uncertainty by a factor \eqn{\sqrt{\phi}}.
#'
#' The practical effect is damped, because the error shrinks only the
#' between-imputation component \eqn{B} of Rubin's total variance
#' \eqn{T = \bar{U} + (1 + 1/m) B}. At a moderate missing-data rate
#' \eqn{\bar{U}} dominates and interval coverage stays close to nominal even
#' with the old code; the correction widens intervals by roughly four per cent
#' at 30 per cent missing values. The error matters where \eqn{B} carries more
#' of \eqn{T}, that is at high missing-data rates and with strongly
#' overdispersed data.
#'
#' \code{pois}, \code{pois.boot}, \code{qpois} and \code{qpois.boot} are
#' aliases for the correspondingly named \code{poisson} and
#' \code{quasipoisson} methods and behave identically.
#'
#' @examples
#' set.seed(1)
#' n <- 200
#' x <- rnorm(n)
#' y <- rpois(n, exp(1 + 0.5 * x))
#' d <- data.frame(y = y, x = x)
#' d$y[sample.int(n, 40)] <- NA
#' imp <- countimp(d, method = c("poisson", ""), m = 5, printFlag = FALSE)
#'
#' @author Kristian Kleinke
#' @name mice.impute.poisson
NULL

#' @describeIn mice.impute.poisson Bayesian Poisson variant
#' @export
mice.impute.poisson <- function(y, ry, x, wy = NULL, EV = FALSE, ...) {
  .countimp_1l_count(y, ry, x, wy = wy, dist = "poisson", bayes = TRUE,
                     EV = EV, ...)
}

#' @describeIn mice.impute.poisson Bootstrap Poisson variant
#' @export
mice.impute.poisson.boot <- function(y, ry, x, wy = NULL, EV = FALSE, ...) {
  .countimp_1l_count(y, ry, x, wy = wy, dist = "poisson", bayes = FALSE,
                     EV = EV, ...)
}

#' @describeIn mice.impute.poisson Bayesian quasi-Poisson variant
#' @export
mice.impute.quasipoisson <- function(y, ry, x, wy = NULL, EV = FALSE, ...) {
  .countimp_1l_count(y, ry, x, wy = wy, dist = "quasipoisson", bayes = TRUE,
                     EV = EV, ...)
}

#' @describeIn mice.impute.poisson Bootstrap quasi-Poisson variant
#' @export
mice.impute.quasipoisson.boot <- function(y, ry, x, wy = NULL, EV = FALSE, ...) {
  .countimp_1l_count(y, ry, x, wy = wy, dist = "quasipoisson", bayes = FALSE,
                     EV = EV, ...)
}

#' @describeIn mice.impute.poisson Bayesian negative binomial variant
#' @export
mice.impute.nb <- function(y, ry, x, wy = NULL, EV = FALSE, ...) {
  .countimp_1l_count(y, ry, x, wy = wy, dist = "negbin", bayes = TRUE,
                     EV = EV, ...)
}

#' @describeIn mice.impute.poisson Bootstrap negative binomial variant
#' @export
mice.impute.nb.boot <- function(y, ry, x, wy = NULL, EV = FALSE, ...) {
  .countimp_1l_count(y, ry, x, wy = wy, dist = "negbin", bayes = FALSE,
                     EV = EV, ...)
}

#' @describeIn mice.impute.poisson Alias for \code{mice.impute.poisson}
#' @export
mice.impute.pois <- function(y, ry, x, wy = NULL, EV = FALSE, ...) {
  mice.impute.poisson(y, ry, x, wy = wy, EV = EV, ...)
}

#' @describeIn mice.impute.poisson Alias for \code{mice.impute.poisson.boot}
#' @export
mice.impute.pois.boot <- function(y, ry, x, wy = NULL, EV = FALSE, ...) {
  mice.impute.poisson.boot(y, ry, x, wy = wy, EV = EV, ...)
}

#' @describeIn mice.impute.poisson Alias for \code{mice.impute.quasipoisson}
#' @export
mice.impute.qpois <- function(y, ry, x, wy = NULL, EV = FALSE, ...) {
  mice.impute.quasipoisson(y, ry, x, wy = wy, EV = EV, ...)
}

#' @describeIn mice.impute.poisson Alias for \code{mice.impute.quasipoisson.boot}
#' @export
mice.impute.qpois.boot <- function(y, ry, x, wy = NULL, EV = FALSE, ...) {
  mice.impute.quasipoisson.boot(y, ry, x, wy = wy, EV = EV, ...)
}
