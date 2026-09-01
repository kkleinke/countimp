#' Zero-truncated count imputation
#'
#' Imputes a count variable whose values are structurally at least one, under a
#' zero-truncated Poisson (\code{ztp}) or zero-truncated negative binomial
#' (\code{ztnb}) model. Both come in a Bayesian variant, which draws the
#' regression coefficients from their sampling distribution before drawing the
#' counts, and a bootstrap variant (suffix \code{.boot}), which bootstraps the
#' observed cases and uses the resulting estimate.
#'
#' @details
#' Use these methods when zero is impossible by design rather than merely
#' unobserved: length of stay among admitted patients, number of purchases
#' among buyers, litter size among mothers who gave birth. Imputing such a
#' variable with \code{"poisson"} or \code{"nb"} produces zeros the population
#' cannot contain, and biases mu downward for every case, because the model has
#' to place probability mass on a value the sampling scheme excluded.
#'
#' The distinction from the hurdle and zero-inflated methods matters and is not
#' a matter of taste:
#'
#' \itemize{
#'   \item \code{ztp}/\code{ztnb} -- zeros are \emph{impossible}. The observed
#'     data contain none, and none may be imputed.
#'   \item \code{hp}/\code{hnb} (hurdle) -- zeros are possible and arise from a
#'     separate process. Observed zeros are present.
#'   \item \code{zip}/\code{zinb} (zero-inflated) -- zeros arise from two
#'     sources, a structural one and the count distribution itself.
#' }
#'
#' Passing observed zeros to \code{ztp} or \code{ztnb} is an error, not a
#' warning: it means the model is the wrong one for these data, and the message
#' names the alternatives.
#'
#' The models are fitted by maximum likelihood with analytic gradients, using
#' only \pkg{stats}. Unlike the multilevel hurdle methods, they need neither
#' \pkg{glmmTMB} nor a C++ toolchain. Accuracy against
#' \code{glmmTMB::truncated_poisson} and \code{glmmTMB::truncated_nbinom2} on
#' the same data is at optimiser tolerance.
#'
#' @param y Incomplete numeric count vector, observed values all >= 1.
#' @param ry Logical vector: \code{TRUE} where \code{y} is observed.
#' @param x Numeric matrix or data frame of predictors, complete.
#' @param wy Logical vector: positions to impute. Defaults to \code{!ry}.
#' @param EV Logical. If \code{TRUE}, imputations flagged as extreme are
#'   redrawn by predictive mean matching. Requires package \pkg{mice}.
#' @param ... Ignored, present for compatibility with the method contract.
#'
#' @return A numeric vector of length \code{sum(wy)} with the imputed counts,
#'   all of them >= 1.
#'
#' @seealso \code{\link{countimp_fit_diag}} to check which count model the
#'   observed data support; \code{\link{mice.impute.hp}} when zeros are present.
#'
#' @examples
#' set.seed(1)
#' n  <- 200
#' x1 <- rnorm(n)
#' mu <- exp(0.8 + 0.5 * x1)
#' y  <- rpois(n, mu)
#' y[y == 0] <- 1                      # zero-truncated by construction
#' ry <- rbinom(n, 1, 0.8) == 1
#' y[!ry] <- NA
#' imp <- mice.impute.ztp(y, ry, data.frame(x1 = x1))
#' range(imp)                          # never zero
#'
#' @name mice.impute.ztp
#' @rdname mice.impute.ztp
#' @export
mice.impute.ztp <- function(y, ry, x, wy = NULL, EV = FALSE, ...) {
  .countimp_1l_ztrunc(y, ry, x, wy = wy, dist = "ztp", bayes = TRUE,
                      EV = EV, ...)
}

#' @rdname mice.impute.ztp
#' @export
mice.impute.ztp.boot <- function(y, ry, x, wy = NULL, EV = FALSE, ...) {
  .countimp_1l_ztrunc(y, ry, x, wy = wy, dist = "ztp", bayes = FALSE,
                      EV = EV, ...)
}

#' @rdname mice.impute.ztp
#' @export
mice.impute.ztnb <- function(y, ry, x, wy = NULL, EV = FALSE, ...) {
  .countimp_1l_ztrunc(y, ry, x, wy = wy, dist = "ztnb", bayes = TRUE,
                      EV = EV, ...)
}

#' @rdname mice.impute.ztp
#' @export
mice.impute.ztnb.boot <- function(y, ry, x, wy = NULL, EV = FALSE, ...) {
  .countimp_1l_ztrunc(y, ry, x, wy = wy, dist = "ztnb", bayes = FALSE,
                      EV = EV, ...)
}
