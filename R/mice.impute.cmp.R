#' Imputation of underdispersed counts by Conway-Maxwell-Poisson regression
#'
#' Imputes missing counts whose variance is SMALLER than their mean. The
#' Poisson model cannot represent that case at all -- it fixes the
#' variance-to-mean ratio at one -- and the negative binomial can only move the
#' ratio upwards, so both produce imputations that are too dispersed. The
#' Conway-Maxwell-Poisson distribution covers under- and overdispersion in one
#' family with a single shape parameter \code{nu}.
#'
#' The regression is in the mean parametrisation: the linear predictor gives
#' \eqn{E[Y]} directly, and the rate parameter is solved per case from the mean
#' and the shape. \code{nu > 1} is underdispersion, \code{nu = 1} is exactly
#' Poisson, \code{nu < 1} is overdispersion. The shape is estimated, not fixed,
#' so a variable that turns out to be equidispersed simply returns
#' \code{nu} near one.
#'
#' \code{mice.impute.cmp} draws the coefficients \emph{and} the shape jointly
#' from their asymptotic posterior; \code{mice.impute.cmp.boot} bootstraps the
#' observed cases and uses the resulting estimates directly. Use
#' \code{countimp_fit_diag} on the observed data to see whether the
#' variance-to-mean ratio warrants this family.
#'
#' @param y Numeric vector of length \code{n}, the incomplete count variable.
#' @param ry Logical vector of length \code{n}, \code{TRUE} where \code{y} is
#'   observed.
#' @param x Numeric matrix or data frame of \code{n} rows with the predictors.
#'   An intercept is added internally.
#' @param wy Logical vector of length \code{n} marking the cases to impute.
#'   Defaults to \code{!ry}.
#' @param EV Logical. If \code{TRUE}, screen the imputations for extreme values
#'   as described in Vink and Van Buuren; see \code{mice.impute.pois}.
#' @param ... Further arguments, ignored.
#'
#' @return Numeric vector of length \code{sum(wy)} with the imputed counts.
#'
#' @section Underdispersion in practice:
#' Counts bounded by design are the common case: items answered correctly out of
#' a fixed number, symptoms endorsed from a checklist, days attended out of a
#' school week. Where the bound is known, \code{\link{mice.impute.bp}} imposes
#' it exactly. Where it is not -- or where the underdispersion comes from
#' something other than a ceiling -- this family models the dispersion rather
#' than the support.
#'
#' @section Cost:
#' The normalising constant has no closed form, so every likelihood evaluation
#' sums over a grid. The grid is a window around the mode whose length follows
#' the standard deviation, which for underdispersed data is narrower than
#' Poisson's; the fit is nevertheless several times slower than
#' \code{mice.impute.pois}. Measured against \code{glmmTMB} with family
#' \code{compois} on the same data, the coefficients, the shape and the
#' log-likelihood agree to four decimals.
#'
#' @references
#' Shmueli, G., Minka, T. P., Kadane, J. B., Borle, S., & Boatwright, P. (2005).
#' A useful distribution for fitting discrete data: revival of the
#' Conway-Maxwell-Poisson distribution. \emph{Journal of the Royal Statistical
#' Society C}, 54(1), 127-142.
#'
#' Huang, A. (2017). Mean-parametrized Conway-Maxwell-Poisson regression models
#' for dispersed counts. \emph{Statistical Modelling}, 17(6), 359-370.
#'
#' Kleinke, K., & Reinecke, J. (2013). Multiple imputation of incomplete
#' zero-inflated count data. \emph{Statistica Neerlandica}, 67(3), 311-336.
#'
#' @seealso \code{\link{mice.impute.bp}} for counts with known bounds,
#'   \code{\link{mice.impute.nb}} for overdispersion,
#'   \code{\link{countimp_fit_diag}} for choosing between them.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 300
#' x <- rnorm(n)
#' ## underdispersed counts: binomial has variance below its mean
#' y <- rbinom(n, 12, plogis(-0.2 + 0.4 * x))
#' d <- data.frame(y = y, x = x)
#' d$y[runif(n) < 0.25] <- NA
#' imp <- countimp(d, method = c(y = "cmp", x = ""), m = 2, printFlag = FALSE)
#' }
#'
#' @rdname mice.impute.cmp
#' @export
mice.impute.cmp <- function(y, ry, x, wy = NULL, EV = FALSE, ...) {
  .countimp_1l_cmp(y, ry, x, wy = wy, bayes = TRUE, EV = EV, ...)
}

#' @rdname mice.impute.cmp
#' @export
mice.impute.cmp.boot <- function(y, ry, x, wy = NULL, EV = FALSE, ...) {
  .countimp_1l_cmp(y, ry, x, wy = wy, bayes = FALSE, EV = EV, ...)
}
