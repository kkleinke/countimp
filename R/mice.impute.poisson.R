#' Multiple Imputation of Poisson Distributed Count Data
#' 
#' Imputes univariate missing data based  on a \code{poisson} GLM following either the Bayesian regression or bootstrap regression (appendix \code{.boot}) MI approach.
#' 
#' A Poisson GLM assumes that the mean of the count variable is equal to its variance (equidispersion assumption). For details, see Zeileis, Kleiber, & Jackman (2008), or Hilbe (2007).
#' The Bayesian method consists of the following steps:
#' \enumerate{
#' \item Fit the model, and find bhat, the posterior mean, and V(bhat), the posterior variance of model parameters b.
#' \item Draw b.star from N(bhat,V(bhat)).
#' \item Compute fitted values using \code{exp(x[!ry, ] \%*\% b.star)}
#' \item Simulate imputations from a Poisson distribution with mean parameter \code{lamda} being the respective fitted value from step 3.
#' }
#' The function uses the standard \code{glm.fit} function, using the \code{poisson} family. 
#' The bootstrap method draws a bootstrap sample from \code{y[ry]} and \code{x[ry,]} and consists of the following steps:
#' \enumerate{
#' \item Fit the model to the bootstrap sample and get model parameters \code{b.star}
#' \item Compute fitted values using \code{exp(x[!ry, ] \%*\% b.star)}
#' \item Simulate imputations from a Poisson distribution.
#' }
#' @param y Numeric vector with incomplete data
#' @param ry Response pattern of \code{y} (\code{TRUE}=observed, \code{FALSE}=missing)
#' @param x matrix with \code{length(y)} rows containing complete covariates
#' @param wy Logical vector of length \code{length(y)}. A \code{TRUE} value indicates locations in \code{y} for which imputations are created. Default is \code{!ry}
#' @param EV should automatic outlier handling of imputed values be enabled?  Default is \code{TRUE}: extreme imputations will be identified. These values will be replaced by imputations obtained by predictive mean matching (function \code{mice.impute.midastouch()})
#' @param ... Other named arguments.
#' @return Numeric vector of length \code{sum(!ry)} with imputations
#' @aliases mice.impute.poisson mice.impute.poisson.boot mice.impute.pois.boot mice.impute.pois
#' @references 
#' \itemize{
#' \item Hilbe, J. M. (2007). \emph{Negative binomial regression}. Cambridge: Cambridge University Press.
#' \item Kleinke, K., & Reinecke, J. (2013). \emph{countimp 1.0 -- A multiple imputation package for incomplete count data} [Technical Report]. University of Bielefeld, Faculty of Sociology, available from \url{www.uni-bielefeld.de/soz/kds/pdf/countimp.pdf}.
#' \item Rubin, D. B. (1987). \emph{Multiple imputation for nonresponse in surveys}. New York: Wiley.
#' \item Zeileis, A., Kleiber, C., & Jackman, S. (2008). Regression models for count data in R. \emph{Journal of Statistical Software}, 27(8), 1–-25.
#' }
#' @importFrom stats coef glm.fit poisson rpois summary.glm
#' @importFrom MASS rnegbin
#' @examples 
#' ## simulate Poisson distributed data
#' set.seed( 1234 )
#' b0 <- 1
#' b1 <- .75
#' b2 <- -.25
#' b3 <- .5
#' N <- 5000
#' x1 <- rnorm(N)
#' x2 <- rnorm(N)
#' x3 <- rnorm(N)
#' lam <- exp( b0 + b1 * x1 + b2 * x2 + b3 * x3 )
#' y <- rpois( N, lam )
#' POIS <- data.frame( y, x1, x2, x3 )
#'
#' ## introduce MAR missingness to simulated data
#' generate.md <- function( data, pos = 1, Z = 2, pmis = .5, strength = c( .5, .5 ) ) 
#' {
#'  total <- round( pmis * nrow(data) )
#'  sm <- which( data[,Z] < mean( data[,Z] ) )
#'  gr <- which( data[,Z] > mean( data[,Z] ) )
#'  sel.sm <- sample( sm, round( strength[1] * total ) )
#'  sel.gr <- sample( gr, round( strength[2] * total ) )
#'  sel <- c( sel.sm, sel.gr )
#'  data[sel,pos] <- NA
#'  return(data)
#' }
#' MPOIS <- generate.md( POIS, pmis = .2, strength = c( .2, .8) )
#'
#' ## impute missing data
#' imp <- countimp( MPOIS, method = c( "poisson" ,"" ,"" ,"" ))
#' @author Kristian Kleinke
#' @describeIn  mice.impute.poisson Bayesian regression variant
#' @export
mice.impute.poisson <-
function (y, ry, x, wy = NULL, EV=TRUE, ...) 
{
  if (is.null(wy)) 
    wy <- !ry
    x <- cbind(1, as.matrix(x))
    fit <- glm.fit(x[ry, ], y[ry], family = poisson(link = log))
    fit.sum <- summary.glm(fit)
    beta <- coef(fit)
    rv <- t(chol(fit.sum$cov.unscaled))
    beta.star <- beta + rv %*% rnorm(ncol(rv))
    p <- exp((x[wy, , drop = FALSE] %*% beta.star))
	im=rpois(length(p),p)
	imputed.values<-im
	if(EV){
	  outliers <- getOutliers(imputed.values, rho = c(0.3, 
	                                                  0.3), FLim = c(0.15, 0.85))
	  nans <- which(is.nan(imputed.values))
	  idx <- c(outliers$iLeft, outliers$iRight, nans)
	  if (length(idx) != 0) {
	    imputed.values[idx] <- NA
	    y[!ry] <- imputed.values
	    R = ry
	    ry <- !is.na(y)
	    new.values <- mice.impute.midastouch(y, ry, x, 
	                                         wy = NULL)
	    imputed.values[idx] <- new.values
	  }}
	return(imputed.values)
}

#' @export
#' @describeIn  mice.impute.poisson Bootstrap variant
mice.impute.poisson.boot <-
  function (y, ry, x, wy = NULL, EV=TRUE, ...) 
  {
    if (is.null(wy)) 
      wy <- !ry
    x <- cbind(1, as.matrix(x))
    xobs<-x[ry,]
    yobs<-y[ry]
    sel<-sample(1:length(yobs),length(yobs),replace=TRUE)
    xast<-xobs[sel,]
    yast<-yobs[sel]
    fit <- glm.fit(xast, yast, family = poisson(link = log))
    fit.sum <- summary.glm(fit)
    beta.star <- coef(fit)
    p <- exp((x[wy, , drop = FALSE] %*% beta.star))
    im=rpois(length(p),p)
    imputed.values<-im
    if(EV){
      outliers <- getOutliers(imputed.values, rho = c(0.3, 
                                                      0.3), FLim = c(0.15, 0.85))
      nans <- which(is.nan(imputed.values))
      idx <- c(outliers$iLeft, outliers$iRight, nans)
      if (length(idx) != 0) {
        imputed.values[idx] <- NA
        y[!ry] <- imputed.values
        R = ry
        ry <- !is.na(y)
        new.values <- mice.impute.midastouch(y, ry, x, 
                                             wy = NULL)
        imputed.values[idx] <- new.values
      }}
    return(imputed.values)
  }

#' @export
#' @describeIn  mice.impute.poisson Identical to \code{mice.impute.poisson}; included for backward compatibility
mice.impute.pois <-
  function (y, ry, x, wy = NULL, EV=TRUE, ...) 
  {
    if (is.null(wy)) 
      wy <- !ry
    x <- cbind(1, as.matrix(x))
    fit <- glm.fit(x[ry, ], y[ry], family = poisson(link = log))
    fit.sum <- summary.glm(fit)
    beta <- coef(fit)
    rv <- t(chol(fit.sum$cov.unscaled))
    beta.star <- beta + rv %*% rnorm(ncol(rv))
    p <- exp((x[wy, , drop = FALSE] %*% beta.star))
    im=rpois(length(p),p)
    imputed.values<-im
    if(EV){
      outliers <- getOutliers(imputed.values, rho = c(0.3, 
                                                      0.3), FLim = c(0.15, 0.85))
      nans <- which(is.nan(imputed.values))
      idx <- c(outliers$iLeft, outliers$iRight, nans)
      if (length(idx) != 0) {
        imputed.values[idx] <- NA
        y[!ry] <- imputed.values
        R = ry
        ry <- !is.na(y)
        new.values <- mice.impute.midastouch(y, ry, x, 
                                             wy = NULL)
        imputed.values[idx] <- new.values
      }}
    return(imputed.values)
  }

#' @export
#' @describeIn  mice.impute.poisson Identical to \code{mice.impute.poisson.boot}; included for backward compatibility
mice.impute.pois.boot <-
  function (y, ry, x, wy = NULL, EV=TRUE, ...) 
  {
    if (is.null(wy)) 
      wy <- !ry
    x <- cbind(1, as.matrix(x))
    xobs<-x[ry,]
    yobs<-y[ry]
    sel<-sample(1:length(yobs),length(yobs),replace=TRUE)
    xast<-xobs[sel,]
    yast<-yobs[sel]
    fit <- glm.fit(xast, yast, family = poisson(link = log))
    fit.sum <- summary.glm(fit)
    beta.star <- coef(fit)
    p <- exp((x[wy, , drop = FALSE] %*% beta.star))
    im=rpois(length(p),p)
    imputed.values<-im
    if(EV){
      outliers <- getOutliers(imputed.values, rho = c(0.3, 
                                                      0.3), FLim = c(0.15, 0.85))
      nans <- which(is.nan(imputed.values))
      idx <- c(outliers$iLeft, outliers$iRight, nans)
      if (length(idx) != 0) {
        imputed.values[idx] <- NA
        y[!ry] <- imputed.values
        R = ry
        ry <- !is.na(y)
        new.values <- mice.impute.midastouch(y, ry, x, 
                                             wy = NULL)
        imputed.values[idx] <- new.values
      }}
    return(imputed.values)
  }