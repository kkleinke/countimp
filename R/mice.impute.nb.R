#' Multiple Imputation of Overdispersed Count Data based on a Negative Binomial Model
#' 
#' Imputes univariate missing data based on a negative binomial model following either the Bayesian regression or bootstrap regression (appendix \code{.boot}) MI approach.
#' 
#' Overdispersed count data (meaning that the variance of the count variable is larger than its the mean) are typically analyzed by a negative binomial (NB) or by a \code{quasipoisson} model. The \code{quasipoisson} model is identical to an ordinary Poisson model, except that it estimates an additional dispersion parameter. For details, see Zeileis, Kleiber, & Jackman (2008), or Hilbe (2007).
#' The Bayesian method consists of the following steps:
#' \enumerate{
#' \item Fit the NB model, and find bhat, the posterior mean, and V(bhat), the posterior variance of model parameters b.
#' \item Draw b.star from N(bhat,V(bhat)).
#' \item Compute fitted values using \code{exp(x[!ry, ] \%*\% b.star)}
#' \item Simulate imputations from a negative binomial distribution to ensure an adequate dispersion of imputed values.
#' }
#' NB imputation uses function \code{glm.nb} from package \pkg{MASS} to fit the model. 
#' The bootstrap method draws a bootstrap sample from \code{y[ry]} and \code{x[ry,]} and consists of the following steps:
#' \enumerate{
#' \item Fit the model to the bootstrap sample and get model parameters \code{b.star}
#' \item Compute fitted values using \code{exp(x[!ry, ] \%*\% b.star)}
#' \item Simulate imputations from a negative binomial distribution to ensure an adequate dispersion of imputed values.
#' }
#' @param y Numeric vector with incomplete data
#' @param ry Response pattern of \code{y} (\code{TRUE}=observed, \code{FALSE}=missing)
#' @param x matrix with \code{length(y)} rows containing complete covariates
#' @param wy Logical vector of length \code{length(y)}. A \code{TRUE} value indicates locations in \code{y} for which imputations are created. Default is \code{!ry}
#' @param EV should automatic outlier handling of imputed values be enabled?  Default is \code{TRUE}: extreme imputations will be identified. These values will be replaced by imputations obtained by predictive mean matching (function \code{mice.impute.midastouch()})
#' @param ... Other named arguments.
#' @return Numeric vector of length \code{sum(!ry)} with imputations
#' @aliases mice.impute.nb.boot 
#' @references 
#' \itemize{
#' \item Hilbe, J. M. (2007). \emph{Negative binomial regression}. Cambridge: Cambridge University Press.
#' \item Kleinke, K., & Reinecke, J. (2013). \emph{countimp 1.0 -- A multiple imputation package for incomplete count data} [Technical Report]. University of Bielefeld, Faculty of Sociology, available from \url{www.uni-bielefeld.de/soz/kds/pdf/countimp.pdf}.
#' \item Rubin, D. B. (1987). \emph{Multiple imputation for nonresponse in surveys}. New York: Wiley.
#' \item Zeileis, A., Kleiber, C., & Jackman, S. (2008). Regression models for count data in R. \emph{Journal of Statistical Software}, 27(8), 1–-25.
#' }
#' @importFrom MASS rnegbin
#' @examples 
#' ## simulate overdespersed count data
#' set.seed( 1234 )
#' b0 <- 1
#' b1 <- .75
#' b2 <- -.25
#' b3 <- .5
#' N <- 5000
#' x1 <- rnorm(N)
#' x2 <- rnorm(N)
#' x3 <- rnorm(N)
#' mu <- exp( b0 + b1 * x1 + b2 * x2 + b3 * x3 )
#' y <- MASS::rnegbin( N, theta = 2, mu )
#' NB <- data.frame( y, x1, x2, x3 )
#' 
#' ## introduce MAR missingness to simulated data
#' total <- round( .2 * N )  ##number of missing data in y
#' sm <- which( NB[,2] < mean( NB[,2] ) )  ##subset: cases with x2<mean(x2)
#' gr <- which( NB[,2] > mean( NB[,2] ) )	##subset: cases with x2>mean(x2)
#' sel.sm <- sample( sm, round( .2 * total ) )	##select cases to set as missing
#' sel.gr <- sample( gr, round( .8 * total ) )	##select cases to set as missing
#' sel <- c( sel.sm,sel.gr )
#' MNB <- NB
#' MNB[sel,1] <- NA	##delete selected data
#' 
#' ## impute missing data
#' imp <- countimp( MNB, method = c( "nb.boot", "", "", "" )) 
#' @author Kristian Kleinke
#' @describeIn  mice.impute.nb Bayesian regression variant
#' @export
mice.impute.nb <-
function (y, ry, x, wy = NULL, EV=TRUE, ...) 
{
  if (is.null(wy)) 
    wy <- !ry
  Y=y[ry]
	X=as.matrix(x[ry,])

	nam=paste("V",1:ncol(X),sep="")
	colnames(X)=nam

	form=paste("Y","~",paste(nam,collapse="+"))
	form=as.formula(form)
	dat=data.frame(Y,X)

	fit <- MASS::glm.nb(form , data=dat)

	fit.sum <- summary(fit)
	beta <- coef(fit)
  rv <- t(chol(fit.sum$cov.unscaled))
	b.star <- beta + rv %*% rnorm(ncol(rv))
  p <- exp((data.matrix(cbind(1,x[wy, , drop = FALSE]))%*% b.star))
	im=MASS::rnegbin(n=length(p),theta=fit.sum$theta,mu=p)
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

#' @describeIn  mice.impute.nb Bootstrap regression variant
#' @export
mice.impute.nb.boot <-
  function (y, ry, x, wy = NULL, EV=TRUE, ...) 
  {
    if (is.null(wy)) 
      wy <- !ry
    Y=y[ry]
    X=as.matrix(x[ry,])
    
    nam=paste("V",1:ncol(X),sep="")
    colnames(X)=nam
    
    form=paste("Y","~",paste(nam,collapse="+"))
    form=as.formula(form)
    dat=data.frame(Y,X)
    
    
    sel<-sample(1:length(Y),length(Y),replace=TRUE)
    datstar<-dat[sel,]
    
    fit <- MASS::glm.nb(form , data=datstar)
    
    fit.sum <- summary(fit)
    b.star <- coef(fit)
    p <- exp((data.matrix(cbind(1,x[wy, , drop = FALSE]))%*% b.star))
    im=MASS::rnegbin(n=length(p),theta=fit.sum$theta,mu=p)
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