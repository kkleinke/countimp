#' Multiple Imputation Inference
#' 
#' Multiple Imputation Inference using Rubin's (1987) rules for MI inference.
#' 
#' Combines \code{m} sets of parmater estimates and standard errors into an overall set of estimate using Rubin's (1987) rules for multiple imputation inference. This function is an adaptation of function \code{mi.inference} from package \pkg{norm}.
#'
#' @param est a list of \code{m} vectors containing parameter estimates (e.g. regression coefficients)
#' @param std.err a list of \code{m} vectors with corresponding standard errors
#' @param confidence confidence level; default is to compute 95\% confidence intervals
#' @return a list of vectors
#' \describe{
#' \item{est}{combined parameter estimates}
#' \item{std.err}{combined standard errors representing both between and  within-imputation variability}
#' \item{t.value}{t-ratio: \code{t.value}=\code{est}/\code{std.err}}
#' \item{df}{degrees of freedom}	
#' \item{signif}{p-values for the two-tailed tests that the respective combined estimate is zero}
#' \item{lower}{lower limits of the MI confidence intervals}
#' \item{upper}{upper limits of the MI confidence intervals}
#' \item{r}{relative increases in variance due to nonresponse}	
#' \item{fminf}{fractions of missing information} 
#' }
#' @references 
#' \itemize{
#' \item Rubin, D. B. (1987). \emph{Multiple imputation for nonresponse in surveys}. New York: Wiley.
#' \item Schafer, J. L. (1997). \emph{Analysis of incomplete multivariate data}. London: Chapman & Hall.
#' }
#' @aliases miinference
#' @examples 
#' ## simulate overdespersed count data
#' b0 <- 1
#' b1 <- .75
#' b2 <- -.25
#' b3 <- .5
#' N <- 5000
#' x1 <- rnorm(N)
#' x2 <- rnorm(N)
#' x3 <- rnorm(N)
#' mu <- exp(b0+b1*x1+b2*x2+b3*x3)
#' y <- MASS::rnegbin( N, theta = 2, mu)
#' NB <- data.frame(y,x1,x2,x3)
#'
#' ## introduce MAR missingness to simulated data
#' total <- round(.2 * N)  ##number of missing data in y     
#' sm <- which( NB[,2] < mean(NB[,2]) )  ##subset: cases with x2<mean(x2)
#' gr <- which( NB[,2] > mean(NB[,2]) )	##subset: cases with x2>mean(x2)
#' sel.sm <- sample( sm, round(.2*total) )	##select cases to set as missing
#' sel.gr <- sample( gr, round(.8*total) )	##select cases to set as missing
#' sel <- c(sel.sm, sel.gr)
#' MNB <- NB
#' MNB[sel,1] <- NA	##delete selected data
#'
#' ## imputation and repeated data analysis
#' imp <- mice( MNB, method=c("nb","","","") )
#' res <- with( imp, glm.nb(y~x1+x2+x3) )

#' ## get MI inferences for dispersion parameter theta
#' EST <- vector( length = 5, mode = "list" )
#' SE <- vector( length = 5, mode = "list" )
#'
#' for (i in 1:5){
#' EST[[i]] <- res$analyses[[i]]$theta 
#' SE[[i]] <- res$analyses[[i]]$SE.theta 
#' }
#'
#' miinference(EST,SE)
#' @export
miinference<-function (est, std.err, confidence = 0.95) 
{
  qstar <- est[[1]]
  for (i in 2:length(est)) {
    qstar <- cbind(qstar, est[[i]])
  }
  qbar <- apply(qstar, 1, mean)
  u <- std.err[[1]]
  for (i in 2:length(std.err)) {
    u <- cbind(u, std.err[[i]])
  }
  names(u) <- names(qstar)
  u <- u^2
  ubar <- apply(u, 1, mean)
  bm <- apply(qstar, 1, var)
  m <- dim(qstar)[2]
  tm <- ubar + ((1 + (1/m)) * bm)
  rem <- (1 + (1/m)) * bm/ubar
  nu <- (m - 1) * (1 + (1/rem))^2
  alpha <- 1 - (1 - confidence)/2
  low <- qbar - qt(alpha, nu) * sqrt(tm)
  up <- qbar + qt(alpha, nu) * sqrt(tm)
  pval <- 2 * (1 - pt(abs(qbar/sqrt(tm)), nu))
  fminf <- (rem + 2/(nu + 3))/(rem + 1)
  result <- list(est = qbar, std.err = sqrt(tm),t.value = qbar/sqrt(tm) ,df = nu, p.value = pval, 
                 lower = low, upper = up, r = rem, fminf = fminf)
  result
}