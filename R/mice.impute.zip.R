#' Multiple Imputation of Flat File Zero-Inflated Count Data
#'
#' The functions impute flat file zero-inflated count data based on a Poisson or negative binomial hurdle model, either using a Bayesian regression or a bootstrap regression approach (appendix: ``\code{.boot}''). Alternatively, a zero-inflated Poisson or NB model can be specified. Hurdle models are mixture models and consist of two model components: the zero model (a binomial GLM), determining, if the observational unit has a zero or non-zero value, and the count model (a zero-truncated Poisson or NB model), determining, what non-zero value the observational unit has.
#' Zero-inflation models are also mixture models and specify a zero model (here a logit model, determining if the observational unit has a ``certain zero'' or not) and a count model (here a Poisson or negative binomial model), determining, what count - both zero and non-zero - the observational unit has. Different sets of covariates (predictors) may be used for the zero and the count models.
#'
#'The functions multiply impute incomplete zero-inflated count data using either the \code{zeroinfl()} function (zero-inflation model) or the \code{hurdle()} function (hurdle model) from package \pkg{pscl} (Zeileis, Kleiber, & Jackman, 2008). 
#'  Model specification details:
#'    \itemize{
#'      \item 0 = variable not included in imputation model
#'      \item 1 = variable will be included in the zero \emph{and} the count model
#'      \item 2 = variable will be included in the count model
#'      \item 3 = variable will be included in the zero model
#'    }
#' The Bayesian regression variants (see Rubin 1987, p. 169-170) consist of the following steps:
#'  \enumerate{
#'  \item Fit the model; find bhat, the posterior mean, and V(bhat), the posterior variance of model parameters b
#'  \item Draw b* from N(bhat,V(bhat))
#'  \item Compute predicted probabilities for observing each count \code{p}
#'  \item Draw imputations from observed counts with selection probabilities \code{p}
#'  }
#' The bootstrap functions draw a bootstrap sample from \code{y[ry]} and \code{x[ry,]} 
#' \enumerate{
#' \item Fit the model to the bootstrap sample
#'  \item Compute predicted probabilities for observing each count \code{p}
#'  \item Draw imputations from observed counts with selection probabilities \code{p}
#' }
#' @param y Numeric vector with incomplete data
#' @param ry Response pattern of \code{y} (\code{TRUE}=observed, \code{FALSE}=missing)
#' @param x matrix with \code{length(y)} rows containing complete covariates
#' @param type vector of length \code{ncol(x)} determining the imputation model; \code{type} is automatically extracted from the \code{predictorMatrix} argument of \code{mice()}.
#' @param wy Logical vector of length \code{length(y)}. A \code{TRUE} value indicates locations in \code{y} for which imputations are created. Default is \code{!ry}
#' @return vector with imputations
#' @aliases mice.impute.hnb mice.impute.hnb.boot mice.impute.hp mice.impute.hp.boot 
#' @aliases mice.impute.zinb mice.impute.zinb.boot mice.impute.zip mice.impute.zip.boot
#' @author Kristian Kleinke
#' @references
#' \itemize{
#' \item Kleinke, K., & Reinecke, J. (2013a). Multiple Imputation of incomplete zero-inflated count data. \emph{Statistica Neerlandica}, available from \url{http://onlinelibrary.wiley.com/doi/10.1111/stan.12009/abstract}.
#' \item Kleinke, K., & Reinecke, J. (2013b). \emph{countimp 1.0 -- A multiple imputation package for incomplete count data} [Technical Report]. University of Bielefeld, Faculty of Sociology, available from \url{www.uni-bielefeld.de/soz/kds/pdf/countimp.pdf}.
#' \item Rubin, D. B. (1987). \emph{Multiple imputation for nonresponse in surveys}. New York: Wiley.
#' \item Zeileis, A., Kleiber, C., & Jackman, S. (2008). Regression models for count data in R. \emph{Journal of Statistical Software}, 27(8), 1--25.
#'} 
#' @examples 
#' ## Example 1:
#' data(crim4w)
#' ini <- countimp(crim4w, maxit=0)
#' meth <- ini$method
#' meth[6:7] <- "hp"
#' meth[8:9] <- "pmm"
#' pred <- ini$predictorMatrix
#' pred[,"id"] <- 0
#' pred["ACRIM",] <- c(0,1,3,2,0,3,3,2,1)
#' imp <- countimp( data = crim4w, method = meth, predictorMatrix = pred )
#'
#' ## Example 2:
#' ## Simulate zero-inflated NB data
#' b0 <- 1
#' b1 <- .3
#' b2 <- .3
#' c0 <- 0
#' c1 <- 2
#' theta <- 1
#' require("pscl")
#' set.seed(1234)
#' N <- 10000
#' x1 <- rnorm(N)
#' x2 <- rnorm(N)
#' x3 <- rnorm(N)
#' mu <- exp( b0 + b1 * x1 + b2 * x2 )
#' yzinb <- rnegbin( N, mu, theta)
#' pzero <- plogis( c1 * x3 )        # zero-infl. prob. depends on x3
#' ## Introduce zero-inflation
#' uni <- runif(N)
#' yzinb[uni < pzero] <- 0
#' zinbdata<-data.frame(yzinb,x1,x2,x3)
#'
#' ## Generate MAR missingness
#' generate.md <- function( data, pos = 1, Z = 2, pmis = .5, strength = c( .5, .5 ) )
#' {
#' total <- round( pmis * nrow(data) )
#'  sm <- which( data[,Z] < mean( data[,Z] ) )
#'  gr <- which( data[,Z] > mean( data[,Z] ) )
#'  sel.sm <- sample( sm, round( strength[1] * total ) )
#'  sel.gr <- sample( gr, round( strength[2] * total ) )
#'  sel <- c( sel.sm, sel.gr )
#'  data[sel,pos] <- NA
#'  return(data)
#' }
#' zinbmdata <- generate.md( zinbdata, pmis = .3, strength = c( .2, .8) )
#'
#' ## Impute missing data
#' ini <- mice( zinbmdata, m = 5, maxit = 0)
#' pred <- ini$predictorMatrix 
#' pred[1,] <- c(0, 2, 2, 3)
#' meth<-ini$method
#' meth[1] <- "zinb"
#' imp.zinb <- countimp( zinbmdata, m = 5, method = meth,
#'             predictorMatrix = pred, seed = 1234, print = FALSE)
#' @import pscl             
#' @importFrom stats C complete.cases contr.treatment cor formula model.frame model.matrix pt qt
#' @export
#' @describeIn mice.impute.2l.zip zero-inflated Poisson model; Bayesian regression variant
mice.impute.zip <-
function(y, ry, x, type, wy = NULL){
  if (is.null(wy)) 
    wy <- !ry
  Y <- y[ry]
  X <- x[ry,]
  X <- data.frame(X)
  nam <- colnames(X)
  b <- which(type==1) # variables used in zero AND count model
  c <- which(type==2) # count model ONLY variables
  z <- which(type==3) # zero model ONLY variables
  zero <- c(b,z); zero <- unique(zero); zero <- sort(zero)
  count <- c(b,c); count <- unique(count); count <- sort(count)
  if (is.integer(zero) && length(zero) == 0L){zeroform="| 1"}else{zeroform=paste("|",paste(nam[zero],collapse="+"))}
  if (is.integer(count) && length(count) == 0L){countform="Y~ 1"}else{countform=paste("Y","~ ",paste(nam[count],collapse="+"))}
  form <-
    as.formula(paste(countform, zeroform))    
  dat <- data.frame(Y,X)
  fit <- pscl::zeroinfl(form,data=dat,dist="poisson",link="logit")
  fit.sum <- summary(fit)
  beta <- coef(fit)
  rv <- t(chol(fit.sum$vcov))
  b.star <- beta+rv %*% rnorm(ncol(rv))
  fit$coefficients$count <-
  b.star[1:length(fit$coefficients$count)]
  fit$coefficients$zero <-
  b.star[(length(fit$coefficients$count)+1):length(b.star)]
  newdata <- data.frame(X=x[wy, , drop = FALSE])
  colnames(newdata) <- nam
  pc <- predict(fit,
  newdata=newdata,type="prob",na.action=na.pass)
  pcvec <- 1:nrow(pc)
  for (i in 1:nrow(pc))
  {
    pcvec[i] <-
    sample(as.numeric(names(pc[i,])),1,pc[i,] ,
    replace=TRUE)
  }
  return(pcvec)
}

#' @export
#' @describeIn mice.impute.2l.zip zero-inflated Poisson model; Bootstrap regression variant
mice.impute.zip.boot <-
  function(y, ry, x, type, wy = NULL){
    if (is.null(wy)) 
      wy <- !ry
    Y <- y[ry]
    X <- x[ry,]
    X <- data.frame(X)
    nam <- colnames(X)
    b <- which(type==1) # variables used in zero AND count model
    c <- which(type==2) # count model ONLY variables
    z <- which(type==3) # zero model ONLY variables
    zero <- c(b,z); zero <- unique(zero); zero <- sort(zero)
    count <- c(b,c); count <- unique(count); count <- sort(count)
    if (is.integer(zero) && length(zero) == 0L){zeroform="| 1"}else{zeroform=paste("|",paste(nam[zero],collapse="+"))}
    if (is.integer(count) && length(count) == 0L){countform="Y~ 1"}else{countform=paste("Y","~ ",paste(nam[count],collapse="+"))}
    form <-
      as.formula(paste(countform, zeroform))    
    dat <- data.frame(Y,X)
    datBS <- dat[sample(1:length(Y),length(Y),replace=TRUE),]
    fit <- pscl::zeroinfl(form,
                    data=datBS,dist="poisson",link="logit")
    newdata <- data.frame(X=x[wy, , drop = FALSE])
    colnames(newdata) <- nam
    pc <- predict(fit,
                  newdata=newdata,type="prob",na.action=na.pass)
    pcvec <- 1:nrow(pc)
    for (i in 1:nrow(pc))
    {
      pcvec[i] <-
        sample(as.numeric(names(pc[i,])),1, pc[i,],
               replace=TRUE)
    }
    return(pcvec)
    }

#' @export
#' @describeIn mice.impute.2l.zip zero-inflated NB model; Bayesian regression variant




mice.impute.zinb <-
  function(y, ry, x, type, wy = NULL){
    if (is.null(wy)) 
      wy <- !ry
    Y <- y[ry]
    X <- x[ry,]
    X <- data.frame(X)
    nam <- colnames(X)
    b <- which(type==1) # variables used in zero AND count model
    c <- which(type==2) # count model ONLY variables
    z <- which(type==3) # zero model ONLY variables
    zero <- c(b,z); zero <- unique(zero); zero <- sort(zero)
    count <- c(b,c); count <- unique(count); count <- sort(count)
    if (is.integer(zero) && length(zero) == 0L){zeroform="| 1"}else{zeroform=paste("|",paste(nam[zero],collapse="+"))}
    if (is.integer(count) && length(count) == 0L){countform="Y~ 1"}else{countform=paste("Y","~ ",paste(nam[count],collapse="+"))}
    form <-
      as.formula(paste(countform, zeroform))
    dat <- data.frame(Y,X)
    fit <- pscl::zeroinfl(form,data=dat,dist="negbin",link="logit")
    fit.sum <- summary(fit)
    beta <- coef(fit)
    rv <- t(chol(fit.sum$vcov))
    b.star <- beta+rv %*% rnorm(ncol(rv))
    fit$coefficients$count <-
      b.star[1:length(fit$coefficients$count)]
    fit$coefficients$zero <-
      b.star[(length(fit$coefficients$count)+1):length(b.star)]
    newdata <- data.frame(X=x[wy, , drop = FALSE])
    colnames(newdata) <- nam
    pc <- predict(fit,
                  newdata=newdata,type="prob",na.action=na.pass)
    pcvec <- 1:nrow(pc)
    for (i in 1:nrow(pc))
    {
      pcvec[i] <-
        sample(as.numeric(names(pc[i,])),1,pc[i,] ,
               replace=TRUE)
    }
    return(pcvec)
    }

#' @export
#' @describeIn mice.impute.2l.zip zero-inflated NB model; Bootstrap regression variant
mice.impute.zinb.boot <-
  function(y, ry, x, type, wy = NULL){
    if (is.null(wy)) 
      wy <- !ry
    Y <- y[ry]
    X <- x[ry,]
    X <- data.frame(X)
    nam <- colnames(X)
    b <- which(type==1) # variables used in zero AND count model
    c <- which(type==2) # count model ONLY variables
    z <- which(type==3) # zero model ONLY variables
    zero <- c(b,z); zero <- unique(zero); zero <- sort(zero)
    count <- c(b,c); count <- unique(count); count <- sort(count)
    if (is.integer(zero) && length(zero) == 0L){zeroform="| 1"}else{zeroform=paste("|",paste(nam[zero],collapse="+"))}
    if (is.integer(count) && length(count) == 0L){countform="Y~ 1"}else{countform=paste("Y","~ ",paste(nam[count],collapse="+"))}
    form <-
      as.formula(paste(countform, zeroform))    
    dat <- data.frame(Y,X)
    datBS <- dat[sample(1:length(Y),length(Y),replace=TRUE),]
    fit <- pscl::zeroinfl(form,
                    data=datBS,dist="negbin",link="logit")
    newdata <- data.frame(X=x[wy, , drop = FALSE])
    colnames(newdata) <- nam
    pc <- predict(fit,
                  newdata=newdata,type="prob",na.action=na.pass)
    pcvec <- 1:nrow(pc)
    for (i in 1:nrow(pc))
    {
      pcvec[i] <-
        sample(as.numeric(names(pc[i,])),1, pc[i,],
               replace=TRUE)
    }
    return(pcvec)
    }

#' @export
#' @describeIn mice.impute.2l.zip hurdle Poisson model; Bayesian regression variant
mice.impute.hp <-
  function (y, ry, x, type, wy = NULL)
  {
    if (is.null(wy)) 
      wy <- !ry
    Y <- y[ry]
    X <- x[ry, ]
    X <- data.frame(X)
    nam <- colnames(X)
    b <- which(type == 1)
    c <- which(type == 2)
    z <- which(type == 3)
    zero <- c(b, z)
    zero <- unique(zero)
    zero <- sort(zero)
    count <- c(b, c)
    count <- unique(count)
    count <- sort(count)
    if (is.integer(zero) && length(zero) == 0L){zeroform="| 1"}else{zeroform=paste("|",paste(nam[zero],collapse="+"))}
    if (is.integer(count) && length(count) == 0L){countform="Y~ 1"}else{countform=paste("Y","~ ",paste(nam[count],collapse="+"))}
    form <-
      as.formula(paste(countform, zeroform))    
    dat <- data.frame(Y, X)
    fit <- pscl::hurdle(form, data = dat, dist = "poisson")
    fit.sum <- summary(fit)
    beta <- coef(fit)
    rv <- t(chol(fit.sum$vcov))
    b.star <- beta + rv %*% rnorm(ncol(rv))
    fit$coefficients$count <- b.star[1:length(fit$coefficients$count)]
    fit$coefficients$zero <- b.star[(length(fit$coefficients$count) +
                                       1):length(b.star)]
    newdata <- data.frame(X = x[wy, , drop = FALSE])
    colnames(newdata) <- nam
    pc <- predict(fit, newdata = newdata, type = "prob", na.action = na.pass)
    pcvec <- 1:nrow(pc)
    for (i in 1:nrow(pc)) {
      pcvec[i] <- sample(as.numeric(names(pc[i, ])), 1, pc[i,
                                                           ], replace = TRUE)
    }
    return(pcvec)
  }

#' @export
#' @describeIn mice.impute.2l.zip hurdle Poisson model; Bootstrap regression variant
mice.impute.hp.boot <-
  function (y, ry, x, type, wy = NULL)
  {
    if (is.null(wy)) 
      wy <- !ry
    Y <- y[ry]
    X <- x[ry, ]
    X <- data.frame(X)
    nam <- colnames(X)
    b <- which(type == 1)
    c <- which(type == 2)
    z <- which(type == 3)
    zero <- c(b, z)
    zero <- unique(zero)
    zero <- sort(zero)
    count <- c(b, c)
    count <- unique(count)
    count <- sort(count)
    if (is.integer(zero) && length(zero) == 0L){zeroform="| 1"}else{zeroform=paste("|",paste(nam[zero],collapse="+"))}
    if (is.integer(count) && length(count) == 0L){countform="Y~ 1"}else{countform=paste("Y","~ ",paste(nam[count],collapse="+"))}
    form <-
      as.formula(paste(countform, zeroform))    
    dat <- data.frame(Y, X)
    datBS <- dat[sample(1:length(Y), length(Y), replace = TRUE),
                 ]
    fit <- pscl::hurdle(form, data = datBS, dist = "poisson")
    newdata <- data.frame(X = x[wy, , drop = FALSE])
    colnames(newdata) <- nam
    pc <- predict(fit, newdata = newdata, type = "prob", na.action = na.pass)
    pcvec <- 1:nrow(pc)
    for (i in 1:nrow(pc)) {
      pcvec[i] <- sample(as.numeric(names(pc[i, ])), 1, pc[i,
                                                           ], replace = TRUE)
    }
    return(pcvec)
  }

#' @export
#' @describeIn mice.impute.2l.zip hurdle NB model; Bayesian regression variant
mice.impute.hnb <-
  function (y, ry, x, type, wy = NULL)
  {
    if (is.null(wy)) 
      wy <- !ry
    Y <- y[ry]
    X <- x[ry, ]
    X <- data.frame(X)
    nam <- colnames(X)
    b <- which(type == 1)
    c <- which(type == 2)
    z <- which(type == 3)
    zero <- c(b, z)
    zero <- unique(zero)
    zero <- sort(zero)
    count <- c(b, c)
    count <- unique(count)
    count <- sort(count)
    if (is.integer(zero) && length(zero) == 0L){zeroform="| 1"}else{zeroform=paste("|",paste(nam[zero],collapse="+"))}
    if (is.integer(count) && length(count) == 0L){countform="Y~ 1"}else{countform=paste("Y","~ ",paste(nam[count],collapse="+"))}
    form <-
      as.formula(paste(countform, zeroform))    
    dat <- data.frame(Y, X)
    fit <- pscl::hurdle(form, data = dat, dist = "negbin")
    fit.sum <- summary(fit)
    beta <- coef(fit)
    rv <- t(chol(fit.sum$vcov))
    b.star <- beta + rv %*% rnorm(ncol(rv))
    fit$coefficients$count <- b.star[1:length(fit$coefficients$count)]
    fit$coefficients$zero <- b.star[(length(fit$coefficients$count) +
                                       1):length(b.star)]
    newdata <- data.frame(X = x[wy, , drop = FALSE])
    colnames(newdata) <- nam
    pc <- predict(fit, newdata = newdata, type = "prob", na.action = na.pass)
    pcvec <- 1:nrow(pc)
    for (i in 1:nrow(pc)) {
      pcvec[i] <- sample(as.numeric(names(pc[i, ])), 1, pc[i,
                                                           ], replace = TRUE)
    }
    return(pcvec)
  }

#' @export
#' @describeIn mice.impute.2l.zip hurdle NB model; Bootstrap regression variant
mice.impute.hnb.boot <-
  function (y, ry, x, type, wy=NULL)
  {
    if (is.null(wy)) 
      wy <- !ry
    Y <- y[ry]
    X <- x[ry, ]
    X <- data.frame(X)
    nam <- colnames(X)
    b <- which(type == 1)
    c <- which(type == 2)
    z <- which(type == 3)
    zero <- c(b, z)
    zero <- unique(zero)
    zero <- sort(zero)
    count <- c(b, c)
    count <- unique(count)
    count <- sort(count)
    if (is.integer(zero) && length(zero) == 0L){zeroform="| 1"}else{zeroform=paste("|",paste(nam[zero],collapse="+"))}
    if (is.integer(count) && length(count) == 0L){countform="Y~ 1"}else{countform=paste("Y","~ ",paste(nam[count],collapse="+"))}
    form <-
      as.formula(paste(countform, zeroform))    
    dat <- data.frame(Y, X)
    datBS <- dat[sample(1:length(Y), length(Y), replace = TRUE),
                 ]
    fit <- pscl::hurdle(form, data = datBS, dist = "negbin")
    
    newdata <- data.frame(X = x[wy, , drop = FALSE])
    colnames(newdata) <- nam
    pc <- predict(fit, newdata = newdata, type = "prob", na.action = na.pass)
    pcvec <- 1:nrow(pc)
    for (i in 1:nrow(pc)) {
      pcvec[i] <- sample(as.numeric(names(pc[i, ])), 1, pc[i,
                                                           ], replace = TRUE)
    }
    return(pcvec)
  }
