#' Multiple Imputation of Incomplete Two-Level Count Data
#'
#' The functions impute multilevel count data based on a two-level Poisson or negative binomial model, either using a Bayesian regression or a bootstrap regression approach (appendix: ``\code{.boot}''). The \code{.noint} variants treat the intercept only as a fixed, but \emph{not} as a random effect. Package \pkg{glmmTMB} is used to fit the model. 
#' 
#'
#' Model specification details:
#'    \itemize{
#'      \item -2 = class variable (only one class variable is allowed!)
#'      \item 0 = variable not included in imputation model
#'      \item 1 = variable will be included as a fixed effect 
#'      \item 2 = variable will be included as a fixed \emph{and} random effect 
#'    }
#' The Bayesian regression variants (see Rubin 1987, p. 169-170) consist of the following steps:
#'  \enumerate{
#'  \item Fit the model; find bhat, the posterior mean, and V(bhat), the posterior variance of model parameters b
#'  \item Draw b* from N(bhat,V(bhat))
#'  \item Obtain fitted values based on b*
#'  \item Draw imputations for the incomplete part from appropriate distribution (Poisson or NB)
#'  }
#' The bootstrap functions draw a bootstrap sample from \code{y[ry]} and \code{x[ry,]} (Note: we resample clusters rather than individual cases) and consist of the following steps:
#' \enumerate{
#' \item Fit the model to the bootstrap sample
#'  \item Obtain fitted values
#'  \item Draw imputations for the incomplete part from appropriate distribution (Poisson or NB)
#' }
#' @param y Numeric vector with incomplete data in long format (i.e. the groups are stacked upon each other)
#' @param ry Response pattern of \code{y} (\code{TRUE}=observed, \code{FALSE}=missing)
#' @param x matrix with \code{length(y)} rows containing complete covariates; also in long format
#' @param type vector of length \code{ncol(x)} identifying fixed, random, and class variables; \code{type} is automatically extracted from the \code{predictorMatrix}; see \pkg{mice}'s user's manual for details about how to specify the imputation model; see also section ``details''.
#' @param intercept \code{TRUE}: model will include intercept as a random effect; \code{FALSE}: intercept will be treated as fixed.
#' @param wy Logical vector of length \code{length(y)}. A \code{TRUE} value indicates locations in \code{y} for which imputations are created. Default is \code{!ry}
#' @param EV should automatic outlier handling of imputed values be enabled?  Default is \code{TRUE}: extreme imputations will be identified. These values will be replaced by imputations obtained by predictive mean matching (function \code{mice.impute.midastouch()})
#' @return Numeric vector of length \code{sum(!ry)} with imputations
#' @export
#' @import glmmTMB extremevalues
#' @aliases mice.impute.2l.poisson mice.impute.2l.poisson.boot mice.impute.2l.poisson.noint mice.impute.2l.poisson.noint.boot   
#' @aliases mice.impute.2l.nb mice.impute.2l.nb.boot mice.impute.2l.nb.noint mice.impute.2l.nb.noint.boot   
#' @aliases mice.impute.2l.nb2 mice.impute.2l.nb2.boot mice.impute.2l.nb2.noint mice.impute.2l.nb2.noint.boot   
#' @author Kristian Kleinke
#' @describeIn mice.impute.2l.poisson Bayesian Poisson regression variant; random intercept
#' @export
mice.impute.2l.poisson<-
  function (y, ry, x, type, intercept = TRUE, wy = NULL, EV=TRUE)
  {
    if (is.null(wy)) 
      wy <- !ry
    Y = y[ry]
    X = x[ry, ]
    if (sum(type == -2) > 1) {
      stop("only one class allowed!")
    }
    grp = which(type == -2)
    nam = colnames(X)
    fixedeff = paste(paste(nam[-grp], collapse = "+"), sep = "")
    if (any(type == 2)) {
      ran = which(type == 2)
      randeff = paste("+", paste(nam[ran], collapse = "+"),
                      sep = "")
      if (!intercept) {
        randeff = paste("(0", randeff, "|", paste(nam[grp],")"),
                        sep = "")
      }
      else {
        randeff = paste("(1", randeff, "|", paste(nam[grp], ")"),
                        sep = "")
      }
    }
    else {
      if (!intercept) {
        randeff <- paste("(0", "|", paste(nam[grp]),")", sep = "")
      }
      else {
        randeff <- paste("(1", "|", paste(nam[grp]),")", sep = "")
      }
    }
    
    form = as.formula(paste("Y", "~", fixedeff, "+", randeff, sep = ""))
    dat = data.frame(Y, X)
    fit <- glmmTMB(formula = form,  data = dat,
                   family=list(family="poisson",link="log"))
    
    fit.sum <- summary(fit)
    beta <- fixef(fit)$cond
    rv <- t(chol(vcov(fit)$cond))
    b.star <- beta + rv %*% rnorm(ncol(rv))
    fitmis <- fit
    fitmis$fit$par[names(fitmis$fit$par)=="beta"] <- b.star
    newdatamis = data.frame(X = x[wy, , drop = FALSE])
    colnames(newdatamis) = nam
    p <- predict(fitmis, newdata = newdatamis, type = "response",
                 na.action = na.pass)
    imputed.values = rpois(length(p), lambda = p)
    if (EV){
    outliers <- extremevalues::getOutliers(imputed.values, rho = c(0.3, 
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
    }
    }
    return(imputed.values)
  }

#' @describeIn mice.impute.2l.poisson Bayesian NB regression variant; fixed intercept
#' @export
mice.impute.2l.poisson.noint<-
  function (y, ry, x, type, intercept = FALSE, wy = NULL, EV=TRUE)
  {
    if (is.null(wy)) 
      wy <- !ry
    Y = y[ry]
    X = x[ry, ]
    if (sum(type == -2) > 1) {
      stop("only one class allowed!")
    }
    grp = which(type == -2)
    nam = colnames(X)
    fixedeff = paste(paste(nam[-grp], collapse = "+"), sep = "")
    if (any(type == 2)) {
      ran = which(type == 2)
      randeff = paste("+", paste(nam[ran], collapse = "+"),
                      sep = "")
      if (!intercept) {
        randeff = paste("(0", randeff, "|", paste(nam[grp],")"),
                        sep = "")
      }
      else {
        randeff = paste("(1", randeff, "|", paste(nam[grp], ")"),
                        sep = "")
      }
    }
    else {
      if (!intercept) {
        randeff <- paste("(0", "|", paste(nam[grp]),")", sep = "")
      }
      else {
        randeff <- paste("(1", "|", paste(nam[grp]),")", sep = "")
      }
    }
    
    form = as.formula(paste("Y", "~", fixedeff, "+", randeff, sep = ""))
    dat = data.frame(Y, X)
    fit <- glmmTMB::glmmTMB(formula = form,  data = dat,
                            family=list(family="poisson",link="log"))
    
    fit.sum <- summary(fit)
    beta <- fixef(fit)$cond
    rv <- t(chol(vcov(fit)$cond))
    b.star <- beta + rv %*% rnorm(ncol(rv))
    fitmis <- fit
    fitmis$fit$par[names(fitmis$fit$par)=="beta"] <- b.star
    newdatamis = data.frame(X = x[wy, , drop = FALSE])
    colnames(newdatamis) = nam
    p <- predict(fitmis, newdata = newdatamis, type = "response",
                 na.action = na.pass)
    im = rpois(length(p), lambda = p)
    imputed.values <- im
    
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

#' @describeIn mice.impute.2l.poisson Bootstrap NB regression variant; random intercept
#' @export
mice.impute.2l.poisson.boot<-
  function (y, ry, x, type, intercept = TRUE, wy = NULL, EV=TRUE)
  {
    if (is.null(wy)) 
      wy <- !ry
    
    if (sum(type == -2) > 1) {
      stop("only one class allowed!")
    }
    grp = which(type == -2)
    
    Y = y[ry]
    X = x[ry, ]
    X = data.frame(X)
    nam=colnames(X)
    
    dat=data.frame(Y,X)
    GRPnam=nam[grp]
    colnames(dat)[colnames(dat)%in%nam[grp]]<-"GRP"
    cls <- sample(unique(dat$GRP), replace=TRUE)
    cls.col <- data.frame(GRP=cls)
    dat <- merge(cls.col, dat, by="GRP")
    colnames(dat)[colnames(dat)=="GRP"] <- GRPnam
    
    
    fixedeff = paste(paste(nam[-grp], collapse = "+"), sep = "")
    if (any(type == 2)) {
      ran = which(type == 2)
      randeff = paste("+", paste(nam[ran], collapse = "+"),
                      sep = "")
      if (!intercept) {
        randeff = paste("(0", randeff, "|", paste(nam[grp],")"),
                        sep = "")
      }
      else {
        randeff = paste("(1", randeff, "|", paste(nam[grp], ")"),
                        sep = "")
      }
    }
    else {
      if (!intercept) {
        randeff <- paste("(0", "|", paste(nam[grp]),")", sep = "")
      }
      else {
        randeff <- paste("(1", "|", paste(nam[grp]),")", sep = "")
      }
    }
    
    form = as.formula(paste("Y", "~", fixedeff, "+", randeff, sep = ""))
    fit <- glmmTMB::glmmTMB(formula = form,  data = dat,
                            family=list(family="poisson",link="log"))
    fit.sum <- summary(fit)
    newdatamis = data.frame(X = x[wy, , drop = FALSE])
    colnames(newdatamis) = nam
    p <- predict(fit, newdata = newdatamis, type = "response",
                 na.action = na.pass, allow.new.levels = TRUE)
    im = rpois(length(p), lambda = p)
    imputed.values <- im
    
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

#' @describeIn mice.impute.2l.poisson Bootstrap NB regression variant; fixed intercept
#' @export
mice.impute.2l.poisson.noint.boot<-
  function (y, ry, x, type, intercept = FALSE, wy = NULL, EV=TRUE)
  {
    if (is.null(wy)) 
      wy <- !ry
    
    if (sum(type == -2) > 1) {
      stop("only one class allowed!")
    }
    grp = which(type == -2)
    
    Y = y[ry]
    X = x[ry, ]
    X = data.frame(X)
    nam=colnames(X)
    
    dat=data.frame(Y,X)
    GRPnam=nam[grp]
    colnames(dat)[colnames(dat)%in%nam[grp]]<-"GRP"
    cls <- sample(unique(dat$GRP), replace=TRUE)
    cls.col <- data.frame(GRP=cls)
    dat <- merge(cls.col, dat, by="GRP")
    colnames(dat)[colnames(dat)=="GRP"] <- GRPnam
    
    
    fixedeff = paste(paste(nam[-grp], collapse = "+"), sep = "")
    if (any(type == 2)) {
      ran = which(type == 2)
      randeff = paste("+", paste(nam[ran], collapse = "+"),
                      sep = "")
      if (!intercept) {
        randeff = paste("(0", randeff, "|", paste(nam[grp],")"),
                        sep = "")
      }
      else {
        randeff = paste("(1", randeff, "|", paste(nam[grp], ")"),
                        sep = "")
      }
    }
    else {
      if (!intercept) {
        randeff <- paste("(0", "|", paste(nam[grp]),")", sep = "")
      }
      else {
        randeff <- paste("(1", "|", paste(nam[grp]),")", sep = "")
      }
    }
    
    form = as.formula(paste("Y", "~", fixedeff, "+", randeff, sep = ""))
    fit <- glmmTMB::glmmTMB(formula = form,  data = dat,
                            family=list(family="poisson",link="log"))
    fit.sum <- summary(fit)
    newdatamis = data.frame(X = x[wy, , drop = FALSE])
    colnames(newdatamis) = nam
    p <- predict(fit, newdata = newdatamis, type = "response",
                 na.action = na.pass, allow.new.levels = TRUE)
    im = rpois(length(p), lambda = p)
    imputed.values <- im
    
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
#### two-level NB imputation


#' @describeIn mice.impute.2l.poisson Bayesian NB regression variant; random intercept
#' @export
mice.impute.2l.nb<-
  function (y, ry, x, type, intercept = TRUE, wy = NULL, EV=TRUE)
  {
    if (is.null(wy)) 
      wy <- !ry
    Y = y[ry]
    X = x[ry, ]
    if (sum(type == -2) > 1) {
      stop("only one class allowed!")
    }
    grp = which(type == -2)
    nam = colnames(X)
    fixedeff = paste(paste(nam[-grp], collapse = "+"), sep = "")
    if (any(type == 2)) {
      ran = which(type == 2)
      randeff = paste("+", paste(nam[ran], collapse = "+"),
                      sep = "")
      if (!intercept) {
        randeff = paste("(0", randeff, "|", paste(nam[grp],")"),
                        sep = "")
      }
      else {
        randeff = paste("(1", randeff, "|", paste(nam[grp], ")"),
                        sep = "")
      }
    }
    else {
      if (!intercept) {
        randeff <- paste("(0", "|", paste(nam[grp]),")", sep = "")
      }
      else {
        randeff <- paste("(1", "|", paste(nam[grp]),")", sep = "")
      }
    }
    
    form = as.formula(paste("Y", "~", fixedeff, "+", randeff, sep = ""))
    dat = data.frame(Y, X)
    fit <- glmmTMB(formula = form,  data = dat,
                   family=list(family="nbinom2",link="log"))
    
    fit.sum <- summary(fit)
    beta <- fixef(fit)$cond
    rv <- t(chol(vcov(fit)$cond))
    b.star <- beta + rv %*% rnorm(ncol(rv))
    fitmis <- fit
    fitmis$fit$par[names(fitmis$fit$par)=="beta"] <- b.star
    newdatamis = data.frame(X = x[wy, , drop = FALSE])
    colnames(newdatamis) = nam
    p <- predict(fitmis, newdata = newdatamis, type = "response",
                 na.action = na.pass)
    imputed.values = rnegbin(length(p), theta = fit.sum$sigma, mu = p)
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
      
      }
    }
    
    return(imputed.values)
  }

#' @describeIn mice.impute.2l.poisson Bayesian NB regression variant; fixed intercept
#' @export
mice.impute.2l.nb.noint<-
  function (y, ry, x, type, intercept = FALSE, wy = NULL, EV=TRUE)
  {
    if (is.null(wy)) 
      wy <- !ry
    Y = y[ry]
    X = x[ry, ]
    if (sum(type == -2) > 1) {
      stop("only one class allowed!")
    }
    grp = which(type == -2)
    nam = colnames(X)
    fixedeff = paste(paste(nam[-grp], collapse = "+"), sep = "")
    if (any(type == 2)) {
      ran = which(type == 2)
      randeff = paste("+", paste(nam[ran], collapse = "+"),
                      sep = "")
      if (!intercept) {
        randeff = paste("(0", randeff, "|", paste(nam[grp],")"),
                        sep = "")
      }
      else {
        randeff = paste("(1", randeff, "|", paste(nam[grp], ")"),
                        sep = "")
      }
    }
    else {
      if (!intercept) {
        randeff <- paste("(0", "|", paste(nam[grp]),")", sep = "")
      }
      else {
        randeff <- paste("(1", "|", paste(nam[grp]),")", sep = "")
      }
    }
    
    form = as.formula(paste("Y", "~", fixedeff, "+", randeff, sep = ""))
    dat = data.frame(Y, X)
    fit <- glmmTMB::glmmTMB(formula = form,  data = dat,
                   family=list(family="nbinom2",link="log"))
    
    fit.sum <- summary(fit)
    beta <- fixef(fit)$cond
    rv <- t(chol(vcov(fit)$cond))
    b.star <- beta + rv %*% rnorm(ncol(rv))
    fitmis <- fit
    fitmis$fit$par[names(fitmis$fit$par)=="beta"] <- b.star
    newdatamis = data.frame(X = x[wy, , drop = FALSE])
    colnames(newdatamis) = nam
    p <- predict(fitmis, newdata = newdatamis, type = "response",
                 na.action = na.pass)
    im = rnegbin(length(p), theta = fit.sum$sigma, mu = p)
    imputed.values <- im
    
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

#' @describeIn mice.impute.2l.poisson Bootstrap NB regression variant; random intercept
#' @export
mice.impute.2l.nb.boot<-
  function (y, ry, x, type, intercept = TRUE, wy = NULL, EV=TRUE)
  {
    if (is.null(wy)) 
      wy <- !ry
    
    if (sum(type == -2) > 1) {
      stop("only one class allowed!")
    }
    grp = which(type == -2)
    
    Y = y[ry]
    X = x[ry, ]
    X = data.frame(X)
    nam=colnames(X)
    
    dat=data.frame(Y,X)
    GRPnam=nam[grp]
    colnames(dat)[colnames(dat)%in%nam[grp]]<-"GRP"
    cls <- sample(unique(dat$GRP), replace=TRUE)
    cls.col <- data.frame(GRP=cls)
    dat <- merge(cls.col, dat, by="GRP")
    colnames(dat)[colnames(dat)=="GRP"] <- GRPnam
    
    
    fixedeff = paste(paste(nam[-grp], collapse = "+"), sep = "")
    if (any(type == 2)) {
      ran = which(type == 2)
      randeff = paste("+", paste(nam[ran], collapse = "+"),
                      sep = "")
      if (!intercept) {
        randeff = paste("(0", randeff, "|", paste(nam[grp],")"),
                        sep = "")
      }
      else {
        randeff = paste("(1", randeff, "|", paste(nam[grp], ")"),
                        sep = "")
      }
    }
    else {
      if (!intercept) {
        randeff <- paste("(0", "|", paste(nam[grp]),")", sep = "")
      }
      else {
        randeff <- paste("(1", "|", paste(nam[grp]),")", sep = "")
      }
    }
    
    form = as.formula(paste("Y", "~", fixedeff, "+", randeff, sep = ""))
    fit <- glmmTMB::glmmTMB(formula = form,  data = dat,
                   family=list(family="nbinom2",link="log"))
    fit.sum <- summary(fit)
    newdatamis = data.frame(X = x[wy, , drop = FALSE])
    colnames(newdatamis) = nam
    p <- predict(fit, newdata = newdatamis, type = "response",
                 na.action = na.pass, allow.new.levels = TRUE)
    im = rnegbin(length(p), theta = fit.sum$sigma, mu = p)
    imputed.values <- im
    
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

#' @describeIn mice.impute.2l.poisson Bootstrap NB regression variant; fixed intercept
#' @export
mice.impute.2l.nb.noint.boot<-
  function (y, ry, x, type, intercept = FALSE, wy = NULL, EV=TRUE)
  {
    if (is.null(wy)) 
      wy <- !ry
    
    if (sum(type == -2) > 1) {
      stop("only one class allowed!")
    }
    grp = which(type == -2)
    
    Y = y[ry]
    X = x[ry, ]
    X = data.frame(X)
    nam=colnames(X)
    
    dat=data.frame(Y,X)
    GRPnam=nam[grp]
    colnames(dat)[colnames(dat)%in%nam[grp]]<-"GRP"
    cls <- sample(unique(dat$GRP), replace=TRUE)
    cls.col <- data.frame(GRP=cls)
    dat <- merge(cls.col, dat, by="GRP")
    colnames(dat)[colnames(dat)=="GRP"] <- GRPnam
    
    
    fixedeff = paste(paste(nam[-grp], collapse = "+"), sep = "")
    if (any(type == 2)) {
      ran = which(type == 2)
      randeff = paste("+", paste(nam[ran], collapse = "+"),
                      sep = "")
      if (!intercept) {
        randeff = paste("(0", randeff, "|", paste(nam[grp],")"),
                        sep = "")
      }
      else {
        randeff = paste("(1", randeff, "|", paste(nam[grp], ")"),
                        sep = "")
      }
    }
    else {
      if (!intercept) {
        randeff <- paste("(0", "|", paste(nam[grp]),")", sep = "")
      }
      else {
        randeff <- paste("(1", "|", paste(nam[grp]),")", sep = "")
      }
    }
    
    form = as.formula(paste("Y", "~", fixedeff, "+", randeff, sep = ""))
    fit <- glmmTMB::glmmTMB(formula = form,  data = dat,
                   family=list(family="nbinom2",link="log"))
    fit.sum <- summary(fit)
    newdatamis = data.frame(X = x[wy, , drop = FALSE])
    colnames(newdatamis) = nam
    p <- predict(fit, newdata = newdatamis, type = "response",
                 na.action = na.pass, allow.new.levels = TRUE)
    im = rnegbin(length(p), theta = fit.sum$sigma, mu = p)
    imputed.values <- im
    
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
#### backward compatibility

#' @describeIn mice.impute.2l.poisson identical to \code{mice.impute.2l.nb}; kept for backward compatibility
#' @export
mice.impute.2l.nb2<-
  function (y, ry, x, type, intercept = TRUE, wy = NULL, EV=TRUE)
  {
    if (is.null(wy)) 
      wy <- !ry
    
    Y = y[ry]
    X = x[ry, ]
    if (sum(type == -2) > 1) {
      stop("only one class allowed!")
    }
    grp = which(type == -2)
    nam = colnames(X)
    fixedeff = paste(paste(nam[-grp], collapse = "+"), sep = "")
    if (any(type == 2)) {
      ran = which(type == 2)
      randeff = paste("+", paste(nam[ran], collapse = "+"),
                      sep = "")
      if (!intercept) {
        randeff = paste("(0", randeff, "|", paste(nam[grp],")"),
                        sep = "")
      }
      else {
        randeff = paste("(1", randeff, "|", paste(nam[grp], ")"),
                        sep = "")
      }
    }
    else {
      if (!intercept) {
        randeff <- paste("(0", "|", paste(nam[grp]),")", sep = "")
      }
      else {
        randeff <- paste("(1", "|", paste(nam[grp]),")", sep = "")
      }
    }
    
    form = as.formula(paste("Y", "~", fixedeff, "+", randeff, sep = ""))
    dat = data.frame(Y, X)
    fit <- glmmTMB::glmmTMB(formula = form,  data = dat,
                   family=list(family="nbinom2",link="log"))
    
    fit.sum <- summary(fit)
    beta <- fixef(fit)$cond
    rv <- t(chol(vcov(fit)$cond))
    b.star <- beta + rv %*% rnorm(ncol(rv))
    fitmis <- fit
    fitmis$fit$par[names(fitmis$fit$par)=="beta"] <- b.star
    newdatamis = data.frame(X = x[wy, , drop = FALSE])
    colnames(newdatamis) = nam
    p <- predict(fitmis, newdata = newdatamis, type = "response",
                 na.action = na.pass)
    im = rnegbin(length(p), theta = fit.sum$sigma, mu = p)
    imputed.values <- im
    
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

#' @describeIn mice.impute.2l.poisson identical to \code{mice.impute.2l.nb.noint}; kept for backward compatibility
#' @export
mice.impute.2l.nb2.noint<-
  function (y, ry, x, type, intercept = FALSE, wy = NULL, EV=TRUE)
  {
    if (is.null(wy)) 
      wy <- !ry
    
    Y = y[ry]
    X = x[ry, ]
    if (sum(type == -2) > 1) {
      stop("only one class allowed!")
    }
    grp = which(type == -2)
    nam = colnames(X)
    fixedeff = paste(paste(nam[-grp], collapse = "+"), sep = "")
    if (any(type == 2)) {
      ran = which(type == 2)
      randeff = paste("+", paste(nam[ran], collapse = "+"),
                      sep = "")
      if (!intercept) {
        randeff = paste("(0", randeff, "|", paste(nam[grp],")"),
                        sep = "")
      }
      else {
        randeff = paste("(1", randeff, "|", paste(nam[grp], ")"),
                        sep = "")
      }
    }
    else {
      if (!intercept) {
        randeff <- paste("(0", "|", paste(nam[grp]),")", sep = "")
      }
      else {
        randeff <- paste("(1", "|", paste(nam[grp]),")", sep = "")
      }
    }
    
    form = as.formula(paste("Y", "~", fixedeff, "+", randeff, sep = ""))
    dat = data.frame(Y, X)
    fit <- glmmTMB::glmmTMB(formula = form,  data = dat,
                   family=list(family="nbinom2",link="log"))
    
    fit.sum <- summary(fit)
    beta <- fixef(fit)$cond
    rv <- t(chol(vcov(fit)$cond))
    b.star <- beta + rv %*% rnorm(ncol(rv))
    fitmis <- fit
    fitmis$fit$par[names(fitmis$fit$par)=="beta"] <- b.star
    newdatamis = data.frame(X = x[wy, , drop = FALSE])
    colnames(newdatamis) = nam
    p <- predict(fitmis, newdata = newdatamis, type = "response",
                 na.action = na.pass)
    im = rnegbin(length(p), theta = fit.sum$sigma, mu = p)
    imputed.values <- im
    
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

#' @describeIn mice.impute.2l.poisson identical to \code{mice.impute.2l.nb.boot}; kept for backward compatibility
#' @export
mice.impute.2l.nb2.boot<-
  function (y, ry, x, type, intercept = TRUE, wy = NULL, EV=TRUE)
  {
    if (is.null(wy)) 
      wy <- !ry
    
    if (sum(type == -2) > 1) {
      stop("only one class allowed!")
    }
    grp = which(type == -2)
    
    Y = y[ry]
    X = x[ry, ]
    X = data.frame(X)
    nam=colnames(X)
    
    dat=data.frame(Y,X)
    GRPnam=nam[grp]
    colnames(dat)[colnames(dat)%in%nam[grp]]<-"GRP"
    cls <- sample(unique(dat$GRP), replace=TRUE)
    cls.col <- data.frame(GRP=cls)
    dat <- merge(cls.col, dat, by="GRP")
    colnames(dat)[colnames(dat)=="GRP"] <- GRPnam
    
    
    fixedeff = paste(paste(nam[-grp], collapse = "+"), sep = "")
    if (any(type == 2)) {
      ran = which(type == 2)
      randeff = paste("+", paste(nam[ran], collapse = "+"),
                      sep = "")
      if (!intercept) {
        randeff = paste("(0", randeff, "|", paste(nam[grp],")"),
                        sep = "")
      }
      else {
        randeff = paste("(1", randeff, "|", paste(nam[grp], ")"),
                        sep = "")
      }
    }
    else {
      if (!intercept) {
        randeff <- paste("(0", "|", paste(nam[grp]),")", sep = "")
      }
      else {
        randeff <- paste("(1", "|", paste(nam[grp]),")", sep = "")
      }
    }
    
    form = as.formula(paste("Y", "~", fixedeff, "+", randeff, sep = ""))
    fit <- glmmTMB::glmmTMB(formula = form,  data = dat,
                   family=list(family="nbinom2",link="log"))
    fit.sum <- summary(fit)
    newdatamis = data.frame(X = x[wy, , drop = FALSE])
    colnames(newdatamis) = nam
    p <- predict(fit, newdata = newdatamis, type = "response",
                 na.action = na.pass, allow.new.levels = TRUE)
    im = rnegbin(length(p), theta = fit.sum$sigma, mu = p)
    imputed.values <- im
    
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

#' @describeIn mice.impute.2l.poisson identical to \code{mice.impute.2l.nb.noint.boot}; kept for backward compatibility
#' @export
mice.impute.2l.nb2.noint.boot<-
  function (y, ry, x, type, intercept = FALSE, wy = NULL, EV=TRUE)
  {
    if (is.null(wy)) 
      wy <- !ry
    
    if (sum(type == -2) > 1) {
      stop("only one class allowed!")
    }
    grp = which(type == -2)
    
    Y = y[ry]
    X = x[ry, ]
    X = data.frame(X)
    nam=colnames(X)
    
    dat=data.frame(Y,X)
    GRPnam=nam[grp]
    colnames(dat)[colnames(dat)%in%nam[grp]]<-"GRP"
    cls <- sample(unique(dat$GRP), replace=TRUE)
    cls.col <- data.frame(GRP=cls)
    dat <- merge(cls.col, dat, by="GRP")
    colnames(dat)[colnames(dat)=="GRP"] <- GRPnam
    
    
    fixedeff = paste(paste(nam[-grp], collapse = "+"), sep = "")
    if (any(type == 2)) {
      ran = which(type == 2)
      randeff = paste("+", paste(nam[ran], collapse = "+"),
                      sep = "")
      if (!intercept) {
        randeff = paste("(0", randeff, "|", paste(nam[grp],")"),
                        sep = "")
      }
      else {
        randeff = paste("(1", randeff, "|", paste(nam[grp], ")"),
                        sep = "")
      }
    }
    else {
      if (!intercept) {
        randeff <- paste("(0", "|", paste(nam[grp]),")", sep = "")
      }
      else {
        randeff <- paste("(1", "|", paste(nam[grp]),")", sep = "")
      }
    }
    
    form = as.formula(paste("Y", "~", fixedeff, "+", randeff, sep = ""))
    fit <- glmmTMB::glmmTMB(formula = form,  data = dat,
                   family=list(family="nbinom2",link="log"))
    fit.sum <- summary(fit)
    newdatamis = data.frame(X = x[wy, , drop = FALSE])
    colnames(newdatamis) = nam
    p <- predict(fit, newdata = newdatamis, type = "response",
                 na.action = na.pass, allow.new.levels = TRUE)
    im = rnegbin(length(p), theta = fit.sum$sigma, mu = p)
    imputed.values <- im
    
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