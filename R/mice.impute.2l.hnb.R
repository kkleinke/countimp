#' Multiple Imputation of Zero-Inflated Two-Level Count Data
#'
#' The functions impute zero-inflated multilevel count data based on a two-level Poisson or negative binomial hurdle model, either using a Bayesian regression or a bootstrap regression approach (appendix: ``\code{.boot}''). The \code{.noint} variants treat the intercept only as a fixed, but \emph{not} as a random effect. It may be specified, if the intercept is excluded from the random part of the zero model (``\code{.noint.zero}''), the count model (``\code{.noint.count}''), or both models (``\code{.noint.both}''). Hurdle models are mixture models and consist of two model components: the zero model (a binomial generalized linear mixed effects model), determining, if the observational unit has a zero or non-zero value, and the count model (a zero-truncated two-level Poisson or NB model), determining, what non-zero value the observational unit has.
#'
#' Model specification details:
#'    \itemize{
#'      \item -2 = class variable (only one class variable is allowed!)
#'      \item 0 = variable not included in imputation model
#'      \item 1 = variable will be included as a fixed effect (zero \emph{and} count model)
#'      \item 2 = variable will be included as a fixed \emph{and} random effect (zero \emph{and} count model)
#'      \item 3 = variable will be included as a fixed effect (count model only)
#'      \item 4 = variable will be included as a fixed \emph{and} random effect (count model only)
#'      \item 5 = variable will be included as a fixed effect (zero model only)
#'      \item 6 = variable will be included as a fixed \emph{and} random effect (zero model only)
#'    }
#' The Bayesian regression variants (see Rubin 1987, p. 169-170) consist of the following steps:
#'  \enumerate{
#'  \item Fit the zero model (a two-level binomial generalized linear mixed effects model), using the \code{glmmTMB} function from package \pkg{glmmTMB}; find bhat, the posterior mean, and V(bhat), the posterior variance of model parameters b
#'  \item Draw b* from N(bhat,V(bhat))
#'  \item Compute predicted probabilities for having a zero vs. non-zero count
#'  \item Draw imputations (zeros and ones) from a binomial distribution with the respective individual probabilities obtained from step 3.
#'  \item Fit the count model (a zero-truncated two-level Poisson or NB model) using the \code{glmmTMB} function from package \pkg{glmmTMB}; find bhat, the posterior mean, and V(bhat), the posterior variance of model parameters b.
#'  \item  Draw b* from N(bhat,V(bhat))
#'  \item Compute predicted values using parameters b* and replace non-zero imputations (from step 4) by a draw from a zero-truncated NB distribution with mean parameter mu being the count predicted for the respective incomplete case.
#'  }
#' The bootstrap functions draw a bootstrap sample from \code{y[ry]} and \code{x[ry,]} (Note: we resample clusters rather than individual cases) and consist of the following steps:
#' \enumerate{
#' \item Fit the zero model to the bootstrap sample
#' \item Compute predicted probabilities for having a zero vs. non-zero count
#' \item Draw imputations from a binomial distribution.
#' \item Fit the count model to the boostrap sample
#' \item Compute predicted counts and draw non-zero imputations (from step 3) from a zero-truncated Poisson or NB distribution.
#' }
#' @param y Numeric vector with incomplete data in long format (i.e. the groups are stacked upon each other)
#' @param ry Response pattern of \code{y} (\code{TRUE}=observed, \code{FALSE}=missing)
#' @param x matrix with \code{length(y)} rows containing complete covariates; also in long format
#' @param type vector of length \code{ncol(x)} identifying fixed, random, and class variables; \code{type} is automatically extracted from the \code{predictorMatrix}; see \pkg{mice}'s user's manual for details about how to specify the imputation model; see also section ``details''.
#' @param intercept.c \code{TRUE}: model will include intercept as a random effect in the count model; \code{FALSE}: count model intercept will be treated as fixed.
#' @param intercept.z \code{TRUE}: model will include intercept as a random effect in the zero model; \code{FALSE}: zero model intercept will be treated as fixed.
#' @param wy Logical vector of length \code{length(y)}. A \code{TRUE} value indicates locations in \code{y} for which imputations are created. Default is \code{!ry}
#' @param EV should automatic outlier handling of imputed values be enabled?  Default is \code{TRUE}: extreme imputations will be identified. These values will be replaced by imputations obtained by predictive mean matching (function \code{mice.impute.midastouch()})
#' @return Numeric vector of length \code{sum(!ry)} with imputations
#' @export
#' @import glmmTMB aster
#' @importFrom stats as.formula na.pass predict rbinom rnorm vcov
#' @aliases mice.impute.2l.hnb mice.impute.2l.hnb.boot mice.impute.2l.hnb.noint.both mice.impute.2l.hnb.noint.both.boot mice.impute.2l.hnb.noint.zero mice.impute.2l.hnb.noint.zero.boot mice.impute.2l.hnb.noint.count mice.impute.2l.hnb.noint.count.boot
#' @aliases mice.impute.2l.hp mice.impute.2l.hp.boot mice.impute.2l.hp.noint.both mice.impute.2l.hp.noint.both.boot mice.impute.2l.hp.noint.zero mice.impute.2l.hp.noint.zero.boot mice.impute.2l.hp.noint.count mice.impute.2l.hp.noint.count.boot
#' @author Kristian Kleinke
#' @describeIn mice.impute.2l.hnb Bayesian regression variant; random intercepts
mice.impute.2l.hnb <-
  function (y, ry, x, type, intercept.c = TRUE, intercept.z = TRUE, wy = NULL, EV=TRUE)
  {
    if (is.null(wy)) 
      wy <- !ry
    if (sum(type == -2) > 1) {
      stop("only one class allowed!")
    }
    
    Y = y[ry]
    X = x[ry, ]
    X = data.frame(X)
    nam = colnames(X)
    grp <- which(type == -2)
    dat=data.frame(Y,X)
    
    fc <- which(type == 3)
    rc <- which(type == 4)
    fz <- which(type == 5)
    rz <- which(type == 6)
    f <- which(type == 1)
    r <- which(type == 2)
    
    ran <- c(r, rz)
    ran <- unique(ran)
    ran <- sort(ran)
    randeff <- paste(nam[ran], collapse = "+")
    if (any(type == 2) | any(type == 6)) {
      randeff = paste("+", paste(nam[ran], collapse = "+"),
                      sep = "")
      if (!intercept.c) {
        randeff = paste("(0", randeff, "|", paste(nam[grp],")"),
                        sep = "")
      }
      else {
        randeff = paste("(1", randeff, "|", paste(nam[grp], ")"),
                        sep = "")
      }
    }
    else {
      if (!intercept.c) {
        randeff <- paste("(0", "|", paste(nam[grp]),")", sep = "")
      }
      else {
        randeff <- paste("(1", "|", paste(nam[grp]),")", sep = "")
      }
    }
    fix <- c(f, r, fz, rz)
    fix <- unique(fix)
    fix <- sort(fix)
    fixedeff <- paste(nam[fix], collapse = "+")
    
    nz <- as.numeric(dat$Y == 0)
    dat$nz = nz
    
    form = as.formula(paste("nz", "~", fixedeff, "+", randeff, sep = ""))
    fit <- glmmTMB::glmmTMB(formula = form,  data = dat,
                   family="binomial")
   
    fit.sum <- summary(fit)
    beta <- fixef(fit)$cond
    rv <- t(chol(vcov(fit)$cond))
    b.star <- beta + rv %*% rnorm(ncol(rv))
    fitmis <- fit
    fitmis$fit$par[names(fitmis$fit$par)=="beta"] <- b.star
    
    newdatamis = data.frame(X = x[wy, , drop = FALSE])
    colnames(newdatamis) = nam
    p <- predict(fitmis, newdata = newdatamis, type = "response",
                 na.action = na.pass, allow.new.levels = TRUE)
    im <- rbinom(n = length(p), size = 1, prob = 1 - p)
    nonzero <- which(im == 1)
    
    ran <- c(r, rc)
    ran <- unique(ran)
    ran <- sort(ran)
    randeff <- paste(nam[ran], collapse = "+")
    if (any(type == 2) | any(type == 4)) {
      randeff = paste("+", paste(nam[ran], collapse = "+"),
                      sep = "")
      if (!intercept.c) {
        randeff = paste("(0", randeff, "|", paste(nam[grp],")"),
                        sep = "")
      }
      else {
        randeff = paste("(1", randeff, "|", paste(nam[grp], ")"),
                        sep = "")
      }
    }
    else {
      if (!intercept.c) {
        randeff <- paste("(0", "|", paste(nam[grp]),")", sep = "")
      }
      else {
        randeff <- paste("(1", "|", paste(nam[grp]),")", sep = "")
      }
    }
    fix <- c(f, r, fc, rc)
    fix <- unique(fix)
    fix <- sort(fix)
    fixedeff = paste(paste(nam[-grp], collapse = "+"), sep = "")
    form = as.formula(paste("Y", "~", fixedeff, "+", randeff, sep = ""))
    fit <- glmmTMB::glmmTMB(formula = form,  data = subset(dat,
                                                  dat[, "Y"] > 0),
                   family=list(family="truncated_nbinom2",link="log"))
    fit.sum=summary(fit)
    
    fit.sum <- summary(fit)
    beta <- fixef(fit)$cond
    rv <- t(chol(vcov(fit)$cond))
    b.star <- beta + rv %*% rnorm(ncol(rv))
    fitmis <- fit
    fitmis$fit$par[names(fitmis$fit$par)=="beta"] <- b.star
    
    newdatamis = data.frame(X = x[wy, , drop = FALSE])
    colnames(newdatamis) = nam
    p <- predict(fitmis, newdata = newdatamis, type = "response",
                 na.action = na.pass, allow.new.levels = TRUE)
    
    imc = aster::rktnb(length(p), size = fit.sum$sigma, k = 0, mu = p, xpred = 1)
    
    im[nonzero] <- imc[nonzero]
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

#' @export
#' @describeIn  mice.impute.2l.hnb bootstrap regression variant; random intercepts
mice.impute.2l.hnb.boot<-
  function (y, ry, x, type, intercept.c = TRUE, intercept.z = TRUE, wy = NULL, EV=TRUE)
  {
    if (is.null(wy)) 
      wy <- !ry
    
    if (sum(type == -2) > 1) {
      stop("only one class allowed!")
    }
    
    Y = y[ry]
    X = x[ry, ]
    X = data.frame(X)
    nam = colnames(X)
    grp <- which(type == -2)
    dat=data.frame(Y,X)
    GRPnam=nam[grp]
    colnames(dat)[colnames(dat)%in%nam[grp]]<-"GRP"
    cls <- sample(unique(dat$GRP),replace=TRUE)
    cls.col <- data.frame(GRP=cls)
    dat <- merge(cls.col,dat,by="GRP")
    colnames(dat)[colnames(dat)=="GRP"] <- GRPnam
    
    fc <- which(type == 3)
    rc <- which(type == 4)
    fz <- which(type == 5)
    rz <- which(type == 6)
    f <- which(type == 1)
    r <- which(type == 2)
    
    ran <- c(r, rz)
    ran <- unique(ran)
    ran <- sort(ran)
    randeff <- paste(nam[ran], collapse = "+")
    if (any(type == 2) | any(type == 6)) {
      randeff = paste("+", paste(nam[ran], collapse = "+"),
                      sep = "")
      if (!intercept.c) {
        randeff = paste("(0", randeff, "|", paste(nam[grp],")"),
                        sep = "")
      }
      else {
        randeff = paste("(1", randeff, "|", paste(nam[grp], ")"),
                        sep = "")
      }
    }
    else {
      if (!intercept.c) {
        randeff <- paste("(0", "|", paste(nam[grp]),")", sep = "")
      }
      else {
        randeff <- paste("(1", "|", paste(nam[grp]),")", sep = "")
      }
    }
    fix <- c(f, r, fz, rz)
    fix <- unique(fix)
    fix <- sort(fix)
    fixedeff <- paste(nam[fix], collapse = "+")
    
    nz <- as.numeric(dat$Y == 0)
    dat$nz = nz
    
    form = as.formula(paste("nz", "~", fixedeff, "+", randeff, sep = ""))
    fit <- glmmTMB::glmmTMB(formula = form,  data = dat,
                   family="binomial")
    
    
    newdatamis = data.frame(X = x[wy, , drop = FALSE])
    colnames(newdatamis) = nam
    p <- predict(fit, newdata = newdatamis, type = "response",
                 na.action = na.pass, allow.new.levels = TRUE)
    im <- rbinom(n = length(p), size = 1, prob = 1 - p)
    nonzero <- which(im == 1)
    
    ran <- c(r, rc)
    ran <- unique(ran)
    ran <- sort(ran)
    randeff <- paste(nam[ran], collapse = "+")
    if (any(type == 2) | any(type == 4)) {
      randeff = paste("+", paste(nam[ran], collapse = "+"),
                      sep = "")
      if (!intercept.c) {
        randeff = paste("(0", randeff, "|", paste(nam[grp],")"),
                        sep = "")
      }
      else {
        randeff = paste("(1", randeff, "|", paste(nam[grp], ")"),
                        sep = "")
      }
    }
    else {
      if (!intercept.c) {
        randeff <- paste("(0", "|", paste(nam[grp]),")", sep = "")
      }
      else {
        randeff <- paste("(1", "|", paste(nam[grp]),")", sep = "")
      }
    }
    fix <- c(f, r, fc, rc)
    fix <- unique(fix)
    fix <- sort(fix)
    fixedeff = paste(paste(nam[-grp], collapse = "+"), sep = "")
    form = as.formula(paste("Y", "~", fixedeff, "+", randeff, sep = ""))
    fit <- glmmTMB::glmmTMB(formula = form,  data = subset(dat,
                                                  dat[, "Y"] > 0),
                   family=list(family="truncated_nbinom2",link="log"))
    fit.sum=summary(fit)
    newdatamis = data.frame(X = x[wy, , drop = FALSE])
    colnames(newdatamis) = nam
    p <- predict(fit, newdata = newdatamis, type = "response",
                 na.action = na.pass, allow.new.levels = TRUE)
    
    imc = aster::rktnb(length(p), size = fit.sum$sigma, k = 0, mu = p, xpred = 1)
    
    im[nonzero] <- imc[nonzero]
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

#' @export
#' @describeIn  mice.impute.2l.hnb bootstrap regression variant; fixed intercepts
mice.impute.2l.hnb.noint.both.boot<-
  function (y, ry, x, type, intercept.c = FALSE, intercept.z = FALSE, wy = NULL, EV=TRUE)
  {
    if (is.null(wy)) 
      wy <- !ry
    
    if (sum(type == -2) > 1) {
      stop("only one class allowed!")
    }
    
    Y = y[ry]
    X = x[ry, ]
    X = data.frame(X)
    nam = colnames(X)
    grp <- which(type == -2)
    dat=data.frame(Y,X)
    GRPnam=nam[grp]
    colnames(dat)[colnames(dat)%in%nam[grp]]<-"GRP"
    cls <- sample(unique(dat$GRP),replace=TRUE)
    cls.col <- data.frame(GRP=cls)
    dat <- merge(cls.col,dat,by="GRP")
    colnames(dat)[colnames(dat)=="GRP"] <- GRPnam
    
    fc <- which(type == 3)
    rc <- which(type == 4)
    fz <- which(type == 5)
    rz <- which(type == 6)
    f <- which(type == 1)
    r <- which(type == 2)
    
    ran <- c(r, rz)
    ran <- unique(ran)
    ran <- sort(ran)
    randeff <- paste(nam[ran], collapse = "+")
    if (any(type == 2) | any(type == 6)) {
      randeff = paste("+", paste(nam[ran], collapse = "+"),
                      sep = "")
      if (!intercept.c) {
        randeff = paste("(0", randeff, "|", paste(nam[grp],")"),
                        sep = "")
      }
      else {
        randeff = paste("(1", randeff, "|", paste(nam[grp], ")"),
                        sep = "")
      }
    }
    else {
      if (!intercept.c) {
        randeff <- paste("(0", "|", paste(nam[grp]),")", sep = "")
      }
      else {
        randeff <- paste("(1", "|", paste(nam[grp]),")", sep = "")
      }
    }
    fix <- c(f, r, fz, rz)
    fix <- unique(fix)
    fix <- sort(fix)
    fixedeff <- paste(nam[fix], collapse = "+")
    
    nz <- as.numeric(dat$Y == 0)
    dat$nz = nz
    
    form = as.formula(paste("nz", "~", fixedeff, "+", randeff, sep = ""))
    fit <- glmmTMB::glmmTMB(formula = form,  data = dat,
                   family="binomial")
    
    
    newdatamis = data.frame(X = x[wy, , drop = FALSE])
    colnames(newdatamis) = nam
    p <- predict(fit, newdata = newdatamis, type = "response",
                 na.action = na.pass, allow.new.levels = TRUE)
    im <- rbinom(n = length(p), size = 1, prob = 1 - p)
    nonzero <- which(im == 1)
    
    ran <- c(r, rc)
    ran <- unique(ran)
    ran <- sort(ran)
    randeff <- paste(nam[ran], collapse = "+")
    if (any(type == 2) | any(type == 4)) {
      randeff = paste("+", paste(nam[ran], collapse = "+"),
                      sep = "")
      if (!intercept.c) {
        randeff = paste("(0", randeff, "|", paste(nam[grp],")"),
                        sep = "")
      }
      else {
        randeff = paste("(1", randeff, "|", paste(nam[grp], ")"),
                        sep = "")
      }
    }
    else {
      if (!intercept.c) {
        randeff <- paste("(0", "|", paste(nam[grp]),")", sep = "")
      }
      else {
        randeff <- paste("(1", "|", paste(nam[grp]),")", sep = "")
      }
    }
    fix <- c(f, r, fc, rc)
    fix <- unique(fix)
    fix <- sort(fix)
    fixedeff = paste(paste(nam[-grp], collapse = "+"), sep = "")
    form = as.formula(paste("Y", "~", fixedeff, "+", randeff, sep = ""))
    fit <- glmmTMB::glmmTMB(formula = form,  data = subset(dat,
                                                  dat[, "Y"] > 0),
                   family=list(family="truncated_nbinom2",link="log"))
    fit.sum=summary(fit)
    newdatamis = data.frame(X = x[wy, , drop = FALSE])
    colnames(newdatamis) = nam
    p <- predict(fit, newdata = newdatamis, type = "response",
                 na.action = na.pass, allow.new.levels = TRUE)
    
    imc = aster::rktnb(length(p), size = fit.sum$sigma, k = 0, mu = p, xpred = 1)
    
    im[nonzero] <- imc[nonzero]
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

#' @export
#' @describeIn  mice.impute.2l.hnb Bayesian regression variant; fixed intercepts
mice.impute.2l.hnb.noint.both <-
  function (y, ry, x, type, intercept.c = FALSE, intercept.z = FALSE, wy = NULL, EV=TRUE)
  {
    if (is.null(wy)) 
      wy <- !ry
    if (sum(type == -2) > 1) {
      stop("only one class allowed!")
    }
    
    Y = y[ry]
    X = x[ry, ]
    X = data.frame(X)
    nam = colnames(X)
    grp <- which(type == -2)
    dat=data.frame(Y,X)
    
    fc <- which(type == 3)
    rc <- which(type == 4)
    fz <- which(type == 5)
    rz <- which(type == 6)
    f <- which(type == 1)
    r <- which(type == 2)
    
    ran <- c(r, rz)
    ran <- unique(ran)
    ran <- sort(ran)
    randeff <- paste(nam[ran], collapse = "+")
    if (any(type == 2) | any(type == 6)) {
      randeff = paste("+", paste(nam[ran], collapse = "+"),
                      sep = "")
      if (!intercept.c) {
        randeff = paste("(0", randeff, "|", paste(nam[grp],")"),
                        sep = "")
      }
      else {
        randeff = paste("(1", randeff, "|", paste(nam[grp], ")"),
                        sep = "")
      }
    }
    else {
      if (!intercept.c) {
        randeff <- paste("(0", "|", paste(nam[grp]),")", sep = "")
      }
      else {
        randeff <- paste("(1", "|", paste(nam[grp]),")", sep = "")
      }
    }
    fix <- c(f, r, fz, rz)
    fix <- unique(fix)
    fix <- sort(fix)
    fixedeff <- paste(nam[fix], collapse = "+")
    
    nz <- as.numeric(dat$Y == 0)
    dat$nz = nz
    
    form = as.formula(paste("nz", "~", fixedeff, "+", randeff, sep = ""))
    fit <- glmmTMB::glmmTMB(formula = form,  data = dat,
                   family="binomial")
    
    fit.sum <- summary(fit)
    beta <- fixef(fit)$cond
    rv <- t(chol(vcov(fit)$cond))
    b.star <- beta + rv %*% rnorm(ncol(rv))
    fitmis <- fit
    fitmis$fit$par[names(fitmis$fit$par)=="beta"] <- b.star
    
    newdatamis = data.frame(X = x[wy, , drop = FALSE])
    colnames(newdatamis) = nam
    p <- predict(fitmis, newdata = newdatamis, type = "response",
                 na.action = na.pass, allow.new.levels = TRUE)
    im <- rbinom(n = length(p), size = 1, prob = 1 - p)
    nonzero <- which(im == 1)
    
    ran <- c(r, rc)
    ran <- unique(ran)
    ran <- sort(ran)
    randeff <- paste(nam[ran], collapse = "+")
    if (any(type == 2) | any(type == 4)) {
      randeff = paste("+", paste(nam[ran], collapse = "+"),
                      sep = "")
      if (!intercept.c) {
        randeff = paste("(0", randeff, "|", paste(nam[grp],")"),
                        sep = "")
      }
      else {
        randeff = paste("(1", randeff, "|", paste(nam[grp], ")"),
                        sep = "")
      }
    }
    else {
      if (!intercept.c) {
        randeff <- paste("(0", "|", paste(nam[grp]),")", sep = "")
      }
      else {
        randeff <- paste("(1", "|", paste(nam[grp]),")", sep = "")
      }
    }
    fix <- c(f, r, fc, rc)
    fix <- unique(fix)
    fix <- sort(fix)
    fixedeff = paste(paste(nam[-grp], collapse = "+"), sep = "")
    form = as.formula(paste("Y", "~", fixedeff, "+", randeff, sep = ""))
    fit <- glmmTMB::glmmTMB(formula = form,  data = subset(dat,
                                                  dat[, "Y"] > 0),
                   family=list(family="truncated_nbinom2",link="log"))
    fit.sum=summary(fit)
    
    fit.sum <- summary(fit)
    beta <- fixef(fit)$cond
    rv <- t(chol(vcov(fit)$cond))
    b.star <- beta + rv %*% rnorm(ncol(rv))
    fitmis <- fit
    fitmis$fit$par[names(fitmis$fit$par)=="beta"] <- b.star
    
    newdatamis = data.frame(X = x[wy, , drop = FALSE])
    colnames(newdatamis) = nam
    p <- predict(fitmis, newdata = newdatamis, type = "response",
                 na.action = na.pass, allow.new.levels = TRUE)
    
    imc = aster::rktnb(length(p), size = fit.sum$sigma, k = 0, mu = p, xpred = 1)
    
    im[nonzero] <- imc[nonzero]
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

#' @export
#' @describeIn  mice.impute.2l.hnb bootstrap regression variant; fixed intercept in count model
mice.impute.2l.hnb.noint.count.boot<-
  function (y, ry, x, type, intercept.c = FALSE, intercept.z = TRUE, wy = NULL, EV=TRUE)
  {
    if (is.null(wy)) 
      wy <- !ry
    
    if (sum(type == -2) > 1) {
      stop("only one class allowed!")
    }
    
    Y = y[ry]
    X = x[ry, ]
    X = data.frame(X)
    nam = colnames(X)
    grp <- which(type == -2)
    dat=data.frame(Y,X)
    GRPnam=nam[grp]
    colnames(dat)[colnames(dat)%in%nam[grp]]<-"GRP"
    cls <- sample(unique(dat$GRP),replace=TRUE)
    cls.col <- data.frame(GRP=cls)
    dat <- merge(cls.col,dat,by="GRP")
    colnames(dat)[colnames(dat)=="GRP"] <- GRPnam
    
    fc <- which(type == 3)
    rc <- which(type == 4)
    fz <- which(type == 5)
    rz <- which(type == 6)
    f <- which(type == 1)
    r <- which(type == 2)
    
    ran <- c(r, rz)
    ran <- unique(ran)
    ran <- sort(ran)
    randeff <- paste(nam[ran], collapse = "+")
    if (any(type == 2) | any(type == 6)) {
      randeff = paste("+", paste(nam[ran], collapse = "+"),
                      sep = "")
      if (!intercept.c) {
        randeff = paste("(0", randeff, "|", paste(nam[grp],")"),
                        sep = "")
      }
      else {
        randeff = paste("(1", randeff, "|", paste(nam[grp], ")"),
                        sep = "")
      }
    }
    else {
      if (!intercept.c) {
        randeff <- paste("(0", "|", paste(nam[grp]),")", sep = "")
      }
      else {
        randeff <- paste("(1", "|", paste(nam[grp]),")", sep = "")
      }
    }
    fix <- c(f, r, fz, rz)
    fix <- unique(fix)
    fix <- sort(fix)
    fixedeff <- paste(nam[fix], collapse = "+")
    
    nz <- as.numeric(dat$Y == 0)
    dat$nz = nz
    
    form = as.formula(paste("nz", "~", fixedeff, "+", randeff, sep = ""))
    fit <- glmmTMB::glmmTMB(formula = form,  data = dat,
                   family="binomial")
    
    
    newdatamis = data.frame(X = x[wy, , drop = FALSE])
    colnames(newdatamis) = nam
    p <- predict(fit, newdata = newdatamis, type = "response",
                 na.action = na.pass, allow.new.levels = TRUE)
    im <- rbinom(n = length(p), size = 1, prob = 1 - p)
    nonzero <- which(im == 1)
    
    ran <- c(r, rc)
    ran <- unique(ran)
    ran <- sort(ran)
    randeff <- paste(nam[ran], collapse = "+")
    if (any(type == 2) | any(type == 4)) {
      randeff = paste("+", paste(nam[ran], collapse = "+"),
                      sep = "")
      if (!intercept.c) {
        randeff = paste("(0", randeff, "|", paste(nam[grp],")"),
                        sep = "")
      }
      else {
        randeff = paste("(1", randeff, "|", paste(nam[grp], ")"),
                        sep = "")
      }
    }
    else {
      if (!intercept.c) {
        randeff <- paste("(0", "|", paste(nam[grp]),")", sep = "")
      }
      else {
        randeff <- paste("(1", "|", paste(nam[grp]),")", sep = "")
      }
    }
    fix <- c(f, r, fc, rc)
    fix <- unique(fix)
    fix <- sort(fix)
    fixedeff = paste(paste(nam[-grp], collapse = "+"), sep = "")
    form = as.formula(paste("Y", "~", fixedeff, "+", randeff, sep = ""))
    fit <- glmmTMB::glmmTMB(formula = form,  data = subset(dat,
                                                  dat[, "Y"] > 0),
                   family=list(family="truncated_nbinom2",link="log"))
    fit.sum=summary(fit)
    newdatamis = data.frame(X = x[wy, , drop = FALSE])
    colnames(newdatamis) = nam
    p <- predict(fit, newdata = newdatamis, type = "response",
                 na.action = na.pass, allow.new.levels = TRUE)
    
    imc = aster::rktnb(length(p), size = fit.sum$sigma, k = 0, mu = p, xpred = 1)
    
    im[nonzero] <- imc[nonzero]
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

#' @export
#' @describeIn  mice.impute.2l.hnb Bayesian regression variant; fixed intercept in count model
mice.impute.2l.hnb.noint.count <-
  function (y, ry, x, type, intercept.c = FALSE, intercept.z = TRUE, wy = NULL, EV=TRUE)
  {
    if (is.null(wy)) 
      wy <- !ry
    if (sum(type == -2) > 1) {
      stop("only one class allowed!")
    }
    
    Y = y[ry]
    X = x[ry, ]
    X = data.frame(X)
    nam = colnames(X)
    grp <- which(type == -2)
    dat=data.frame(Y,X)
    
    fc <- which(type == 3)
    rc <- which(type == 4)
    fz <- which(type == 5)
    rz <- which(type == 6)
    f <- which(type == 1)
    r <- which(type == 2)
    
    ran <- c(r, rz)
    ran <- unique(ran)
    ran <- sort(ran)
    randeff <- paste(nam[ran], collapse = "+")
    if (any(type == 2) | any(type == 6)) {
      randeff = paste("+", paste(nam[ran], collapse = "+"),
                      sep = "")
      if (!intercept.c) {
        randeff = paste("(0", randeff, "|", paste(nam[grp],")"),
                        sep = "")
      }
      else {
        randeff = paste("(1", randeff, "|", paste(nam[grp], ")"),
                        sep = "")
      }
    }
    else {
      if (!intercept.c) {
        randeff <- paste("(0", "|", paste(nam[grp]),")", sep = "")
      }
      else {
        randeff <- paste("(1", "|", paste(nam[grp]),")", sep = "")
      }
    }
    fix <- c(f, r, fz, rz)
    fix <- unique(fix)
    fix <- sort(fix)
    fixedeff <- paste(nam[fix], collapse = "+")
    
    nz <- as.numeric(dat$Y == 0)
    dat$nz = nz
    
    form = as.formula(paste("nz", "~", fixedeff, "+", randeff, sep = ""))
    fit <- glmmTMB::glmmTMB(formula = form,  data = dat,
                   family="binomial")
    
    fit.sum <- summary(fit)
    beta <- fixef(fit)$cond
    rv <- t(chol(vcov(fit)$cond))
    b.star <- beta + rv %*% rnorm(ncol(rv))
    fitmis <- fit
    fitmis$fit$par[names(fitmis$fit$par)=="beta"] <- b.star
    
    newdatamis = data.frame(X = x[wy, , drop = FALSE])
    colnames(newdatamis) = nam
    p <- predict(fitmis, newdata = newdatamis, type = "response",
                 na.action = na.pass, allow.new.levels = TRUE)
    im <- rbinom(n = length(p), size = 1, prob = 1 - p)
    nonzero <- which(im == 1)

    
    ran <- c(r, rc)
    ran <- unique(ran)
    ran <- sort(ran)
    randeff <- paste(nam[ran], collapse = "+")
    if (any(type == 2) | any(type == 4)) {
      randeff = paste("+", paste(nam[ran], collapse = "+"),
                      sep = "")
      if (!intercept.c) {
        randeff = paste("(0", randeff, "|", paste(nam[grp],")"),
                        sep = "")
      }
      else {
        randeff = paste("(1", randeff, "|", paste(nam[grp], ")"),
                        sep = "")
      }
    }
    else {
      if (!intercept.c) {
        randeff <- paste("(0", "|", paste(nam[grp]),")", sep = "")
      }
      else {
        randeff <- paste("(1", "|", paste(nam[grp]),")", sep = "")
      }
    }
    fix <- c(f, r, fc, rc)
    fix <- unique(fix)
    fix <- sort(fix)
    fixedeff = paste(paste(nam[-grp], collapse = "+"), sep = "")
    form = as.formula(paste("Y", "~", fixedeff, "+", randeff, sep = ""))
    fit <- glmmTMB::glmmTMB(formula = form,  data = subset(dat,
                                                  dat[, "Y"] > 0),
                   family=list(family="truncated_nbinom2",link="log"))
    fit.sum=summary(fit)
    
    fit.sum <- summary(fit)
    beta <- fixef(fit)$cond
    rv <- t(chol(vcov(fit)$cond))
    b.star <- beta + rv %*% rnorm(ncol(rv))
    fitmis <- fit
    fitmis$fit$par[names(fitmis$fit$par)=="beta"] <- b.star
    
    newdatamis = data.frame(X = x[wy, , drop = FALSE])
    colnames(newdatamis) = nam
    p <- predict(fitmis, newdata = newdatamis, type = "response",
                 na.action = na.pass, allow.new.levels = TRUE)
    
    imc = aster::rktnb(length(p), size = fit.sum$sigma, k = 0, mu = p, xpred = 1)
    
    im[nonzero] <- imc[nonzero]
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
#' @export
#' @describeIn  mice.impute.2l.hnb bootstrap regression variant; fixed interceot in zero model
mice.impute.2l.hnb.noint.zero.boot<-
  function (y, ry, x, type, intercept.c = TRUE, intercept.z = FALSE, wy = NULL, EV=TRUE)
  {
    if (is.null(wy)) 
      wy <- !ry
    
    if (sum(type == -2) > 1) {
      stop("only one class allowed!")
    }
    
    Y = y[ry]
    X = x[ry, ]
    X = data.frame(X)
    nam = colnames(X)
    grp <- which(type == -2)
    dat=data.frame(Y,X)
    GRPnam=nam[grp]
    colnames(dat)[colnames(dat)%in%nam[grp]]<-"GRP"
    cls <- sample(unique(dat$GRP),replace=TRUE)
    cls.col <- data.frame(GRP=cls)
    dat <- merge(cls.col,dat,by="GRP")
    colnames(dat)[colnames(dat)=="GRP"] <- GRPnam
    
    fc <- which(type == 3)
    rc <- which(type == 4)
    fz <- which(type == 5)
    rz <- which(type == 6)
    f <- which(type == 1)
    r <- which(type == 2)
    
    ran <- c(r, rz)
    ran <- unique(ran)
    ran <- sort(ran)
    randeff <- paste(nam[ran], collapse = "+")
    if (any(type == 2) | any(type == 6)) {
      randeff = paste("+", paste(nam[ran], collapse = "+"),
                      sep = "")
      if (!intercept.c) {
        randeff = paste("(0", randeff, "|", paste(nam[grp],")"),
                        sep = "")
      }
      else {
        randeff = paste("(1", randeff, "|", paste(nam[grp], ")"),
                        sep = "")
      }
    }
    else {
      if (!intercept.c) {
        randeff <- paste("(0", "|", paste(nam[grp]),")", sep = "")
      }
      else {
        randeff <- paste("(1", "|", paste(nam[grp]),")", sep = "")
      }
    }
    fix <- c(f, r, fz, rz)
    fix <- unique(fix)
    fix <- sort(fix)
    fixedeff <- paste(nam[fix], collapse = "+")
    
    nz <- as.numeric(dat$Y == 0)
    dat$nz = nz
    
    form = as.formula(paste("nz", "~", fixedeff, "+", randeff, sep = ""))
    fit <- glmmTMB::glmmTMB(formula = form,  data = dat,
                   family="binomial")
    
    
    newdatamis = data.frame(X = x[wy, , drop = FALSE])
    colnames(newdatamis) = nam
    p <- predict(fit, newdata = newdatamis, type = "response",
                 na.action = na.pass, allow.new.levels = TRUE)
    im <- rbinom(n = length(p), size = 1, prob = 1 - p)
    nonzero <- which(im == 1)
    
    
    ran <- c(r, rc)
    ran <- unique(ran)
    ran <- sort(ran)
    randeff <- paste(nam[ran], collapse = "+")
    if (any(type == 2) | any(type == 4)) {
      randeff = paste("+", paste(nam[ran], collapse = "+"),
                      sep = "")
      if (!intercept.c) {
        randeff = paste("(0", randeff, "|", paste(nam[grp],")"),
                        sep = "")
      }
      else {
        randeff = paste("(1", randeff, "|", paste(nam[grp], ")"),
                        sep = "")
      }
    }
    else {
      if (!intercept.c) {
        randeff <- paste("(0", "|", paste(nam[grp]),")", sep = "")
      }
      else {
        randeff <- paste("(1", "|", paste(nam[grp]),")", sep = "")
      }
    }
    fix <- c(f, r, fc, rc)
    fix <- unique(fix)
    fix <- sort(fix)
    fixedeff = paste(paste(nam[-grp], collapse = "+"), sep = "")
    form = as.formula(paste("Y", "~", fixedeff, "+", randeff, sep = ""))
    fit <- glmmTMB::glmmTMB(formula = form,  data = subset(dat,
                                                  dat[, "Y"] > 0),
                   family=list(family="truncated_nbinom2",link="log"))
    fit.sum=summary(fit)
    newdatamis = data.frame(X = x[wy, , drop = FALSE])
    colnames(newdatamis) = nam
    p <- predict(fit, newdata = newdatamis, type = "response",
                 na.action = na.pass, allow.new.levels = TRUE)
    
    imc = aster::rktnb(length(p), size = fit.sum$sigma, k = 0, mu = p, xpred = 1)
    
    im[nonzero] <- imc[nonzero]
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
#' @export
#' @describeIn  mice.impute.2l.hnb Bayesian regression variant; fixed intercept in zero model
mice.impute.2l.hnb.noint.zero <-
  function (y, ry, x, type, intercept.c = TRUE, intercept.z = FALSE, wy = NULL, EV=TRUE)
  {
    if (is.null(wy)) 
      wy <- !ry
    if (sum(type == -2) > 1) {
      stop("only one class allowed!")
    }
    
    Y = y[ry]
    X = x[ry, ]
    X = data.frame(X)
    nam = colnames(X)
    grp <- which(type == -2)
    dat=data.frame(Y,X)
    
    fc <- which(type == 3)
    rc <- which(type == 4)
    fz <- which(type == 5)
    rz <- which(type == 6)
    f <- which(type == 1)
    r <- which(type == 2)
    
    ran <- c(r, rz)
    ran <- unique(ran)
    ran <- sort(ran)
    randeff <- paste(nam[ran], collapse = "+")
    if (any(type == 2) | any(type == 6)) {
      randeff = paste("+", paste(nam[ran], collapse = "+"),
                      sep = "")
      if (!intercept.c) {
        randeff = paste("(0", randeff, "|", paste(nam[grp],")"),
                        sep = "")
      }
      else {
        randeff = paste("(1", randeff, "|", paste(nam[grp], ")"),
                        sep = "")
      }
    }
    else {
      if (!intercept.c) {
        randeff <- paste("(0", "|", paste(nam[grp]),")", sep = "")
      }
      else {
        randeff <- paste("(1", "|", paste(nam[grp]),")", sep = "")
      }
    }
    fix <- c(f, r, fz, rz)
    fix <- unique(fix)
    fix <- sort(fix)
    fixedeff <- paste(nam[fix], collapse = "+")
    
    nz <- as.numeric(dat$Y == 0)
    dat$nz = nz
    
    form = as.formula(paste("nz", "~", fixedeff, "+", randeff, sep = ""))
    fit <- glmmTMB::glmmTMB(formula = form,  data = dat,
                   family="binomial")
    
    fit.sum <- summary(fit)
    beta <- fixef(fit)$cond
    rv <- t(chol(vcov(fit)$cond))
    b.star <- beta + rv %*% rnorm(ncol(rv))
    fitmis <- fit
    fitmis$fit$par[names(fitmis$fit$par)=="beta"] <- b.star
    
    newdatamis = data.frame(X = x[wy, , drop = FALSE])
    colnames(newdatamis) = nam
    p <- predict(fitmis, newdata = newdatamis, type = "response",
                 na.action = na.pass, allow.new.levels = TRUE)
    im <- rbinom(n = length(p), size = 1, prob = 1 - p)
    nonzero <- which(im == 1)
    
    
    ran <- c(r, rc)
    ran <- unique(ran)
    ran <- sort(ran)
    randeff <- paste(nam[ran], collapse = "+")
    if (any(type == 2) | any(type == 4)) {
      randeff = paste("+", paste(nam[ran], collapse = "+"),
                      sep = "")
      if (!intercept.c) {
        randeff = paste("(0", randeff, "|", paste(nam[grp],")"),
                        sep = "")
      }
      else {
        randeff = paste("(1", randeff, "|", paste(nam[grp], ")"),
                        sep = "")
      }
    }
    else {
      if (!intercept.c) {
        randeff <- paste("(0", "|", paste(nam[grp]),")", sep = "")
      }
      else {
        randeff <- paste("(1", "|", paste(nam[grp]),")", sep = "")
      }
    }
    fix <- c(f, r, fc, rc)
    fix <- unique(fix)
    fix <- sort(fix)
    fixedeff = paste(paste(nam[-grp], collapse = "+"), sep = "")
    form = as.formula(paste("Y", "~", fixedeff, "+", randeff, sep = ""))
    fit <- glmmTMB::glmmTMB(formula = form,  data = subset(dat,
                                                  dat[, "Y"] > 0),
                   family=list(family="truncated_nbinom2",link="log"))
    fit.sum=summary(fit)
    
    fit.sum <- summary(fit)
    beta <- fixef(fit)$cond
    rv <- t(chol(vcov(fit)$cond))
    b.star <- beta + rv %*% rnorm(ncol(rv))
    fitmis <- fit
    fitmis$fit$par[names(fitmis$fit$par)=="beta"] <- b.star
    
    newdatamis = data.frame(X = x[wy, , drop = FALSE])
    colnames(newdatamis) = nam
    p <- predict(fitmis, newdata = newdatamis, type = "response",
                 na.action = na.pass, allow.new.levels = TRUE)
    
    imc = aster::rktnb(length(p), size = fit.sum$sigma, k = 0, mu = p, xpred = 1)
    
    im[nonzero] <- imc[nonzero]
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


### new 2018-04-28
#' @export
#' @describeIn mice.impute.2l.hnb Bayesian regression variant; random intercepts
mice.impute.2l.hp <-
  function (y, ry, x, type, intercept.c = TRUE, intercept.z = TRUE, wy = NULL, EV=TRUE)
  {
    if (is.null(wy)) 
      wy <- !ry
    if (sum(type == -2) > 1) {
      stop("only one class allowed!")
    }
    
    Y = y[ry]
    X = x[ry, ]
    X = data.frame(X)
    nam = colnames(X)
    grp <- which(type == -2)
    dat=data.frame(Y,X)
    
    fc <- which(type == 3)
    rc <- which(type == 4)
    fz <- which(type == 5)
    rz <- which(type == 6)
    f <- which(type == 1)
    r <- which(type == 2)
    
    ran <- c(r, rz)
    ran <- unique(ran)
    ran <- sort(ran)
    randeff <- paste(nam[ran], collapse = "+")
    if (any(type == 2) | any(type == 6)) {
      randeff = paste("+", paste(nam[ran], collapse = "+"),
                      sep = "")
      if (!intercept.c) {
        randeff = paste("(0", randeff, "|", paste(nam[grp],")"),
                        sep = "")
      }
      else {
        randeff = paste("(1", randeff, "|", paste(nam[grp], ")"),
                        sep = "")
      }
    }
    else {
      if (!intercept.c) {
        randeff <- paste("(0", "|", paste(nam[grp]),")", sep = "")
      }
      else {
        randeff <- paste("(1", "|", paste(nam[grp]),")", sep = "")
      }
    }
    fix <- c(f, r, fz, rz)
    fix <- unique(fix)
    fix <- sort(fix)
    fixedeff <- paste(nam[fix], collapse = "+")
    
    nz <- as.numeric(dat$Y == 0)
    dat$nz = nz
    
    form = as.formula(paste("nz", "~", fixedeff, "+", randeff, sep = ""))
    fit <- glmmTMB::glmmTMB(formula = form,  data = dat,
                            family="binomial")
    
    fit.sum <- summary(fit)
    beta <- fixef(fit)$cond
    rv <- t(chol(vcov(fit)$cond))
    b.star <- beta + rv %*% rnorm(ncol(rv))
    fitmis <- fit
    fitmis$fit$par[names(fitmis$fit$par)=="beta"] <- b.star
    
    newdatamis = data.frame(X = x[wy, , drop = FALSE])
    colnames(newdatamis) = nam
    p <- predict(fitmis, newdata = newdatamis, type = "response",
                 na.action = na.pass, allow.new.levels = TRUE)
    im <- rbinom(n = length(p), size = 1, prob = 1 - p)
    nonzero <- which(im == 1)
    
    ran <- c(r, rc)
    ran <- unique(ran)
    ran <- sort(ran)
    randeff <- paste(nam[ran], collapse = "+")
    if (any(type == 2) | any(type == 4)) {
      randeff = paste("+", paste(nam[ran], collapse = "+"),
                      sep = "")
      if (!intercept.c) {
        randeff = paste("(0", randeff, "|", paste(nam[grp],")"),
                        sep = "")
      }
      else {
        randeff = paste("(1", randeff, "|", paste(nam[grp], ")"),
                        sep = "")
      }
    }
    else {
      if (!intercept.c) {
        randeff <- paste("(0", "|", paste(nam[grp]),")", sep = "")
      }
      else {
        randeff <- paste("(1", "|", paste(nam[grp]),")", sep = "")
      }
    }
    fix <- c(f, r, fc, rc)
    fix <- unique(fix)
    fix <- sort(fix)
    fixedeff = paste(paste(nam[-grp], collapse = "+"), sep = "")
    form = as.formula(paste("Y", "~", fixedeff, "+", randeff, sep = ""))
    fit <- glmmTMB::glmmTMB(formula = form,  data = subset(dat,
                                                           dat[, "Y"] > 0),
                            family=list(family="truncated_poisson",link="log"))
    fit.sum=summary(fit)
    
    fit.sum <- summary(fit)
    beta <- fixef(fit)$cond
    rv <- t(chol(vcov(fit)$cond))
    b.star <- beta + rv %*% rnorm(ncol(rv))
    fitmis <- fit
    fitmis$fit$par[names(fitmis$fit$par)=="beta"] <- b.star
    
    newdatamis = data.frame(X = x[wy, , drop = FALSE])
    colnames(newdatamis) = nam
    p <- predict(fitmis, newdata = newdatamis, type = "response",
                 na.action = na.pass, allow.new.levels = TRUE)
    
    imc = aster::rktp(n=length(p), k = 0, mu = p, xpred = 1)
    
    im[nonzero] <- imc[nonzero]
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

#' @export
#' @describeIn  mice.impute.2l.hnb bootstrap regression variant; random intercepts
mice.impute.2l.hp.boot<-
  function (y, ry, x, type, intercept.c = TRUE, intercept.z = TRUE, wy = NULL, EV=TRUE)
  {
    if (is.null(wy)) 
      wy <- !ry
    
    if (sum(type == -2) > 1) {
      stop("only one class allowed!")
    }
    
    Y = y[ry]
    X = x[ry, ]
    X = data.frame(X)
    nam = colnames(X)
    grp <- which(type == -2)
    dat=data.frame(Y,X)
    GRPnam=nam[grp]
    colnames(dat)[colnames(dat)%in%nam[grp]]<-"GRP"
    cls <- sample(unique(dat$GRP),replace=TRUE)
    cls.col <- data.frame(GRP=cls)
    dat <- merge(cls.col,dat,by="GRP")
    colnames(dat)[colnames(dat)=="GRP"] <- GRPnam
    
    fc <- which(type == 3)
    rc <- which(type == 4)
    fz <- which(type == 5)
    rz <- which(type == 6)
    f <- which(type == 1)
    r <- which(type == 2)
    
    ran <- c(r, rz)
    ran <- unique(ran)
    ran <- sort(ran)
    randeff <- paste(nam[ran], collapse = "+")
    if (any(type == 2) | any(type == 6)) {
      randeff = paste("+", paste(nam[ran], collapse = "+"),
                      sep = "")
      if (!intercept.c) {
        randeff = paste("(0", randeff, "|", paste(nam[grp],")"),
                        sep = "")
      }
      else {
        randeff = paste("(1", randeff, "|", paste(nam[grp], ")"),
                        sep = "")
      }
    }
    else {
      if (!intercept.c) {
        randeff <- paste("(0", "|", paste(nam[grp]),")", sep = "")
      }
      else {
        randeff <- paste("(1", "|", paste(nam[grp]),")", sep = "")
      }
    }
    fix <- c(f, r, fz, rz)
    fix <- unique(fix)
    fix <- sort(fix)
    fixedeff <- paste(nam[fix], collapse = "+")
    
    nz <- as.numeric(dat$Y == 0)
    dat$nz = nz
    
    form = as.formula(paste("nz", "~", fixedeff, "+", randeff, sep = ""))
    fit <- glmmTMB::glmmTMB(formula = form,  data = dat,
                            family="binomial")
    
    
    newdatamis = data.frame(X = x[wy, , drop = FALSE])
    colnames(newdatamis) = nam
    p <- predict(fit, newdata = newdatamis, type = "response",
                 na.action = na.pass, allow.new.levels = TRUE)
    im <- rbinom(n = length(p), size = 1, prob = 1 - p)
    nonzero <- which(im == 1)
    
    ran <- c(r, rc)
    ran <- unique(ran)
    ran <- sort(ran)
    randeff <- paste(nam[ran], collapse = "+")
    if (any(type == 2) | any(type == 4)) {
      randeff = paste("+", paste(nam[ran], collapse = "+"),
                      sep = "")
      if (!intercept.c) {
        randeff = paste("(0", randeff, "|", paste(nam[grp],")"),
                        sep = "")
      }
      else {
        randeff = paste("(1", randeff, "|", paste(nam[grp], ")"),
                        sep = "")
      }
    }
    else {
      if (!intercept.c) {
        randeff <- paste("(0", "|", paste(nam[grp]),")", sep = "")
      }
      else {
        randeff <- paste("(1", "|", paste(nam[grp]),")", sep = "")
      }
    }
    fix <- c(f, r, fc, rc)
    fix <- unique(fix)
    fix <- sort(fix)
    fixedeff = paste(paste(nam[-grp], collapse = "+"), sep = "")
    form = as.formula(paste("Y", "~", fixedeff, "+", randeff, sep = ""))
    fit <- glmmTMB::glmmTMB(formula = form,  data = subset(dat,
                                                           dat[, "Y"] > 0),
                            family=list(family="truncated_poisson",link="log"))
    fit.sum=summary(fit)
    newdatamis = data.frame(X = x[wy, , drop = FALSE])
    colnames(newdatamis) = nam
    p <- predict(fit, newdata = newdatamis, type = "response",
                 na.action = na.pass, allow.new.levels = TRUE)
    
    imc = aster::rktp(n=length(p),  k = 0, mu = p, xpred = 1)
    
    im[nonzero] <- imc[nonzero]
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

#' @export
#' @describeIn  mice.impute.2l.hnb bootstrap regression variant; fixed intercepts
mice.impute.2l.hp.noint.both.boot<-
  function (y, ry, x, type, intercept.c = FALSE, intercept.z = FALSE, wy = NULL, EV=TRUE)
  {
    if (is.null(wy)) 
      wy <- !ry
    
    if (sum(type == -2) > 1) {
      stop("only one class allowed!")
    }
    
    Y = y[ry]
    X = x[ry, ]
    X = data.frame(X)
    nam = colnames(X)
    grp <- which(type == -2)
    dat=data.frame(Y,X)
    GRPnam=nam[grp]
    colnames(dat)[colnames(dat)%in%nam[grp]]<-"GRP"
    cls <- sample(unique(dat$GRP),replace=TRUE)
    cls.col <- data.frame(GRP=cls)
    dat <- merge(cls.col,dat,by="GRP")
    colnames(dat)[colnames(dat)=="GRP"] <- GRPnam
    
    fc <- which(type == 3)
    rc <- which(type == 4)
    fz <- which(type == 5)
    rz <- which(type == 6)
    f <- which(type == 1)
    r <- which(type == 2)
    
    ran <- c(r, rz)
    ran <- unique(ran)
    ran <- sort(ran)
    randeff <- paste(nam[ran], collapse = "+")
    if (any(type == 2) | any(type == 6)) {
      randeff = paste("+", paste(nam[ran], collapse = "+"),
                      sep = "")
      if (!intercept.c) {
        randeff = paste("(0", randeff, "|", paste(nam[grp],")"),
                        sep = "")
      }
      else {
        randeff = paste("(1", randeff, "|", paste(nam[grp], ")"),
                        sep = "")
      }
    }
    else {
      if (!intercept.c) {
        randeff <- paste("(0", "|", paste(nam[grp]),")", sep = "")
      }
      else {
        randeff <- paste("(1", "|", paste(nam[grp]),")", sep = "")
      }
    }
    fix <- c(f, r, fz, rz)
    fix <- unique(fix)
    fix <- sort(fix)
    fixedeff <- paste(nam[fix], collapse = "+")
    
    nz <- as.numeric(dat$Y == 0)
    dat$nz = nz
    
    form = as.formula(paste("nz", "~", fixedeff, "+", randeff, sep = ""))
    fit <- glmmTMB::glmmTMB(formula = form,  data = dat,
                            family="binomial")
    
    
    newdatamis = data.frame(X = x[wy, , drop = FALSE])
    colnames(newdatamis) = nam
    p <- predict(fit, newdata = newdatamis, type = "response",
                 na.action = na.pass, allow.new.levels = TRUE)
    im <- rbinom(n = length(p), size = 1, prob = 1 - p)
    nonzero <- which(im == 1)
    
    ran <- c(r, rc)
    ran <- unique(ran)
    ran <- sort(ran)
    randeff <- paste(nam[ran], collapse = "+")
    if (any(type == 2) | any(type == 4)) {
      randeff = paste("+", paste(nam[ran], collapse = "+"),
                      sep = "")
      if (!intercept.c) {
        randeff = paste("(0", randeff, "|", paste(nam[grp],")"),
                        sep = "")
      }
      else {
        randeff = paste("(1", randeff, "|", paste(nam[grp], ")"),
                        sep = "")
      }
    }
    else {
      if (!intercept.c) {
        randeff <- paste("(0", "|", paste(nam[grp]),")", sep = "")
      }
      else {
        randeff <- paste("(1", "|", paste(nam[grp]),")", sep = "")
      }
    }
    fix <- c(f, r, fc, rc)
    fix <- unique(fix)
    fix <- sort(fix)
    fixedeff = paste(paste(nam[-grp], collapse = "+"), sep = "")
    form = as.formula(paste("Y", "~", fixedeff, "+", randeff, sep = ""))
    fit <- glmmTMB::glmmTMB(formula = form,  data = subset(dat,
                                                           dat[, "Y"] > 0),
                            family=list(family="truncated_poisson",link="log"))
    fit.sum=summary(fit)
    newdatamis = data.frame(X = x[wy, , drop = FALSE])
    colnames(newdatamis) = nam
    p <- predict(fit, newdata = newdatamis, type = "response",
                 na.action = na.pass, allow.new.levels = TRUE)
    
    imc = aster::rktp(n=length(p), k = 0, mu = p, xpred = 1)
    
    im[nonzero] <- imc[nonzero]
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

#' @export
#' @describeIn  mice.impute.2l.hnb Bayesian regression variant; fixed intercepts
mice.impute.2l.hp.noint.both <-
  function (y, ry, x, type, intercept.c = FALSE, intercept.z = FALSE, wy = NULL, EV=TRUE)
  {
    if (is.null(wy)) 
      wy <- !ry
    if (sum(type == -2) > 1) {
      stop("only one class allowed!")
    }
    
    Y = y[ry]
    X = x[ry, ]
    X = data.frame(X)
    nam = colnames(X)
    grp <- which(type == -2)
    dat=data.frame(Y,X)
    
    fc <- which(type == 3)
    rc <- which(type == 4)
    fz <- which(type == 5)
    rz <- which(type == 6)
    f <- which(type == 1)
    r <- which(type == 2)
    
    ran <- c(r, rz)
    ran <- unique(ran)
    ran <- sort(ran)
    randeff <- paste(nam[ran], collapse = "+")
    if (any(type == 2) | any(type == 6)) {
      randeff = paste("+", paste(nam[ran], collapse = "+"),
                      sep = "")
      if (!intercept.c) {
        randeff = paste("(0", randeff, "|", paste(nam[grp],")"),
                        sep = "")
      }
      else {
        randeff = paste("(1", randeff, "|", paste(nam[grp], ")"),
                        sep = "")
      }
    }
    else {
      if (!intercept.c) {
        randeff <- paste("(0", "|", paste(nam[grp]),")", sep = "")
      }
      else {
        randeff <- paste("(1", "|", paste(nam[grp]),")", sep = "")
      }
    }
    fix <- c(f, r, fz, rz)
    fix <- unique(fix)
    fix <- sort(fix)
    fixedeff <- paste(nam[fix], collapse = "+")
    
    nz <- as.numeric(dat$Y == 0)
    dat$nz = nz
    
    form = as.formula(paste("nz", "~", fixedeff, "+", randeff, sep = ""))
    fit <- glmmTMB::glmmTMB(formula = form,  data = dat,
                            family="binomial")
    
    fit.sum <- summary(fit)
    beta <- fixef(fit)$cond
    rv <- t(chol(vcov(fit)$cond))
    b.star <- beta + rv %*% rnorm(ncol(rv))
    fitmis <- fit
    fitmis$fit$par[names(fitmis$fit$par)=="beta"] <- b.star
    
    newdatamis = data.frame(X = x[wy, , drop = FALSE])
    colnames(newdatamis) = nam
    p <- predict(fitmis, newdata = newdatamis, type = "response",
                 na.action = na.pass, allow.new.levels = TRUE)
    im <- rbinom(n = length(p), size = 1, prob = 1 - p)
    nonzero <- which(im == 1)
    
    ran <- c(r, rc)
    ran <- unique(ran)
    ran <- sort(ran)
    randeff <- paste(nam[ran], collapse = "+")
    if (any(type == 2) | any(type == 4)) {
      randeff = paste("+", paste(nam[ran], collapse = "+"),
                      sep = "")
      if (!intercept.c) {
        randeff = paste("(0", randeff, "|", paste(nam[grp],")"),
                        sep = "")
      }
      else {
        randeff = paste("(1", randeff, "|", paste(nam[grp], ")"),
                        sep = "")
      }
    }
    else {
      if (!intercept.c) {
        randeff <- paste("(0", "|", paste(nam[grp]),")", sep = "")
      }
      else {
        randeff <- paste("(1", "|", paste(nam[grp]),")", sep = "")
      }
    }
    fix <- c(f, r, fc, rc)
    fix <- unique(fix)
    fix <- sort(fix)
    fixedeff = paste(paste(nam[-grp], collapse = "+"), sep = "")
    form = as.formula(paste("Y", "~", fixedeff, "+", randeff, sep = ""))
    fit <- glmmTMB::glmmTMB(formula = form,  data = subset(dat,
                                                           dat[, "Y"] > 0),
                            family=list(family="truncated_poisson",link="log"))
    fit.sum=summary(fit)
    
    fit.sum <- summary(fit)
    beta <- fixef(fit)$cond
    rv <- t(chol(vcov(fit)$cond))
    b.star <- beta + rv %*% rnorm(ncol(rv))
    fitmis <- fit
    fitmis$fit$par[names(fitmis$fit$par)=="beta"] <- b.star
    
    newdatamis = data.frame(X = x[wy, , drop = FALSE])
    colnames(newdatamis) = nam
    p <- predict(fitmis, newdata = newdatamis, type = "response",
                 na.action = na.pass, allow.new.levels = TRUE)
    
    imc = aster::rktp(n=length(p), k = 0, mu = p, xpred = 1)
    
    im[nonzero] <- imc[nonzero]
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

#' @export
#' @describeIn  mice.impute.2l.hnb bootstrap regression variant; fixed intercept in count model
mice.impute.2l.hp.noint.count.boot<-
  function (y, ry, x, type, intercept.c = FALSE, intercept.z = TRUE, wy = NULL, EV=TRUE)
  {
    if (is.null(wy)) 
      wy <- !ry
    
    if (sum(type == -2) > 1) {
      stop("only one class allowed!")
    }
    
    Y = y[ry]
    X = x[ry, ]
    X = data.frame(X)
    nam = colnames(X)
    grp <- which(type == -2)
    dat=data.frame(Y,X)
    GRPnam=nam[grp]
    colnames(dat)[colnames(dat)%in%nam[grp]]<-"GRP"
    cls <- sample(unique(dat$GRP),replace=TRUE)
    cls.col <- data.frame(GRP=cls)
    dat <- merge(cls.col,dat,by="GRP")
    colnames(dat)[colnames(dat)=="GRP"] <- GRPnam
    
    fc <- which(type == 3)
    rc <- which(type == 4)
    fz <- which(type == 5)
    rz <- which(type == 6)
    f <- which(type == 1)
    r <- which(type == 2)
    
    ran <- c(r, rz)
    ran <- unique(ran)
    ran <- sort(ran)
    randeff <- paste(nam[ran], collapse = "+")
    if (any(type == 2) | any(type == 6)) {
      randeff = paste("+", paste(nam[ran], collapse = "+"),
                      sep = "")
      if (!intercept.c) {
        randeff = paste("(0", randeff, "|", paste(nam[grp],")"),
                        sep = "")
      }
      else {
        randeff = paste("(1", randeff, "|", paste(nam[grp], ")"),
                        sep = "")
      }
    }
    else {
      if (!intercept.c) {
        randeff <- paste("(0", "|", paste(nam[grp]),")", sep = "")
      }
      else {
        randeff <- paste("(1", "|", paste(nam[grp]),")", sep = "")
      }
    }
    fix <- c(f, r, fz, rz)
    fix <- unique(fix)
    fix <- sort(fix)
    fixedeff <- paste(nam[fix], collapse = "+")
    
    nz <- as.numeric(dat$Y == 0)
    dat$nz = nz
    
    form = as.formula(paste("nz", "~", fixedeff, "+", randeff, sep = ""))
    fit <- glmmTMB::glmmTMB(formula = form,  data = dat,
                            family="binomial")
    
    
    newdatamis = data.frame(X = x[wy, , drop = FALSE])
    colnames(newdatamis) = nam
    p <- predict(fit, newdata = newdatamis, type = "response",
                 na.action = na.pass, allow.new.levels = TRUE)
    im <- rbinom(n = length(p), size = 1, prob = 1 - p)
    nonzero <- which(im == 1)
    
    ran <- c(r, rc)
    ran <- unique(ran)
    ran <- sort(ran)
    randeff <- paste(nam[ran], collapse = "+")
    if (any(type == 2) | any(type == 4)) {
      randeff = paste("+", paste(nam[ran], collapse = "+"),
                      sep = "")
      if (!intercept.c) {
        randeff = paste("(0", randeff, "|", paste(nam[grp],")"),
                        sep = "")
      }
      else {
        randeff = paste("(1", randeff, "|", paste(nam[grp], ")"),
                        sep = "")
      }
    }
    else {
      if (!intercept.c) {
        randeff <- paste("(0", "|", paste(nam[grp]),")", sep = "")
      }
      else {
        randeff <- paste("(1", "|", paste(nam[grp]),")", sep = "")
      }
    }
    fix <- c(f, r, fc, rc)
    fix <- unique(fix)
    fix <- sort(fix)
    fixedeff = paste(paste(nam[-grp], collapse = "+"), sep = "")
    form = as.formula(paste("Y", "~", fixedeff, "+", randeff, sep = ""))
    fit <- glmmTMB::glmmTMB(formula = form,  data = subset(dat,
                                                           dat[, "Y"] > 0),
                            family=list(family="truncated_poisson",link="log"))
    fit.sum=summary(fit)
    newdatamis = data.frame(X = x[wy, , drop = FALSE])
    colnames(newdatamis) = nam
    p <- predict(fit, newdata = newdatamis, type = "response",
                 na.action = na.pass, allow.new.levels = TRUE)
    
    imc = aster::rktp(n=length(p),  k = 0, mu = p, xpred = 1)
    
    im[nonzero] <- imc[nonzero]
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

#' @export
#' @describeIn  mice.impute.2l.hnb Bayesian regression variant; fixed intercept in count model
mice.impute.2l.hp.noint.count <-
  function (y, ry, x, type, intercept.c = FALSE, intercept.z = TRUE, wy = NULL, EV=TRUE)
  {
    if (is.null(wy)) 
      wy <- !ry
    if (sum(type == -2) > 1) {
      stop("only one class allowed!")
    }
    
    Y = y[ry]
    X = x[ry, ]
    X = data.frame(X)
    nam = colnames(X)
    grp <- which(type == -2)
    dat=data.frame(Y,X)
    
    fc <- which(type == 3)
    rc <- which(type == 4)
    fz <- which(type == 5)
    rz <- which(type == 6)
    f <- which(type == 1)
    r <- which(type == 2)
    
    ran <- c(r, rz)
    ran <- unique(ran)
    ran <- sort(ran)
    randeff <- paste(nam[ran], collapse = "+")
    if (any(type == 2) | any(type == 6)) {
      randeff = paste("+", paste(nam[ran], collapse = "+"),
                      sep = "")
      if (!intercept.c) {
        randeff = paste("(0", randeff, "|", paste(nam[grp],")"),
                        sep = "")
      }
      else {
        randeff = paste("(1", randeff, "|", paste(nam[grp], ")"),
                        sep = "")
      }
    }
    else {
      if (!intercept.c) {
        randeff <- paste("(0", "|", paste(nam[grp]),")", sep = "")
      }
      else {
        randeff <- paste("(1", "|", paste(nam[grp]),")", sep = "")
      }
    }
    fix <- c(f, r, fz, rz)
    fix <- unique(fix)
    fix <- sort(fix)
    fixedeff <- paste(nam[fix], collapse = "+")
    
    nz <- as.numeric(dat$Y == 0)
    dat$nz = nz
    
    form = as.formula(paste("nz", "~", fixedeff, "+", randeff, sep = ""))
    fit <- glmmTMB::glmmTMB(formula = form,  data = dat,
                            family="binomial")
    
    fit.sum <- summary(fit)
    beta <- fixef(fit)$cond
    rv <- t(chol(vcov(fit)$cond))
    b.star <- beta + rv %*% rnorm(ncol(rv))
    fitmis <- fit
    fitmis$fit$par[names(fitmis$fit$par)=="beta"] <- b.star
    
    newdatamis = data.frame(X = x[wy, , drop = FALSE])
    colnames(newdatamis) = nam
    p <- predict(fitmis, newdata = newdatamis, type = "response",
                 na.action = na.pass, allow.new.levels = TRUE)
    im <- rbinom(n = length(p), size = 1, prob = 1 - p)
    nonzero <- which(im == 1)
    
    
    ran <- c(r, rc)
    ran <- unique(ran)
    ran <- sort(ran)
    randeff <- paste(nam[ran], collapse = "+")
    if (any(type == 2) | any(type == 4)) {
      randeff = paste("+", paste(nam[ran], collapse = "+"),
                      sep = "")
      if (!intercept.c) {
        randeff = paste("(0", randeff, "|", paste(nam[grp],")"),
                        sep = "")
      }
      else {
        randeff = paste("(1", randeff, "|", paste(nam[grp], ")"),
                        sep = "")
      }
    }
    else {
      if (!intercept.c) {
        randeff <- paste("(0", "|", paste(nam[grp]),")", sep = "")
      }
      else {
        randeff <- paste("(1", "|", paste(nam[grp]),")", sep = "")
      }
    }
    fix <- c(f, r, fc, rc)
    fix <- unique(fix)
    fix <- sort(fix)
    fixedeff = paste(paste(nam[-grp], collapse = "+"), sep = "")
    form = as.formula(paste("Y", "~", fixedeff, "+", randeff, sep = ""))
    fit <- glmmTMB::glmmTMB(formula = form,  data = subset(dat,
                                                           dat[, "Y"] > 0),
                            family=list(family="truncated_poisson",link="log"))
    fit.sum=summary(fit)
    
    fit.sum <- summary(fit)
    beta <- fixef(fit)$cond
    rv <- t(chol(vcov(fit)$cond))
    b.star <- beta + rv %*% rnorm(ncol(rv))
    fitmis <- fit
    fitmis$fit$par[names(fitmis$fit$par)=="beta"] <- b.star
    
    newdatamis = data.frame(X = x[wy, , drop = FALSE])
    colnames(newdatamis) = nam
    p <- predict(fitmis, newdata = newdatamis, type = "response",
                 na.action = na.pass, allow.new.levels = TRUE)
    
    imc = aster::rktp(length(p),  k = 0, mu = p, xpred = 1)
    
    im[nonzero] <- imc[nonzero]
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

#' @export
#' @describeIn  mice.impute.2l.hnb bootstrap regression variant; fixed interceot in zero model
mice.impute.2l.hp.noint.zero.boot<-
  function (y, ry, x, type, intercept.c = TRUE, intercept.z = FALSE, wy = NULL, EV=TRUE)
  {
    if (is.null(wy)) 
      wy <- !ry
    
    if (sum(type == -2) > 1) {
      stop("only one class allowed!")
    }
    
    Y = y[ry]
    X = x[ry, ]
    X = data.frame(X)
    nam = colnames(X)
    grp <- which(type == -2)
    dat=data.frame(Y,X)
    GRPnam=nam[grp]
    colnames(dat)[colnames(dat)%in%nam[grp]]<-"GRP"
    cls <- sample(unique(dat$GRP),replace=TRUE)
    cls.col <- data.frame(GRP=cls)
    dat <- merge(cls.col,dat,by="GRP")
    colnames(dat)[colnames(dat)=="GRP"] <- GRPnam
    
    fc <- which(type == 3)
    rc <- which(type == 4)
    fz <- which(type == 5)
    rz <- which(type == 6)
    f <- which(type == 1)
    r <- which(type == 2)
    
    ran <- c(r, rz)
    ran <- unique(ran)
    ran <- sort(ran)
    randeff <- paste(nam[ran], collapse = "+")
    if (any(type == 2) | any(type == 6)) {
      randeff = paste("+", paste(nam[ran], collapse = "+"),
                      sep = "")
      if (!intercept.c) {
        randeff = paste("(0", randeff, "|", paste(nam[grp],")"),
                        sep = "")
      }
      else {
        randeff = paste("(1", randeff, "|", paste(nam[grp], ")"),
                        sep = "")
      }
    }
    else {
      if (!intercept.c) {
        randeff <- paste("(0", "|", paste(nam[grp]),")", sep = "")
      }
      else {
        randeff <- paste("(1", "|", paste(nam[grp]),")", sep = "")
      }
    }
    fix <- c(f, r, fz, rz)
    fix <- unique(fix)
    fix <- sort(fix)
    fixedeff <- paste(nam[fix], collapse = "+")
    
    nz <- as.numeric(dat$Y == 0)
    dat$nz = nz
    
    form = as.formula(paste("nz", "~", fixedeff, "+", randeff, sep = ""))
    fit <- glmmTMB::glmmTMB(formula = form,  data = dat,
                            family="binomial")
    
    
    newdatamis = data.frame(X = x[wy, , drop = FALSE])
    colnames(newdatamis) = nam
    p <- predict(fit, newdata = newdatamis, type = "response",
                 na.action = na.pass, allow.new.levels = TRUE)
    im <- rbinom(n = length(p), size = 1, prob = 1 - p)
    nonzero <- which(im == 1)
    
    
    ran <- c(r, rc)
    ran <- unique(ran)
    ran <- sort(ran)
    randeff <- paste(nam[ran], collapse = "+")
    if (any(type == 2) | any(type == 4)) {
      randeff = paste("+", paste(nam[ran], collapse = "+"),
                      sep = "")
      if (!intercept.c) {
        randeff = paste("(0", randeff, "|", paste(nam[grp],")"),
                        sep = "")
      }
      else {
        randeff = paste("(1", randeff, "|", paste(nam[grp], ")"),
                        sep = "")
      }
    }
    else {
      if (!intercept.c) {
        randeff <- paste("(0", "|", paste(nam[grp]),")", sep = "")
      }
      else {
        randeff <- paste("(1", "|", paste(nam[grp]),")", sep = "")
      }
    }
    fix <- c(f, r, fc, rc)
    fix <- unique(fix)
    fix <- sort(fix)
    fixedeff = paste(paste(nam[-grp], collapse = "+"), sep = "")
    form = as.formula(paste("Y", "~", fixedeff, "+", randeff, sep = ""))
    fit <- glmmTMB::glmmTMB(formula = form,  data = subset(dat,
                                                           dat[, "Y"] > 0),
                            family=list(family="truncated_poisson",link="log"))
    fit.sum=summary(fit)
    newdatamis = data.frame(X = x[wy, , drop = FALSE])
    colnames(newdatamis) = nam
    p <- predict(fit, newdata = newdatamis, type = "response",
                 na.action = na.pass, allow.new.levels = TRUE)
    
    imc = aster::rktp(length(p),  k = 0, mu = p, xpred = 1)
    
    im[nonzero] <- imc[nonzero]
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

#' @export
#' @describeIn  mice.impute.2l.hnb Bayesian regression variant; fixed intercept in zero model
mice.impute.2l.hp.noint.zero <-
  function (y, ry, x, type, intercept.c = TRUE, intercept.z = FALSE, wy = NULL, EV=TRUE)
  {
    if (is.null(wy)) 
      wy <- !ry
    if (sum(type == -2) > 1) {
      stop("only one class allowed!")
    }
    
    Y = y[ry]
    X = x[ry, ]
    X = data.frame(X)
    nam = colnames(X)
    grp <- which(type == -2)
    dat=data.frame(Y,X)
    
    fc <- which(type == 3)
    rc <- which(type == 4)
    fz <- which(type == 5)
    rz <- which(type == 6)
    f <- which(type == 1)
    r <- which(type == 2)
    
    ran <- c(r, rz)
    ran <- unique(ran)
    ran <- sort(ran)
    randeff <- paste(nam[ran], collapse = "+")
    if (any(type == 2) | any(type == 6)) {
      randeff = paste("+", paste(nam[ran], collapse = "+"),
                      sep = "")
      if (!intercept.c) {
        randeff = paste("(0", randeff, "|", paste(nam[grp],")"),
                        sep = "")
      }
      else {
        randeff = paste("(1", randeff, "|", paste(nam[grp], ")"),
                        sep = "")
      }
    }
    else {
      if (!intercept.c) {
        randeff <- paste("(0", "|", paste(nam[grp]),")", sep = "")
      }
      else {
        randeff <- paste("(1", "|", paste(nam[grp]),")", sep = "")
      }
    }
    fix <- c(f, r, fz, rz)
    fix <- unique(fix)
    fix <- sort(fix)
    fixedeff <- paste(nam[fix], collapse = "+")
    
    nz <- as.numeric(dat$Y == 0)
    dat$nz = nz
    
    form = as.formula(paste("nz", "~", fixedeff, "+", randeff, sep = ""))
    fit <- glmmTMB::glmmTMB(formula = form,  data = dat,
                            family="binomial")
    
    fit.sum <- summary(fit)
    beta <- fixef(fit)$cond
    rv <- t(chol(vcov(fit)$cond))
    b.star <- beta + rv %*% rnorm(ncol(rv))
    fitmis <- fit
    fitmis$fit$par[names(fitmis$fit$par)=="beta"] <- b.star
    
    newdatamis = data.frame(X = x[wy, , drop = FALSE])
    colnames(newdatamis) = nam
    p <- predict(fitmis, newdata = newdatamis, type = "response",
                 na.action = na.pass, allow.new.levels = TRUE)
    im <- rbinom(n = length(p), size = 1, prob = 1 - p)
    nonzero <- which(im == 1)
    
    
    ran <- c(r, rc)
    ran <- unique(ran)
    ran <- sort(ran)
    randeff <- paste(nam[ran], collapse = "+")
    if (any(type == 2) | any(type == 4)) {
      randeff = paste("+", paste(nam[ran], collapse = "+"),
                      sep = "")
      if (!intercept.c) {
        randeff = paste("(0", randeff, "|", paste(nam[grp],")"),
                        sep = "")
      }
      else {
        randeff = paste("(1", randeff, "|", paste(nam[grp], ")"),
                        sep = "")
      }
    }
    else {
      if (!intercept.c) {
        randeff <- paste("(0", "|", paste(nam[grp]),")", sep = "")
      }
      else {
        randeff <- paste("(1", "|", paste(nam[grp]),")", sep = "")
      }
    }
    fix <- c(f, r, fc, rc)
    fix <- unique(fix)
    fix <- sort(fix)
    fixedeff = paste(paste(nam[-grp], collapse = "+"), sep = "")
    form = as.formula(paste("Y", "~", fixedeff, "+", randeff, sep = ""))
    fit <- glmmTMB::glmmTMB(formula = form,  data = subset(dat,
                                                           dat[, "Y"] > 0),
                            family=list(family="truncated_poisson",link="log"))
    fit.sum=summary(fit)
    
    fit.sum <- summary(fit)
    beta <- fixef(fit)$cond
    rv <- t(chol(vcov(fit)$cond))
    b.star <- beta + rv %*% rnorm(ncol(rv))
    fitmis <- fit
    fitmis$fit$par[names(fitmis$fit$par)=="beta"] <- b.star
    
    newdatamis = data.frame(X = x[wy, , drop = FALSE])
    colnames(newdatamis) = nam
    p <- predict(fitmis, newdata = newdatamis, type = "response",
                 na.action = na.pass, allow.new.levels = TRUE)
    
    imc = aster::rktp(length(p), k = 0, mu = p, xpred = 1)
    
    im[nonzero] <- imc[nonzero]
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