#' Two-level predictive mean matching
#'
#' The function imputes an incomplete variable based on a normal linear mixed effects model. The model is estimated using function \code{glmmPQL()} from package \pkg{MASS}. Matching is done by \code{.pmm.match} from package \pkg{mice}. 
#'
#'  Model specification / allowed entries in \code{mice}'s \code{predictorMatrix}:
#'    \itemize{
#'      \item 0 = variable not included in imputation model
#'      \item 1 = fixed effect
#'      \item 2 = fixex and random effect
#'      \item -2 = class variable
#'    }
#' @param y Numeric vector with incomplete data
#' @param ry Response pattern of \code{y} (\code{TRUE}=observed, \code{FALSE}=missing)
#' @param x matrix with \code{length(y)} rows containing complete covariates
#' @param type vector of length \code{ncol(x)} determining the imputation model; \code{type} is automatically extracted from the \code{predictorMatrix} argument of \code{mice()}.
#' @param wy Logical vector of length \code{length(y)}. A \code{TRUE} value indicates locations in \code{y} for which imputations are created. Default is \code{!ry}
#' @param intercept Logical. shall the intercept be included as a fixed and random effect?. \code{TRUE} = yes; \code{FALSE} = fixed effect only.
#' @param donors The size of the donor pool; default is 5.
#' @param ... additional arguments passed down from the main mice call 
#' @return vector with imputations
#' @author Kristian Kleinke
#' @references
#' \itemize{
#' \item Kleinke, K. (2016, September). Multiple Imputation of Multilevel Data by "Two-Level Predictive Mean Matching". Paper presented at the 50th Congress of the German Psychological Society (DGPs), Leipzig, Germany.
#'} 
#' @import MASS mice            
#' @export
mice.impute.2l.pmm <-
  function (y, ry, x, type, intercept = TRUE, donors=5, wy = NULL, ...){
    if (is.null(wy)) 
      wy <- !ry
    Y <- y[ry]
    X <- x[ry, ]
    nam <- colnames(X) 
    if (sum(type == -2) > 1) {
      stop("only one class allowed!")
    }
    grp <- which(type == -2)
    fixedeff <- paste("1+", paste(nam[-grp], collapse = "+"), 
                      sep = "")
    if (any(type == 2)) {
      ran <- which(type == 2)
      randeff <- paste("+", paste(nam[ran], collapse = "+"), 
                       sep = "")
      if (!intercept) {
        randeff <- paste("~0", randeff, "|", paste(nam[grp]), 
                         sep = "")
      }
      else {
        randeff <- paste("~1", randeff, "|", paste(nam[grp]), 
                         sep = "")
      }
    }
    else {
      if (!intercept) {
        randeff <- paste("~0", "|", paste(nam[grp]), sep = "")
      }
      else {
        randeff <- paste("~1", "|", paste(nam[grp]), sep = "")
      }
    }
    randeff <- as.formula(randeff)
    fixedeff <- as.formula(paste("Y", "~", fixedeff, sep = ""))
    dat <- data.frame(Y, X)
    fit <- MASS::glmmPQL(fixed = fixedeff, data = dat, random = randeff, 
                   family = "gaussian", control = list(opt =
                                                         "optim"), na.action = na.omit)
    fit.sum <- summary(fit)
    beta <- fit$coefficients$fixed
    rv <- t(chol(vcov(fit)))
    b.star <- beta + rv %*% rnorm(ncol(rv))
    fitmis <- fit
    fitmis$coefficients$fixed <- b.star
    newdatamis <- data.frame(X = x[!ry, ])
    newdataobs <- data.frame(X = x[ry, ])
    colnames(newdatamis) <- nam
    colnames(newdataobs) <- nam
    yhatmis <- predict(fitmis, newdata = newdatamis, type = "response", 
                       na.action = na.pass)
    yhatobs <- predict(fit, newdata = newdataobs, type = "response", 
                       na.action = na.pass)
    return(apply(as.array(yhatmis), 1, .pmm.match, yhat = yhatobs, 
                 y = y[ry], donors = donors, ...))
  }	

.pmm.match <-
function (z, yhat = yhat, y = y, donors = 5, ...) 
{
  d <- abs(yhat - z)
  f <- d > 0
  a1 <- ifelse(any(f), min(d[f]), 1)
  d <- d + runif(length(d), 0, a1/10^10)
  if (donors == 1) 
    return(y[which.min(d)])
  donors <- min(donors, length(d))
  donors <- max(donors, 1)
  ds <- sort.int(d, partial = donors)
  m <- sample(y[d <= ds[donors]], 1)
  return(m)
}
