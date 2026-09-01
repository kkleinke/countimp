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
#' @import MASS
#' @export
mice.impute.2l.pmm <-
  function (y, ry, x, type, intercept = TRUE, donors=5, wy = NULL, ...){
    if (is.null(wy)) 
      wy <- !ry
    Y <- y[ry]
    X <- x[ry, ]
    nam <- colnames(X) 
    ## One grouping level only, and here that is a property of the backend:
    ## this method fits through MASS::glmmPQL/nlme, which takes a single
    ## grouping factor. The glmmTMB-based two-level families accept several
    ## since 3.0.0 (a second -2 makes a three-level model); pmm does not.
    if (sum(type == -2) > 1) {
      stop("countimp: mice.impute.2l.pmm handles one grouping level, but `type` ",
           "codes ", sum(type == -2), ".\n",
           "  This method fits through glmmPQL (nlme), which takes a single ",
           "grouping factor.\n  Use a model-based two-level method (2l.poisson, ",
           "2l.nb, 2l.zip, ...) for more\n  than one level.", call. = FALSE)
    }
    grp <- which(type == -2)
    if (length(grp) == 0L)
      stop("2l.pmm needs a class variable: code one predictor as -2 in `type`.",
           call. = FALSE)
    fixedeff <- paste("1+", paste(nam[-grp], collapse = "+"), 
                      sep = "")
    ## The four inline branches this replaces built "~0 | grp" when
    ## intercept = FALSE and no predictor was coded 2 -- an empty random part.
    ## nlme then aborts inside model.matrix.reStruct() with "'data' must be of
    ## a vector type, was 'NULL'", which points at neither cause. The shared
    ## helper raises the same guarded error as the consolidated engines (B53).
    randeff <- as.formula(.countimp_randeff(
      grp.name  = nam[grp],
      slopes    = if (any(type == 2)) nam[which(type == 2)] else character(0),
      intercept = intercept,
      part      = "2l.pmm",
      style     = "nlme"))
    fixedeff <- as.formula(paste("Y", "~", fixedeff, sep = ""))
    dat <- data.frame(Y, X)
    ## verbose = FALSE: glmmPQL()'s default is TRUE, which cats "Iteration i"
    ## for each of its niter = 10 PQL steps. Inside countimp() that is once per
    ## imputation per iteration per incomplete variable -- m = 5, maxit = 10 and
    ## three variables produce ~300 lines that say nothing about the imputation
    ## and bury real warnings. The consolidated engines do not call glmmPQL, so
    ## this was the only source of that noise.
    fit <- MASS::glmmPQL(fixed = fixedeff, data = dat, random = randeff, 
                   family = .countimp_family("gaussian"), control = list(opt =
                                                         "optim"), na.action = na.omit,
                   verbose = FALSE)
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
    ## Clusters that have no observed y at all are unknown to the fit.
    ## glmmPQL returns NA for them; the NA then reaches .pmm.match(), where
    ## sort.int(d, partial = donors) aborts with "'partial' index outside
    ## bounds" -- a message that blames the donor pool. .countimp_predict_2l()
    ## draws the unobserved cluster effect from N(0, tau^2) instead; see
    ## predict2l.R for why u_j = 0 is not an acceptable substitute.
    obs.levels <- .countimp_obs_levels(dat, nam[grp])
    yhatmis <- .countimp_predict_2l(fitmis, newdatamis, type = "response",
                                    grp = nam[grp], obs_levels = obs.levels)
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
