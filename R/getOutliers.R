# note: functions were imported from package extremevalues 2.3.2 2016-01-05
getOutliers<-function (y, rho = c(1, 1), FLim = c(0.1, 0.9), distribution = "normal") 
{
  if (!is.vector(y)) 
    stop("First argument is not of type vector")
  if (sum(y < 0) > 0 & !(distribution == "normal")) 
    stop("First argument contains nonpositive values")
  if (sum(rho <= 0, na.rm = TRUE) > 0) 
    stop("Values of rho must be positive")
  if (FLim[2] <= FLim[1] | sum(FLim < 0 | FLim > 1) > 0) 
    stop("Invalid range in FLim: 0<=FLim[1]<FLim[2]<=1")
  if (!distribution %in% c("lognormal", "pareto", "exponential", 
                           "weibull", "normal")) 
    stop("Invalid distribution (lognormal, pareto, exponential, weibull, normal).")
  Y <- y
  y <- sort(y)
  N <- length(y)
  P <- (1:N)/(N + 1)
  Lambda <- P >= FLim[1] & P <= FLim[2]
  y <- y[Lambda]
  p <- P[Lambda]
  out <- switch(distribution, lognormal = getLognormalLimit(y, 
                                                            p, N, rho), pareto = getParetoLimit(y, p, N, rho), exponential = getExponentialLimit(y, 
                                                                                                                                                 p, N, rho), weibull = getWeibullLimit(y, p, N, rho), 
                normal = getNormalLimit(y, p, N, rho))
  out$method <- "Method I"
  out$distribution = distribution
  out$iRight = which(Y > out$limit[2])
  out$iLeft = which(Y < out$limit[1])
  out$nOut = c(Left = length(out$iLeft), Right = length(out$iRight))
  out$yMin <- y[1]
  out$yMax <- tail(y, 1)
  out$rho = c(Left = rho[1], Right = rho[2])
  out$Fmin = FLim[1]
  out$Fmax = FLim[2]
  return(out)
}

getNormalLimit<-function (y, p, N, rho) 
{
  param <- fitNormal(y, p)
  ell <- c(Left = -Inf, Right = Inf)
  if (!is.na(rho[1])) 
    ell[1] <- sqrt(2) * param$sigma * invErf(2 * rho[1]/N - 
                                               1) + param$mu
  if (!is.na(rho[2])) 
    ell[2] <- sqrt(2) * param$sigma * invErf(1 - 2 * rho[2]/N) + 
    param$mu
  return(list(mu = param$mu, sigma = param$sigma, nFit = length(y), 
              R2 = param$R2, limit = ell))
}

invErf <- function (x) 
{
  if (sum(x >= 1) > 0 | sum(x <= -1) > 0) 
    stop("Argument must be between -1 and 1")
  return(qnorm((1 + x)/2)/sqrt(2))
}

fitNormal<-function (y, p) 
{
  if (!is.vector(y)) 
    stop("First argument is not of type vector")
  if (!is.vector(p)) 
    stop("First argument is not of type vector")
  if (sum(p <= 0) > 0 | sum(p >= 1) > 0) 
    stop("Second argument contains values out of range (0,1)")
  if (length(y) != length(p)) 
    stop("First and second argument have different length")
  N <- length(y)
  Y <- as.matrix(y, nrow = N)
  p <- as.matrix(p, nrow = N)
  A <- matrix(0, nrow = N, ncol = 2)
  A[, 1] <- 1 + double(N)
  A[, 2] <- sqrt(2) * invErf(2 * p - 1)
  param <- solve(t(A) %*% A) %*% t(A) %*% Y
  r2 <- 1 - var(A %*% param - y)/var(y)
  return(list(mu = param[1], sigma = param[2], R2 = r2))
}