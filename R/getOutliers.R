# note: functions were imported from package extremevalues 2.3.2 2016-01-05

## countimp 2.0.8: outlier handling is deprecated and no longer the default.
## Rationale: replacing large imputations by predictive mean matching draws
## truncates the right tail of the count distribution, which for overdispersed
## counts is part of the target distribution rather than an artefact. Measured
## consequence (n = 500, 35% MAR, R = 1000 replications, NB model):
##   EV = TRUE : rel. bias -14.3%, coverage of the nominal 95% CI = 79.2%
##   EV = FALSE: rel. bias  -0.3%, coverage                       = 94.1%
.countimp_state <- new.env(parent = emptyenv())
.countimp_state$ev_warned <- FALSE
.countimp_state$beta_route_ok <- FALSE      # see .countimp_check_beta_route() in predict2l.R
.countimp_state$theta_fixed_noted <- FALSE  # see .countimp_note_theta_fixed() in rtrunc.R
.countimp_state$underdisp_warned <- FALSE   # see .countimp_warn_underdisp() in rtrunc.R
.countimp_state$joint_noted <- FALSE       # see .countimp_draw_2l()
.countimp_state$joint_zi_noted <- FALSE    # see .countimp_draw_2l_zi()
.countimp_state$mixed_scale_warned <- FALSE # see .countimp_1l_censored()

.countimp_warn_ev <- function() {
  if (!isTRUE(.countimp_state$ev_warned)) {
    warning("EV = TRUE (automatic outlier handling) is deprecated and statistically ",
            "not recommended. It truncates the upper tail of the count distribution, ",
            "which biases point estimates downwards (approx. -14% in simulation) and ",
            "reduces coverage of the nominal 95% CI to about 79%. Use EV = FALSE ",
            "(the default since version 2.0.8). This warning is shown once per session.",
            call. = FALSE)
    .countimp_state$ev_warned <- TRUE
  }
  invisible(NULL)
}

getOutliers<-function (y, rho = c(1, 1), FLim = c(0.1, 0.9), distribution = "normal") 
{
  .countimp_warn_ev()
  if (!is.vector(y)) 
    stop("First argument is not of type vector")
  if (sum(y < 0) > 0 & !(distribution == "normal")) 
    stop("First argument contains nonpositive values")
  if (sum(rho <= 0, na.rm = TRUE) > 0) 
    stop("Values of rho must be positive")
  if (FLim[2] <= FLim[1] | sum(FLim < 0 | FLim > 1) > 0) 
    stop("Invalid range in FLim: 0<=FLim[1]<FLim[2]<=1")
  ## Only the normal branch is carried inside countimp. The lognormal, pareto,
  ## exponential and weibull branches of extremevalues::getOutliers call
  ## fitting helpers (getLognormalLimit and friends) that were NOT copied
  ## along with this function -- calling them would fail with "could not find
  ## function". Rejecting them here turns a confusing runtime error into a
  ## statement of what this copy supports, and points at the original.
  if (!identical(distribution, "normal"))
    stop("countimp's internal getOutliers() supports distribution = \"normal\" ",
         "only.\n",
         "  Requested: \"", distribution, "\".\n",
         "  For the other distributions use extremevalues::getOutliers()\n",
         "  directly (install.packages(\"extremevalues\")).", call. = FALSE)
  Y <- y
  y <- sort(y)
  N <- length(y)
  P <- (1:N)/(N + 1)
  Lambda <- P >= FLim[1] & P <= FLim[2]
  y <- y[Lambda]
  p <- P[Lambda]
  out <- getNormalLimit(y, p, N, rho)
  out$method <- "Method I"
  out$distribution = distribution
  out$iRight = which(Y > out$limit[2])
  out$iLeft = which(Y < out$limit[1])
  out$nOut = c(Left = length(out$iLeft), Right = length(out$iRight))
  out$yMin <- y[1]
  out$yMax <- utils::tail(y, 1)
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
