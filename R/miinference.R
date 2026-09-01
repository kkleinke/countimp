#' Multiple Imputation Inference
#' 
#' Multiple Imputation Inference using Rubin's (1987) rules for MI inference.
#' 
#' Combines \code{m} sets of parameter estimates and standard errors into an overall set of estimate using Rubin's (1987) rules for multiple imputation inference. This function is an adaptation of function \code{mi.inference} from package \pkg{norm}.
#'
#' @param est one of three things: (a) a list of \code{m} vectors of parameter
#'   estimates, in which case \code{std.err} must be given as well; (b) a list
#'   of the \code{m} fitted models, from which estimates, standard errors and
#'   \code{dfcom} are taken automatically; (c) an object with an
#'   \code{analyses} component holding those fits, which is what \code{with()}
#'   returns. Forms (b) and (c) are detected by structure, so they work
#'   without \pkg{mice} being installed.
#' @param std.err a list of \code{m} vectors with corresponding standard
#'   errors. Only needed for input form (a); ignored otherwise.
#' @param confidence confidence level; default is to compute 95\% confidence intervals
#' @param dfcom the residual degrees of freedom of the complete-data analysis.
#'   This single argument selects the degrees-of-freedom convention, because
#'   Rubin's (1987) formula is the limiting case of Barnard and Rubin (1999)
#'   as \code{dfcom} goes to infinity:
#'   \describe{
#'     \item{\code{dfcom = Inf}}{Rubin (1987). The default, so that results
#'       from earlier versions of this function reproduce exactly. The formula
#'       has no upper bound and can therefore return degrees of freedom above
#'       the complete-data degrees of freedom.}
#'     \item{\code{dfcom} finite}{Barnard and Rubin (1999), which bounds the
#'       degrees of freedom by the complete-data degrees of freedom. This
#'       reproduces the \code{df} reported by \code{mice::pool}.}
#'   }
#'   When \code{est} is given as fitted models and \code{dfcom} is not
#'   supplied, it is read off the first fit via \code{df.residual}, so the
#'   convenient call also gets the bounded degrees of freedom. The convention
#'   actually used is reported in the \code{df.method} element of the result.
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
#' \item{m}{number of imputations that were combined}
#' \item{dfcom}{the complete-data degrees of freedom that were used}
#' \item{df.method}{\code{"rubin1987"} or \code{"barnard.rubin"}, depending on
#'   whether \code{dfcom} was infinite or finite}
#' }
#' @references 
#' \itemize{
#' \item Rubin, D. B. (1987). \emph{Multiple imputation for nonresponse in surveys}. New York: Wiley.
#' \item Barnard, J., & Rubin, D. B. (1999). Small-sample degrees of freedom with multiple imputation. \emph{Biometrika}, 86(4), 948--955.
#' \item Schafer, J. L. (1997). \emph{Analysis of incomplete multivariate data}. London: Chapman & Hall.
#' }
#' @aliases miinference
#' @examples 
#' ## simulate overdispersed count data
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
#' imp <- countimp( MNB, method=c("nb","","",""), printFlag = FALSE )
#' res <- with( imp, MASS::glm.nb(y~x1+x2+x3) )
#'
#' ## pooling the regression coefficients needs no further arguments: the
#' ## estimates and standard errors are read from the fits
#' miinference( res )

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
miinference <- function (est, std.err = NULL, confidence = 0.95, dfcom = Inf)
{
  ## Degrees of freedom: one formula, two conventions -------------------------
  ## Rubin (1987, eq. 3.1.6) and Barnard & Rubin (1999) are not competing
  ## methods -- Rubin's is the limiting case of Barnard-Rubin as the
  ## complete-data df go to infinity. So `dfcom` is the only switch needed:
  ##   dfcom = Inf     -> Rubin (1987), the historical default of this
  ##                      function; kept as the default so that numbers
  ##                      published with earlier versions still reproduce.
  ##   dfcom = finite  -> Barnard-Rubin, which bounds the df by the
  ##                      complete-data df. Rubin's formula has no upper
  ##                      bound and can return df above the complete-data df
  ##                      (in a check with n = 150 that happened in 7% of
  ##                      replications), which cannot be right.
  ## When `est` is a list of fitted models, dfcom is read off the fits, so the
  ## convenient call is also the better-behaved one.

  ## Accept fitted models -----------------------------------------------------
  ## Three input forms: a list of m fitted models, an object with an
  ## $analyses list (what with() returns), or the classic pair of numeric
  ## lists. Detected by structure, not by class, so this does not depend on
  ## any particular imputation package being installed.
  auto_dfcom <- NULL
  if (is.list(est) && !is.null(est$analyses) && is.list(est$analyses))
    est <- est$analyses
  is_numlist <- is.list(est) && length(est) > 0L &&
    all(vapply(est, function(z) is.numeric(z) && is.null(dim(z)), NA))
  if (is.list(est) && !is_numlist) {
    fits <- est
    ## A list of NULLs reaches this branch (NULL is not numeric), and coef()
    ## accepts NULL without complaint -- the failure then surfaced from
    ## vcov(NULL) as "no applicable method for 'vcov'", which names neither
    ## the argument nor the cause. The usual way to get here is an extraction
    ## that silently returned nothing, e.g. lapply(res$analyses, `[[`, "theta")
    ## on fits that do not carry `theta`, so say that.
    empty <- vapply(fits, function(f) is.null(f) || length(f) == 0L, NA)
    if (any(empty))
      stop("`est` contains ", sum(empty), " empty element(s) (NULL or ",
           "length 0). If these came from extracting a component of fitted ",
           "models, check that the component exists in every fit.",
           call. = FALSE)
    ok <- vapply(fits, function(f)
      !inherits(try(stats::coef(f), silent = TRUE), "try-error") &&
      !inherits(try(stats::vcov(f), silent = TRUE), "try-error"), NA)
    if (!all(ok))
      stop("`est` must be a list of numeric vectors, a list of fitted ",
           "models, or an object with an `analyses` component.", call. = FALSE)
    est     <- lapply(fits, stats::coef)
    std.err <- lapply(fits, function(f) sqrt(diag(as.matrix(stats::vcov(f)))))
    auto_dfcom <- suppressWarnings(
      tryCatch(as.numeric(stats::df.residual(fits[[1L]])),
               error = function(e) NULL))
    if (missing(dfcom) && !is.null(auto_dfcom) &&
        length(auto_dfcom) == 1L && is.finite(auto_dfcom) && auto_dfcom > 0)
      dfcom <- auto_dfcom
  }

  ## Input checks -------------------------------------------------------------
  ## The old code went straight to `for (i in 2:length(est))`, which counts
  ## DOWN when length(est) == 1 and then failed with "subscript out of bounds".
  ## m = 1 is not an error the user can debug from that message, so name it.
  if (!is.list(est) || !is.list(std.err))
    stop("`est` and `std.err` must each be a list of m vectors.", call. = FALSE)
  m <- length(est)
  if (m != length(std.err))
    stop(sprintf("`est` has %d elements but `std.err` has %d.", m, length(std.err)),
         call. = FALSE)
  if (m < 2)
    stop("Multiple-imputation inference needs m >= 2 sets of estimates; got m = ",
         m, ". With a single imputation there is no between-imputation\n",
         "  variance to pool, so Rubin's rules do not apply.", call. = FALSE)
  k <- length(est[[1]])
  ## Zero-length estimates passed every check above and produced a table with
  ## no rows and no message -- the user asked for pooled inference and got an
  ## empty object back. Name it instead.
  if (k == 0L)
    stop("`est` has no parameters to pool (elements have length 0).",
         call. = FALSE)
  if (any(vapply(est, length, 1L) != k) || any(vapply(std.err, length, 1L) != k))
    stop("All elements of `est` and `std.err` must have the same length.",
         call. = FALSE)
  if (!is.numeric(dfcom) || length(dfcom) != 1L || is.na(dfcom) || dfcom <= 0)
    stop("`dfcom` must be a single positive number, or Inf for Rubin (1987) ",
         "degrees of freedom.", call. = FALSE)

  qstar <- do.call(cbind, est)
  u     <- do.call(cbind, std.err)^2
  qbar  <- rowMeans(qstar)
  ubar  <- rowMeans(u)                     # within-imputation variance
  bm    <- apply(qstar, 1, stats::var)     # between-imputation variance
  tm    <- ubar + (1 + 1/m) * bm           # total variance
  rem   <- (1 + 1/m) * bm/ubar             # relative increase in variance
  lambda <- (1 + 1/m) * bm/tm              # fraction of missing information

  nu_rubin <- (m - 1)/lambda^2             # Rubin (1987); == (m-1)*(1+1/rem)^2
  if (is.infinite(dfcom)) {
    nu <- nu_rubin
  } else {
    tmp <- (1 - lambda) * (1 + dfcom) * dfcom
    nu  <- (m - 1) * tmp/((dfcom + 3) * (m - 1) + lambda^2 * tmp)
  }

  alpha <- 1 - (1 - confidence)/2
  low   <- qbar - stats::qt(alpha, nu) * sqrt(tm)
  up    <- qbar + stats::qt(alpha, nu) * sqrt(tm)
  pval  <- 2 * (1 - stats::pt(abs(qbar/sqrt(tm)), nu))
  fminf <- (rem + 2/(nu + 3))/(rem + 1)

  list(est = qbar, std.err = sqrt(tm), t.value = qbar/sqrt(tm), df = nu,
       p.value = pval, lower = low, upper = up, r = rem, fminf = fminf,
       m = m, dfcom = dfcom,
       df.method = if (is.infinite(dfcom)) "rubin1987" else "barnard.rubin")
}
