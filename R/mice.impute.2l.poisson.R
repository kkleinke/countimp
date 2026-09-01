#' Multiple Imputation of Incomplete Two-Level Count Data
#'
#' The functions impute multilevel count data based on a two-level Poisson or
#' negative binomial model, either using a Bayesian regression variant or a
#' cluster bootstrap variant. Since version 3.0.0 all twelve exported names
#' share one engine, \code{.countimp_2l_count()}; they differ only in the
#' count family, the intercept treatment and the draw.
#'
#' For counts with more zeros than a Poisson or NB law can produce, see
#' \code{\link{mice.impute.2l.zinb}} (zero inflation: two ways to observe a
#' zero) or \code{\link{mice.impute.2l.hnb}} (hurdle: one gate decides zero
#' versus positive, and the positive part cannot be zero).
#'
#' Model specification details:
#'    \itemize{
#'      \item -2 = grouping (class) variable. More than one may be coded: two
#'        entries give a three-level model, e.g. pupils in classes in schools.
#'        The FIRST one carries the random slopes and the \code{.noint} flag;
#'        further levels enter as plain random intercepts.
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
#' The bootstrap functions draw a bootstrap sample from \code{y[ry]} and
#' \code{x[ry,]}. Note that whole clusters are resampled rather than
#' individual rows, so that the two-level structure is preserved:
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
#' @param EV should automatic outlier handling of imputed values be enabled? Default is \code{FALSE} since version 2.0.8. Setting \code{EV = TRUE} identifies extreme imputations and replaces them by predictive mean matching draws (\code{mice.impute.midastouch()}). This is \strong{not recommended}: for count data the upper tail is part of the target distribution, not an artefact. In a simulation study (n = 500, 35 \% MAR, R = 1000) \code{EV = TRUE} biased the regression slope by -14.3 \% and reduced the coverage of the nominal 95 \% confidence interval to 79.2 \% (versus 94.1 \% with \code{EV = FALSE}). The option is retained only for reproducing earlier analyses.
#' @return Numeric vector of length \code{sum(!ry)} with imputations
#' @export
#' @import glmmTMB
#' @importFrom MASS rnegbin
#' @aliases mice.impute.2l.poisson.boot mice.impute.2l.poisson.noint mice.impute.2l.poisson.noint.boot
#' @aliases mice.impute.2l.nb mice.impute.2l.nb.boot mice.impute.2l.nb.noint mice.impute.2l.nb.noint.boot
#' @aliases mice.impute.2l.nb2 mice.impute.2l.nb2.boot mice.impute.2l.nb2.noint mice.impute.2l.nb2.noint.boot
#' @author Kristian Kleinke
#' @describeIn mice.impute.2l.poisson Bayesian Poisson regression variant; random intercept
mice.impute.2l.poisson <-
  function(y, ry, x, type, intercept = TRUE, wy = NULL, EV = FALSE) {
    .countimp_2l_count(y, ry, x, type, family = "poisson", draw = "bayes",
                       intercept = intercept, wy = wy, EV = EV)
  }

#' @describeIn mice.impute.2l.poisson Bayesian Poisson regression variant; fixed intercept
#' @export
mice.impute.2l.poisson.noint <-
  function(y, ry, x, type, intercept = FALSE, wy = NULL, EV = FALSE) {
    .countimp_2l_count(y, ry, x, type, family = "poisson", draw = "bayes",
                       intercept = intercept, wy = wy, EV = EV)
  }

#' @describeIn mice.impute.2l.poisson Bootstrap Poisson regression variant; random intercept
#' @export
mice.impute.2l.poisson.boot <-
  function(y, ry, x, type, intercept = TRUE, wy = NULL, EV = FALSE) {
    .countimp_2l_count(y, ry, x, type, family = "poisson", draw = "boot",
                       intercept = intercept, wy = wy, EV = EV)
  }

#' @describeIn mice.impute.2l.poisson Bootstrap Poisson regression variant; fixed intercept
#' @export
mice.impute.2l.poisson.noint.boot <-
  function(y, ry, x, type, intercept = FALSE, wy = NULL, EV = FALSE) {
    .countimp_2l_count(y, ry, x, type, family = "poisson", draw = "boot",
                       intercept = intercept, wy = wy, EV = EV)
  }

#' @describeIn mice.impute.2l.poisson Bayesian NB regression variant; random intercept
#' @export
mice.impute.2l.nb <-
  function(y, ry, x, type, intercept = TRUE, wy = NULL, EV = FALSE) {
    .countimp_2l_count(y, ry, x, type, family = "nbinom2", draw = "bayes",
                       intercept = intercept, wy = wy, EV = EV)
  }

#' @describeIn mice.impute.2l.poisson Bayesian NB regression variant; fixed intercept
#' @export
mice.impute.2l.nb.noint <-
  function(y, ry, x, type, intercept = FALSE, wy = NULL, EV = FALSE) {
    .countimp_2l_count(y, ry, x, type, family = "nbinom2", draw = "bayes",
                       intercept = intercept, wy = wy, EV = EV)
  }

#' @describeIn mice.impute.2l.poisson Bootstrap NB regression variant; random intercept
#' @export
mice.impute.2l.nb.boot <-
  function(y, ry, x, type, intercept = TRUE, wy = NULL, EV = FALSE) {
    .countimp_2l_count(y, ry, x, type, family = "nbinom2", draw = "boot",
                       intercept = intercept, wy = wy, EV = EV)
  }

#' @describeIn mice.impute.2l.poisson Bootstrap NB regression variant; fixed intercept
#' @export
mice.impute.2l.nb.noint.boot <-
  function(y, ry, x, type, intercept = FALSE, wy = NULL, EV = FALSE) {
    .countimp_2l_count(y, ry, x, type, family = "nbinom2", draw = "boot",
                       intercept = intercept, wy = wy, EV = EV)
  }

## The 2l.nb2* names are documented as identical to their 2l.nb* originals
## and kept for backward compatibility. Up to 3.0.0 they were full copies;
## as real assignments they cannot drift apart from the originals.

#' @describeIn mice.impute.2l.poisson identical to \code{mice.impute.2l.nb}; kept for backward compatibility
#' @export
mice.impute.2l.nb2 <- mice.impute.2l.nb

#' @describeIn mice.impute.2l.poisson identical to \code{mice.impute.2l.nb.noint}; kept for backward compatibility
#' @export
mice.impute.2l.nb2.noint <- mice.impute.2l.nb.noint

#' @describeIn mice.impute.2l.poisson identical to \code{mice.impute.2l.nb.boot}; kept for backward compatibility
#' @export
mice.impute.2l.nb2.boot <- mice.impute.2l.nb.boot

#' @describeIn mice.impute.2l.poisson identical to \code{mice.impute.2l.nb.noint.boot}; kept for backward compatibility
#' @export
mice.impute.2l.nb2.noint.boot <- mice.impute.2l.nb.noint.boot
