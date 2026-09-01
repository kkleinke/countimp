#' Multiple Imputation of Zero-Inflated Two-Level Count Data
#'
#' The functions impute zero-inflated multilevel count data based on a two-level Poisson or negative binomial zero-inflation model, either using a Bayesian regression or a bootstrap regression approach (appendix: ``\code{.boot}''). The \code{.noint} variants treat the intercept only as a fixed, but \emph{not} as a random effect. It may be specified, if the intercept is excluded from the random part of the zero model (``\code{.noint.zero}''), the count model (``\code{.noint.count}''), or both models (``\code{.noint.both}''). Zero-inflation models are mixture models and consist of two model components: the zero model (a binomial generalized linear mixed effects model), determining, if the observational unit belongs to zero-inflation process (certain zeros) or to the count process, and the count model (a two-level Poisson or NB model), determining, what count (zero or non-zero) the observational unit has.
#'
#' Model specification details:
#'    \itemize{
#'      \item -2 = grouping (class) variable. More than one may be coded: two
#'        entries give a three-level model, e.g. pupils in classes in schools.
#'        The FIRST one carries the random slopes and the \code{.noint} flag;
#'        further levels enter as plain random intercepts.
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
#'  \item Fit the zero-inflation model -- \emph{one} \code{glmmTMB} model with both a conditional (count) part and a zero-inflation part -- to the observed data.
#'  \item Draw b* = (beta*, u*, betazi*, log theta*) from the joint Laplace posterior of \emph{all} model parameters. Because the two model parts are estimated jointly and are correlated, one draw must move both linear predictors; drawing them separately would discard that correlation.
#'  \item Compute pi*, the individual probability of belonging to the zero-inflation process, and mu*, the individual count rate, from b*.
#'  \item For each incomplete case draw an indicator from Bernoulli(pi*). Cases in the zero-inflation process are imputed as zero.
#'  \item Impute the remaining cases by a draw from the (untruncated) Poisson or NB count distribution with rate mu*. Note that this draw may itself return a zero: under zero inflation a count-process unit may legitimately have a zero count. This is what distinguishes a zero-inflation model from a hurdle model, for which \code{\link{mice.impute.2l.hnb}} and \code{\link{mice.impute.2l.hp}} should be used instead.
#'  }
#' The bootstrap functions resample \emph{clusters} (not individual cases) from \code{y[ry]} and \code{x[ry,]}, then fit the same zero-inflation model to the bootstrap sample and carry out steps 3 to 5 at the bootstrap estimates. The resampling itself supplies the parameter uncertainty, so no posterior draw is added.
#' @param y Numeric vector with incomplete data in long format (i.e. the groups are stacked upon each other)
#' @param ry Response pattern of \code{y} (\code{TRUE}=observed, \code{FALSE}=missing)
#' @param x matrix with \code{length(y)} rows containing complete covariates; also in long format
#' @param type vector of length \code{ncol(x)} identifying fixed, random, and class variables; \code{type} is automatically extracted from the \code{predictorMatrix}; see \pkg{mice}'s user's manual for details about how to specify the imputation model; see also section ``details''.
#' @param intercept.c \code{TRUE}: model will include intercept as a random effect in the count model; \code{FALSE}: count model intercept will be treated as fixed.
#' @param intercept.z \code{TRUE}: model will include intercept as a random effect in the zero model; \code{FALSE}: zero model intercept will be treated as fixed.
#' @param wy Logical vector of length \code{length(y)}. A \code{TRUE} value indicates locations in \code{y} for which imputations are created. Default is \code{!ry}
#' @param EV should automatic outlier handling of imputed values be enabled? Default is \code{FALSE} since version 2.0.8. Setting \code{EV = TRUE} identifies extreme imputations and replaces them by predictive mean matching draws (\code{mice.impute.midastouch()}). This is \strong{not recommended}: for count data the upper tail is part of the target distribution, not an artefact. In a simulation study (n = 500, 35 \% MAR, R = 1000) \code{EV = TRUE} biased the regression slope by -14.3 \% and reduced the coverage of the nominal 95 \% confidence interval to 79.2 \% (versus 94.1 \% with \code{EV = FALSE}). The option is retained only for reproducing earlier analyses.
#' @return Numeric vector of length \code{sum(!ry)} with imputations
#' @export
#' @import glmmTMB stats
#' @importFrom stats as.formula na.pass predict rbinom rnorm vcov
#' @aliases mice.impute.2l.zinb mice.impute.2l.zinb.boot mice.impute.2l.zinb.noint.both mice.impute.2l.zinb.noint.both.boot mice.impute.2l.zinb.noint.zero mice.impute.2l.zinb.noint.zero.boot mice.impute.2l.zinb.noint.count mice.impute.2l.zinb.noint.count.boot
#' @aliases mice.impute.2l.zip mice.impute.2l.zip.boot mice.impute.2l.zip.noint.both mice.impute.2l.zip.noint.both.boot mice.impute.2l.zip.noint.zero mice.impute.2l.zip.noint.zero.boot mice.impute.2l.zip.noint.count mice.impute.2l.zip.noint.count.boot
#' @section Corrected in 3.0.0:
#'   Up to 2.6.0 these methods fitted a \emph{hurdle} model (a binomial GLMM
#'   on 1\{y = 0\} plus an untruncated count GLMM on the positive subset) but
#'   drew from it as if the count part were unrestricted. Two errors compounded:
#'   the count model was fitted to positive-only data without accounting for the
#'   truncation, and a case assigned to the count process could still be imputed
#'   as zero. Both push in the same direction -- too many zeros, conditional
#'   mean too small. In simulation (R = 200, m = 5, n = 625, 38-40 \% MAR) the
#'   coverage of the nominal 95 \% interval for the marginal covariate effect
#'   was 88.5 \% under a true zero-inflation process and 87.0 \% under a true
#'   hurdle process. These methods now fit a genuine zero-inflation model, which
#'   restores coverage to 91.5 \% under a zero-inflation process; for hurdle
#'   data use \code{\link{mice.impute.2l.hnb}} / \code{\link{mice.impute.2l.hp}}
#'   (92.0 \%). Imputations from earlier versions are not comparable.
#'   The replication script is \code{skripte/b30_zi_comparison.R}.
#' @section Corrected in 2.6.0:
#'   Two defects affected every two-level hurdle and ZI method up to 2.5.0.
#'   (1) \code{intercept.z} was declared but never used: both model parts were
#'   governed by \code{intercept.c}, so \code{.noint.zero} had no effect and
#'   \code{.noint.count} removed the random intercept from both parts.
#'   (2) In the count model the predictor set built from \code{type} was
#'   overwritten by all variables except the class variable, so codes 3-6 did
#'   not restrict a predictor to one model part. Both are fixed; imputations
#'   from earlier versions are not comparable where codes 3-6 were used.
#' @author Kristian Kleinke
#' @describeIn mice.impute.2l.zinb Bayesian regression variant; random intercepts
mice.impute.2l.zinb <-
  function(y, ry, x, type, intercept.c = TRUE, intercept.z = TRUE,
           wy = NULL, EV = FALSE)
    .countimp_2l_zi(y, ry, x, type, family = "nbinom2", draw = "bayes",
                    intercept.c = intercept.c, intercept.z = intercept.z,
                    wy = wy, EV = EV)


#' @export
#' @describeIn mice.impute.2l.zinb bootstrap regression variant; random intercepts
mice.impute.2l.zinb.boot <-
  function(y, ry, x, type, intercept.c = TRUE, intercept.z = TRUE,
           wy = NULL, EV = FALSE)
    .countimp_2l_zi(y, ry, x, type, family = "nbinom2", draw = "boot",
                    intercept.c = intercept.c, intercept.z = intercept.z,
                    wy = wy, EV = EV)


#' @export
#' @describeIn mice.impute.2l.zinb Bayesian regression variant; fixed intercepts
mice.impute.2l.zinb.noint.both <-
  function(y, ry, x, type, intercept.c = FALSE, intercept.z = FALSE,
           wy = NULL, EV = FALSE)
    .countimp_2l_zi(y, ry, x, type, family = "nbinom2", draw = "bayes",
                    intercept.c = intercept.c, intercept.z = intercept.z,
                    wy = wy, EV = EV)


#' @export
#' @describeIn mice.impute.2l.zinb bootstrap regression variant; fixed intercepts
mice.impute.2l.zinb.noint.both.boot <-
  function(y, ry, x, type, intercept.c = FALSE, intercept.z = FALSE,
           wy = NULL, EV = FALSE)
    .countimp_2l_zi(y, ry, x, type, family = "nbinom2", draw = "boot",
                    intercept.c = intercept.c, intercept.z = intercept.z,
                    wy = wy, EV = EV)


#' @export
#' @describeIn mice.impute.2l.zinb Bayesian regression variant; fixed intercept in count model
mice.impute.2l.zinb.noint.count <-
  function(y, ry, x, type, intercept.c = FALSE, intercept.z = TRUE,
           wy = NULL, EV = FALSE)
    .countimp_2l_zi(y, ry, x, type, family = "nbinom2", draw = "bayes",
                    intercept.c = intercept.c, intercept.z = intercept.z,
                    wy = wy, EV = EV)


#' @export
#' @describeIn mice.impute.2l.zinb bootstrap regression variant; fixed intercept in count model
mice.impute.2l.zinb.noint.count.boot <-
  function(y, ry, x, type, intercept.c = FALSE, intercept.z = TRUE,
           wy = NULL, EV = FALSE)
    .countimp_2l_zi(y, ry, x, type, family = "nbinom2", draw = "boot",
                    intercept.c = intercept.c, intercept.z = intercept.z,
                    wy = wy, EV = EV)


#' @export
#' @describeIn mice.impute.2l.zinb Bayesian regression variant; fixed intercept in zero model
mice.impute.2l.zinb.noint.zero <-
  function(y, ry, x, type, intercept.c = TRUE, intercept.z = FALSE,
           wy = NULL, EV = FALSE)
    .countimp_2l_zi(y, ry, x, type, family = "nbinom2", draw = "bayes",
                    intercept.c = intercept.c, intercept.z = intercept.z,
                    wy = wy, EV = EV)


#' @export
#' @describeIn mice.impute.2l.zinb bootstrap regression variant; fixed intercept in zero model
mice.impute.2l.zinb.noint.zero.boot <-
  function(y, ry, x, type, intercept.c = TRUE, intercept.z = FALSE,
           wy = NULL, EV = FALSE)
    .countimp_2l_zi(y, ry, x, type, family = "nbinom2", draw = "boot",
                    intercept.c = intercept.c, intercept.z = intercept.z,
                    wy = wy, EV = EV)


#' @export
#' @describeIn mice.impute.2l.zinb Bayesian regression variant; random intercepts
mice.impute.2l.zip <-
  function(y, ry, x, type, intercept.c = TRUE, intercept.z = TRUE,
           wy = NULL, EV = FALSE)
    .countimp_2l_zi(y, ry, x, type, family = "poisson", draw = "bayes",
                    intercept.c = intercept.c, intercept.z = intercept.z,
                    wy = wy, EV = EV)


#' @export
#' @describeIn mice.impute.2l.zinb bootstrap regression variant; random intercepts
mice.impute.2l.zip.boot <-
  function(y, ry, x, type, intercept.c = TRUE, intercept.z = TRUE,
           wy = NULL, EV = FALSE)
    .countimp_2l_zi(y, ry, x, type, family = "poisson", draw = "boot",
                    intercept.c = intercept.c, intercept.z = intercept.z,
                    wy = wy, EV = EV)


#' @export
#' @describeIn mice.impute.2l.zinb Bayesian regression variant; fixed intercepts
mice.impute.2l.zip.noint.both <-
  function(y, ry, x, type, intercept.c = FALSE, intercept.z = FALSE,
           wy = NULL, EV = FALSE)
    .countimp_2l_zi(y, ry, x, type, family = "poisson", draw = "bayes",
                    intercept.c = intercept.c, intercept.z = intercept.z,
                    wy = wy, EV = EV)


#' @export
#' @describeIn mice.impute.2l.zinb bootstrap regression variant; fixed intercepts
mice.impute.2l.zip.noint.both.boot <-
  function(y, ry, x, type, intercept.c = FALSE, intercept.z = FALSE,
           wy = NULL, EV = FALSE)
    .countimp_2l_zi(y, ry, x, type, family = "poisson", draw = "boot",
                    intercept.c = intercept.c, intercept.z = intercept.z,
                    wy = wy, EV = EV)


#' @export
#' @describeIn mice.impute.2l.zinb Bayesian regression variant; fixed intercept in count model
mice.impute.2l.zip.noint.count <-
  function(y, ry, x, type, intercept.c = FALSE, intercept.z = TRUE,
           wy = NULL, EV = FALSE)
    .countimp_2l_zi(y, ry, x, type, family = "poisson", draw = "bayes",
                    intercept.c = intercept.c, intercept.z = intercept.z,
                    wy = wy, EV = EV)


#' @export
#' @describeIn mice.impute.2l.zinb bootstrap regression variant; fixed intercept in count model
mice.impute.2l.zip.noint.count.boot <-
  function(y, ry, x, type, intercept.c = FALSE, intercept.z = TRUE,
           wy = NULL, EV = FALSE)
    .countimp_2l_zi(y, ry, x, type, family = "poisson", draw = "boot",
                    intercept.c = intercept.c, intercept.z = intercept.z,
                    wy = wy, EV = EV)


#' @export
#' @describeIn mice.impute.2l.zinb Bayesian regression variant; fixed intercept in zero model
mice.impute.2l.zip.noint.zero <-
  function(y, ry, x, type, intercept.c = TRUE, intercept.z = FALSE,
           wy = NULL, EV = FALSE)
    .countimp_2l_zi(y, ry, x, type, family = "poisson", draw = "bayes",
                    intercept.c = intercept.c, intercept.z = intercept.z,
                    wy = wy, EV = EV)


#' @export
#' @describeIn mice.impute.2l.zinb bootstrap regression variant; fixed intercept in zero model
mice.impute.2l.zip.noint.zero.boot <-
  function(y, ry, x, type, intercept.c = TRUE, intercept.z = FALSE,
           wy = NULL, EV = FALSE)
    .countimp_2l_zi(y, ry, x, type, family = "poisson", draw = "boot",
                    intercept.c = intercept.c, intercept.z = intercept.z,
                    wy = wy, EV = EV)
