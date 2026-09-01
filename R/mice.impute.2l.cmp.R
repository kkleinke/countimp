#' Multiple Imputation of Underdispersed Two-Level Count Data
#'
#' Imputes clustered counts whose conditional variance is \emph{smaller} than
#' their mean, under a two-level Conway-Maxwell-Poisson model with a random
#' intercept (and random slopes where the type codes ask for them). This is the
#' two-level counterpart of \code{\link{mice.impute.cmp}}: the Poisson model
#' fixes the variance-to-mean ratio at one and the negative binomial can only
#' raise it, so both overstate the spread of underdispersed counts.
#'
#' @details
#' Counts bounded by design are the usual source: items correct out of a fixed
#' number, symptoms endorsed from a checklist, days attended in a school week --
#' all of them frequently clustered within schools, wards or classes.
#'
#' Model specification follows the other two-level methods:
#'    \itemize{
#'      \item -2 = grouping (class) variable. More than one may be coded: two
#'        entries give a three-level model, e.g. pupils in classes in schools.
#'        The FIRST one carries the random slopes and the \code{.noint} flag;
#'        further levels enter as plain random intercepts.
#'      \item 0 = variable not included in the imputation model
#'      \item 1 = variable included as a fixed effect
#'      \item 2 = variable included as a fixed \emph{and} random effect
#'    }
#'
#' The fit uses \pkg{glmmTMB}, as every two-level method in this package does --
#' random effects need the Laplace approximation. Two parts stay countimp's own,
#' and both are measured rather than assumed:
#' \itemize{
#'   \item \strong{the draw.} \pkg{glmmTMB}'s \code{simulate()} ignores
#'     \code{newdata} and returns values for the fitting cases, which is not
#'     what imputation needs; imputations come from countimp's own exact
#'     sampler.
#'   \item \strong{the shape.} \pkg{glmmTMB} reports \eqn{1/\nu} through
#'     \code{sigma()}. countimp scores both readings against the data with its
#'     own likelihood and uses whichever fits, so a future change of that
#'     convention cannot silently invert the dispersion of the imputations. A
#'     disagreement is recorded in \code{\link{countimp_diagnostics}}.
#' }
#'
#' \strong{Cost.} The COM-Poisson normalising constant has no closed form, so
#' both the fit and the draw sum over a grid. Expect these methods to be several
#' times slower than \code{\link{mice.impute.2l.poisson}}; use
#' \code{\link{countimp_fit_diag}} on the observed data first to see whether the
#' conditional variance-to-mean ratio warrants the family at all.
#'
#' @param y Numeric vector with incomplete data in long format (the groups are
#'   stacked upon each other).
#' @param ry Response pattern of \code{y} (\code{TRUE} = observed,
#'   \code{FALSE} = missing).
#' @param x Matrix with \code{length(y)} rows containing complete covariates,
#'   also in long format.
#' @param type Vector of length \code{ncol(x)} identifying fixed, random and
#'   class variables; extracted automatically from the
#'   \code{predictorMatrix}. See section Details.
#' @param intercept \code{TRUE}: the intercept enters as a random effect;
#'   \code{FALSE}: the intercept is treated as fixed.
#' @param wy Logical vector of length \code{length(y)}; \code{TRUE} marks the
#'   positions to impute. Defaults to \code{!ry}.
#' @param EV Should automatic outlier handling of imputed values be enabled?
#'   Default \code{FALSE}, and \strong{not recommended} -- see
#'   \code{\link{mice.impute.2l.poisson}} for the measured cost.
#'
#' @return Numeric vector of length \code{sum(wy)} with imputations.
#'
#' @references
#' Shmueli, G., Minka, T. P., Kadane, J. B., Borle, S., & Boatwright, P. (2005).
#' A useful distribution for fitting discrete data: revival of the
#' Conway-Maxwell-Poisson distribution. \emph{Journal of the Royal Statistical
#' Society C}, 54(1), 127-142.
#'
#' Huang, A. (2017). Mean-parametrized Conway-Maxwell-Poisson regression models
#' for dispersed counts. \emph{Statistical Modelling}, 17(6), 359-370.
#'
#' @seealso \code{\link{mice.impute.cmp}} for the single-level case,
#'   \code{\link{mice.impute.2l.poisson}} for equidispersed and overdispersed
#'   clustered counts, \code{\link{countimp_fit_diag}} for choosing between
#'   them.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' ngrp <- 20; nj <- 15; n <- ngrp * nj
#' id <- rep(seq_len(ngrp), each = nj)
#' x1 <- rnorm(n)
#' u  <- rnorm(ngrp, 0, 0.3)
#' ## binomial counts are underdispersed: variance below the mean
#' y  <- rbinom(n, 10, plogis(-0.4 + 0.5 * x1 + u[id]))
#' d  <- data.frame(y = y, x1 = x1, id = id)
#' d$y[sample(n, 60)] <- NA
#' pred <- matrix(c(0, 0, 0,  1, 0, 0,  -2, 0, 0), nrow = 3, byrow = TRUE,
#'                dimnames = list(names(d), names(d)))
#' pred["y", ] <- c(0, 1, -2)
#' imp <- countimp(d, method = c(y = "2l.cmp", x1 = "", id = ""),
#'                 predictorMatrix = pred, m = 2, maxit = 1, printFlag = FALSE)
#' }
#'
#' @author Kristian Kleinke
#' @aliases mice.impute.2l.cmp.noint mice.impute.2l.cmp.boot
#' @aliases mice.impute.2l.cmp.noint.boot
#' @export
#' @describeIn mice.impute.2l.cmp Bayesian variant; random intercept
mice.impute.2l.cmp <-
  function(y, ry, x, type, intercept = TRUE, wy = NULL, EV = FALSE) {
    .countimp_2l_cmp(y, ry, x, type, draw = "bayes", intercept = intercept,
                     wy = wy, EV = EV)
  }

#' @describeIn mice.impute.2l.cmp Bayesian variant; fixed intercept
#' @export
mice.impute.2l.cmp.noint <-
  function(y, ry, x, type, intercept = FALSE, wy = NULL, EV = FALSE) {
    .countimp_2l_cmp(y, ry, x, type, draw = "bayes", intercept = intercept,
                     wy = wy, EV = EV)
  }

#' @describeIn mice.impute.2l.cmp Cluster bootstrap variant; random intercept
#' @export
mice.impute.2l.cmp.boot <-
  function(y, ry, x, type, intercept = TRUE, wy = NULL, EV = FALSE) {
    .countimp_2l_cmp(y, ry, x, type, draw = "boot", intercept = intercept,
                     wy = wy, EV = EV)
  }

#' @describeIn mice.impute.2l.cmp Cluster bootstrap variant; fixed intercept
#' @export
mice.impute.2l.cmp.noint.boot <-
  function(y, ry, x, type, intercept = FALSE, wy = NULL, EV = FALSE) {
    .countimp_2l_cmp(y, ry, x, type, draw = "boot", intercept = intercept,
                     wy = wy, EV = EV)
  }
