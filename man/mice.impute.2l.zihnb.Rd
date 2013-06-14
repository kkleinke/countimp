\name{mice.impute.2l.zihnb}
\alias{mice.impute.2l.zihnb}
\alias{mice.impute.2l.zihnb.noint.zero}
\alias{mice.impute.2l.zihnb.noint.count}
\alias{mice.impute.2l.zihnb.noint.both}
\alias{mice.impute.2l.zihnb.boot}
\alias{mice.impute.2l.zihnb.noint.zero.boot}
\alias{mice.impute.2l.zihnb.noint.count.boot}
\alias{mice.impute.2l.zihnb.noint.both.boot}
\title{
Multiple Imputation of Zero-Inflated Two-Level Count Data
}
\description{
The functions impute zero-inflated (and overdispersed) multilevel count data based on a two-level Negative Binomial hurdle model, either using a Bayesian regression (``\code{.zihnb}'') or a bootstrap approach (``\code{.boot}''). The \code{.noint} variants treat the intercept only as a fixed, but \emph{not} as a random effect. It may be specified, if the intercept is excluded from the random effects only in the zero model (``\code{.noint.zero}''), the count model (``\code{.noint.count}''), or both models (``\code{.noint.both}''). Hurdle models are mixture models and consist of two models: the zero model (a binomial GLMM), determining, if the observational unit has a zero or non-zero value, and the count model (a zero-truncated two-level NB model) determining, what kind of non-zero value the observational unit has.
}
\usage{
mice.impute.2l.zihnb(y,ry,x,type,intercept.c=TRUE,intercept.z=TRUE)
mice.impute.2l.zihnb.noint.zero(y,ry,x,type,intercept.c=TRUE,intercept.z=FALSE)
mice.impute.2l.zihnb.noint.count(y,ry,x,type,intercept.c=FALSE,intercept.z=TRUE)
mice.impute.2l.zihnb.noint.both(y,ry,x,type,intercept.c=FALSE,intercept.z=FALSE)
mice.impute.2l.zihnb.boot(y,ry,x,type,intercept.c=TRUE,intercept.z=TRUE)
mice.impute.2l.zihnb.noint.zero.boot(y,ry,x,type,intercept.c=TRUE,intercept.z=FALSE)
mice.impute.2l.zihnb.noint.count.boot(y,ry,x,type,intercept.c=FALSE,intercept.z=TRUE)
mice.impute.2l.zihnb.noint.both.boot(y,ry,x,type,intercept.c=FALSE,intercept.z=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{y}{
Numeric vector with incomplete data in long format (i.e. the groups are stacked upon each other)
}
  \item{ry}{
Response pattern of \code{y} (\code{TRUE}=observed, \code{FALSE}=missing)
}
  \item{x}{
matrix with \code{length(y)} rows containing complete covariates; also in long format
}
  \item{type}{
vector of length \code{ncol(x)} determining the imputation model; \code{type} is automatically extracted from the \code{predictorMatrix} argument of \code{mice()}; see \pkg{mice}'s user's manual for details about how to specify the imputation model; see also the ``Details'' section. 
}
  \item{intercept.c}{\code{TRUE}: model will include intercept as a random effect in the count model; \code{FALSE}: model will not use intercept as a random effect.
}
\item{intercept.z}{\code{TRUE}: model will include intercept as a random effect in the zero model; \code{FALSE}: model will not use intercept as a random effect.
}
}
\details{
Model specification details: 
\itemize{
\item \code{y,ry,x}  are obtained from \code{mice()}

\item \code{type} is automatically extracted from the respective row of the \code{predictorMatrix}. Allowed entries in the \code{predictorMatrix} are: 
\itemize{
\item -2 = class variable (only one class variable is allowed!)
\item 0 = variable not included in imputation model
\item 1 = variable will be included as a fixed effect (zero \emph{and} count model)
\item 2 = variable will be included as a fixed \emph{and} random effect (zero \emph{and} count model)
\item 3 = variable will be included as a fixed effect (count model only)
\item 4 = variable will be included as a fixed \emph{and} random effect (count model only)
\item 5 = variable will be included as a fixed effect (zero model only)
\item 6 = variable will be included as a fixed \emph{and} random effect (zero model only)
}
}

The Bayesian regression variants (see Rubin 1987, p. 169-170) consist of the following steps:

1. Fit the zero model (a two-level binomial generalized linear mixed effects model), using the \code{glmmPQL} function from package \pkg{MASS} and find bhat, the posterior mean, and V(bhat), the posterior variance of model parameters b

2. Draw b* from N(bhat,V(bhat))

3. Compute predicted probabilities for having a zero vs. non-zero count

4. Draw imputations from a Binomial distribution and ``remember'' cases, who are supposed to get a non-zero count later on

5. Fit the count model (a zero-truncated 2L NB model) using the \code{glmmadmb} function from package \pkg{glmmADMB} and the \code{truncnbinom} family; find bhat, the posterior mean, and V(bhat), the posterior variance of model parameters b.

6. Draw b* from N(bhat,V(bhat))

7. Compute predicted counts and draw non-zero imputations (from step 4) from a zero-truncated NB distribution.

The bootstrap functions draw a bootstrap sample from \code{y[ry]} and \code{x[ry,]} and consist of the following steps:

1. Fit the zero model to the bootstrap sample

2. Compute predicted probabilities for having a zero vs. non-zero count

3. Draw imputations from a Binomial distribution and ``remember'' cases, who are supposed to get a non-zero count

4. Fit the count model to the boostrap sample

5. Compute predicted counts and draw non-zero imputations (from step 3) from a zero-truncated NB distribution.


}
\value{
\item{imp}{Numeric vector of length \code{sum(!ry)} with imputations}
}
\references{
Hilbe, J. M. (2007). \emph{Negative binomial regression}. Cambridge: Cambridge University Press.

Kleinke, K., & Reinecke, J. (2013). \emph{countimp 1.0 -- A multiple imputation package for incomplete count data} [Technical Report]. University of Bielefeld, Faculty of Sociology, available from \url{www.uni-bielefeld.de/soz/kds/pdf/countimp.pdf}.

Rubin, D. B. (1987). \emph{Multiple imputation for nonresponse in surveys}. New York: Wiley.

Zeileis, A., Kleiber, C., & Jackman, S. (2008). Regression models for count data in R. \emph{Journal of Statistical Software}, 27(8), 1--25.
}
\author{
Kristian Kleinke and Jost Reinecke
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data ( "MZINB.data" )
ini <- mice( MZINB.data, maxit = 0 )
pred <- ini$predictorMatrix
pred[1,] <- c( 0, 0, 6, 0, -2, 0, 1)
meth <- ini$method
meth[1] <- "2l.zihnb"
imp.x1.re.zero.only <- mice( MZINB.data, maxit=1, method=meth, pred=pred,
   seed=1234, m=5, print=FALSE)
result.x1.re.zero.only <- do.mira( imp.x1.re.zero.only, 
	DV = "Y", 
	fixedeff = c("X1+Z1", "Z1"), 
	randeff = c( "X1", "1"), fam = "truncnbinom", 
	grp = "GRP", id = "ID" )
summary(result.x1.re.zero.only)
}