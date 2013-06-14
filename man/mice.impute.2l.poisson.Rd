\name{mice.impute.2l.poisson}
\alias{mice.impute.2l.poisson}
\alias{mice.impute.2l.poisson.boot}
\alias{mice.impute.2l.poisson.noint}
\alias{mice.impute.2l.poisson.noint.boot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multiple Imputation of Incomplete Two-Level Poisson Data
}
\description{
The functions impute multilevel Poisson data based on a two-level generalized linear mixed effects Poisson model, either using a Bayesian regression approach (funxtions \code{mice.impute.2l.poisson} and \code{mice.impute.2l.poisson.noint}) or a bootstrap approach (functions with ``\code{.boot}'' in the function name). The \code{.noint} variants treat the intercept only as a fixed, but \emph{not} as a random effect. 
}
\usage{
mice.impute.2l.poisson(y, ry, x, type, intercept = TRUE)
mice.impute.2l.poisson.noint(y, ry, x, type, intercept = FALSE)
mice.impute.2l.poisson.boot(y, ry, x, type, intercept = TRUE)
mice.impute.2l.poisson.noint.boot(y, ry, x, type, intercept = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{y}{
Numeric vector with incomplete data in long format (i.e. the groups are stacked upon each other).
}
  \item{ry}{
Response pattern of \code{y} (\code{TRUE}=observed, \code{FALSE}=missing)
}
  \item{x}{
matrix with \code{length(y)} rows containing complete covariates, also in long format
}
  \item{type}{
vector of length \code{ncol(x)} determining the imputation model; \code{type} is automatically extracted from the \code{predictorMatrix} argument of \code{mice()}; see \pkg{mice}'s user's manual for details about how to specify the imputation model; see also the ``Details'' section. 
}
  \item{intercept}{\code{TRUE}: model will include intercept as a random effect; \code{FALSE}: model will not use intercept as a random effect.
%%     ~~Describe \code{intercept} here~~
}
}
\details{
The functions have been designed to multiply impute incomplete ``ordinary'' two-level count data (i.e. assuming that the data by and large follow the Poisson distribution) using a glmmPQL Poisson model (generalized linear mixed model with multivariate normal random effects, using Penalized Quasi-Likelihood; see Schall (1991) for details). We use the \code{glmmPQL()} function from package \pkg{MASS} to fit this model. 

Model specification details: 
\itemize{
\item \code{y,ry,x}  are obtained from \code{mice()}

\item \code{type} is automatically extracted from the respective row of the \code{predictorMatrix}. Allowed entries in the \code{predictorMatrix} are: \itemize{
\item 0 = variable not included in imputation model
\item 1 = variable will be included as a fixed effect
\item 2 = variable will be included both as a fixed \emph{and} random effect
\item -2 = class variable (only one class variable is allowed!)
}
}%%  ~~ If necessary, more details than the description above ~~

The Bayesian regression variants (see Rubin 1987, p. 169--170) consist of the following steps:

1. Fit the glmmPQL model, and find bhat, the posterior mean, and V(bhat), the posterior variance of model parameters b

2. Draw b* from N(bhat,V(bhat))

3. Compute predicted counts \code{p}

4. Draw imputations from Poisson distribution with mean \code{p}


The bootstrap functions draw a bootstrap sample from \code{y[ry]} and \code{x[ry,]} and consist of the following steps:

1. Fit the glmmPQL model to the bootstrap sample and get model parameters b

2. Compute predicted counts \code{p}

3. Draw imputations from Poisson distribution with mean \code{p}

}
\value{
\item{imp}{Numeric vector of length \code{sum(!ry)} with imputations}
}
\references{
Kleinke, K. & Reinecke, J. (2013a). \emph{Multiple Imputation of Multilevel Count Data}. Prepared for inclusion in ``Improving Survey Methods'', New York: Routledge. 

Kleinke, K., & Reinecke, J. (2013b). \emph{countimp 1.0 -- A multiple imputation package for incomplete count data} [Technical Report]. University of Bielefeld, Faculty of Sociology, available from \url{www.uni-bielefeld.de/soz/kds/pdf/countimp.pdf}.

Rubin, D. B. (1987). \emph{Multiple imputation for nonresponse in surveys}. New York: Wiley.

Schall, R. (1991). Estimation in generalized linear models with random effects. \emph{Biometrika}, 78, 719--727.
}
\author{
Kristian Kleinke and Jost Reinecke
}
\note{
The functions require packages \pkg{MASS} and \pkg{nlme}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data ( "NB.data" )

ini <- mice( NB.data, maxit = 0 )
pred <- ini$predictorMatrix
pred[1,] <-c ( 0, 0, 2, 0, -2, 0, 1 )
meth <- ini$method
meth[1] <- "2l.poisson"

imp <- mice( NB.data, maxit = 1, method = meth, pred = pred, 
seed = 1234, print = FALSE )

result<-do.mira( imp = imp, DV = "Y", 
fixedeff = "X1+Z1", randeff = "X1", 
grp = "GRP", id = "ID", fam = "poisson")
summary(result)
}