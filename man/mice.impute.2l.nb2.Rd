\name{mice.impute.2l.nb2}
\alias{mice.impute.2l.nb2}
\alias{mice.impute.2l.nb2.boot}
\alias{mice.impute.2l.nb2.noint}
\alias{mice.impute.2l.nb2.noint.boot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multiple Imputation of Incomplete Overdispersed Two-Level Count Data
}
\description{
The functions impute overdispersed two-level count data based on a generalized linear mixed effects model (a NB2 model, see Hilbe, 2007, for details), either using a Bayesian regression (``\code{.nb2}'') or a bootstrap approach (``\code{.boot}''). The \code{.noint} variants treat the intercept only as a fixed, but \emph{not} as a random effect. 

}
\usage{
mice.impute.2l.nb2(y, ry, x, type, intercept = TRUE)
mice.impute.2l.nb2.boot(y, ry, x, type, intercept = TRUE)
mice.impute.2l.nb2.noint(y, ry, x, type, intercept = FALSE)
mice.impute.2l.nb2.noint.boot(y, ry, x, type, intercept = FALSE)
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
  \item{intercept}{\code{TRUE}: model will include intercept as a random effect; \code{FALSE}: model will not use intercept as a random effect.
%%     ~~Describe \code{intercept} here~~
}
}

\details{
The functions have been designed to multiply impute incomplete overdispersed two-level count data (assuming that the data by and large follow the Negative Binomial distribution). Overdispersion means that the variance of the data is larger in comparison to the mean -- which violates the assumption of equi-dispersion of the classical Poisson model. To get an estimate of the dispersion parameter, we use the \code{glmmadmb()} function of package \pkg{glmmADMB}. Note that \pkg{glmmADMB} is not on CRAN! The package can be obtained from \url{http://glmmadmb.r-forge.r-project.org}. Note also that the current version of \code{glmmadmb()} does not produce estimates of random effects covariances. Until this issue is addressed, we use \code{glmmadmb()} only to get an estimate for the dispersion parameter \eqn{\alpha} and fit the NB2 model using the \code{glmmPQL} function from package \pkg{MASS} with family \code{negative.binomial(}\eqn{\alpha}\code{)}. Unfortunately, we cannot use \code{glmmPQL()} directly to estimate the whole model, as the \code{negative.binomial} family in \pkg{MASS} requires an estimate of the dispersion parameter as input.

Model specification details: 
\itemize{
\item \code{y,ry,x}  are obtained from \code{mice()}

\item \code{type} is automatically extracted from the respective row of the \code{predictorMatrix}. Allowed entries in the \code{predictorMatrix} are: 
\itemize{
\item 0 = variable not included in imputation model
\item 1 = variable will be included as a fixed effect
\item 2 = variable will be included both as a fixed \emph{and} random effect
\item -2 = class variable (only one class variable is allowed!)
}
}%%  ~~ If necessary, more details than the description above ~~

The Bayesian regression variants (see Rubin 1987, p. 169--170) consist of the following steps:

1. Fit the NB2 model (first using \code{glmmadmb()} to get an estimate of \eqn{\alpha}, then using \code{glmmPQL()}, as described above), and find bhat, the posterior mean, and V(bhat), the posterior variance of model parameters b

2. Draw b* from N(bhat,V(bhat))

3. Compute predicted counts \code{p}

4. Draw imputations from Negative Binomial distribution with mean \code{p} and dispersion \eqn{\alpha}


The bootstrap functions draw a bootstrap sample from \code{y[ry]} and \code{x[ry,]} and consist of the following steps:

1. Fit the model to the bootstrap sample and get model parameters b

2. Compute predicted counts \code{p}

3. Draw imputations from Negative Binomial distribution with mean \code{p} and dispersion \eqn{\alpha}

}
\value{
\item{imp}{Numeric vector of length \code{sum(!ry)} with imputations}
}
\references{
Hilbe, J. M. (2007). \emph{Negative binomial regression}. Cambridge: Cambridge University Press.

Kleinke, K., & Reinecke, J. (2013). \emph{countimp 1.0 -- A multiple imputation package for incomplete count data} [Technical Report]. University of Bielefeld, Faculty of Sociology, available from \url{www.uni-bielefeld.de/soz/kds/pdf/countimp.pdf}.

Rubin, D. B. (1987). \emph{Multiple imputation for nonresponse in surveys}. New York: Wiley.
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
data( "NB.data" )
head( NB.data )

ini <- mice( NB.data, maxit = 0 )
pred <- ini$predictorMatrix
pred[1,] <-c ( 0, 0, 2, 0, -2, 0, 1 )
meth <- ini$method
meth[1] <- "2l.nb2.boot"


# impute missing data
# note: imputation may take a while ...
imp <- mice( NB.data, maxit = 1, method = meth, pred = pred,
        print = FALSE, seed = 1234 )

# do automated repeated data analysis and pooling
result<-do.mira( imp = imp, DV = "Y", 
        fixedeff = "X1+Z1", randeff = "X1", 
        grp = "GRP", id = "ID", fam = "nbinom")
summary(result)
}