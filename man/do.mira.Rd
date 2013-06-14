\name{do.mira}
\alias{do.mira}
\alias{summary.do.mira}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Automated Repeated Data Analysis for two-level GLMMs
}
\description{
The function automates repeated data analysis and MI inference for two-level GLMMS. The generic summary function extracts the following output: (a) the pooled fixed effects coefficients (produced by \code{miinference()}), (b)
the pooled random effects standard deviations (mean of the \code{m} SDs), (c) the pooled random effects correlations (mean of the \code{m} correlations)
.
}
\usage{
do.mira(imp,DV,fixedeff,randeff,grp,id,fam="poisson")
\method{summary}{do.mira}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{imp}{
an object of class \code{mids} created by \code{mice()} using any of our two-level imputation functions} 

\item{DV}{string, dependend variable in the two-level count model}

\item{fixedeff}{
string of length 1 or 2 of fixed effects in R formula language, see details
}
\item{randeff}{
string of length 1 or 2 of random effects in R formula language, see details
}
\item{grp}{class variable; string}
\item{id}{participant identifier; string}
\item{fam}{
GLMM family, default is \code{"poisson"} for a Poisson GLMM; other possible entries are: \code{"nbinom"} for a negative binomial two-level model, or \code{"truncnbinom"} for a two-level zero-inflated hurdle NB model.}
\item{object}{an object of class \code{"do.mira"}}
\item{\dots}{ not used}
}

\details{
The \code{fam} argument determines, which model is fit. 

\code{fam="poisson"} fits a two-level generalized linear mixed effects model using the \code{glmmPQL} function from package \pkg{MASS}. 

\code{fam="negbin"} fits a two-level generalized linear mixed effects NB model using first \code{glmmadmb()} from package \pkg{glmmADMB} to get an estimate for the dispersion parameter of the NB model, then fitting the NB model using \code{glmmPQL()} with family \code{negative.binomial(theta)} and the previously estimated dispersion parameter as input. We use this two-step approach, as -- currently -- \code{glmmadmb()} does not produce estimates of random effects covariances. 

\code{fam="truncnbinom"} fits a two-level zero-inflated hurdle NB model using the \code{glmmadmb()} function. At the moment, the function does not produce estimates of random effect covariances. However, it is being worked on that issue and we hope that this function will be available with the next update. If \code{fam="truncnbinom"} is selected, \code{fixedeff} and \code{randeff} can take a string vector with two entries: the first with the model formula for the zero model, the second entry with the model formula for the count model. See the example for details.
}
\value{
a list containing the following elements:
\item{FIT}{the fitted model objects}
\item{COEF.RA}{fixed effects of the repeated analyses}
\item{SE.RA}{standard errors of the coefficients in \code{COEF.RA}}
\item{MIINFERENCE}{Results of the MI inference produced by \code{miinference()} for the fixed effects} 
\item{R.SD}{Random effects standard deviation(s)}
\item{R.COR}{Random effects correlation(s), list with \code{imp$m} elements}
\item{PR.SD}{Pooled random effects standard deviation(s); the mean of the \code{imp$m} estimates}
\item{PR.COR}{Pooled random effects correlation(s); the mean of the \code{imp$m} estimates}
}
\references{
Kleinke, K., & Reinecke, J. (2013). \emph{countimp 1.0 -- A multiple imputation package for incomplete count data} [Technical Report]. University of Bielefeld, Faculty of Sociology, available from \url{www.uni-bielefeld.de/soz/kds/pdf/countimp.pdf}.
}

\author{
Kristian Kleinke and Jost Reinecke
}
\note{

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