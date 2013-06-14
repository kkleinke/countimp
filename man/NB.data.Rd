\name{NB.data}
\alias{NB.data}
\docType{data}
\title{
Example Dataset Containing Simulated Incomplete Two-Level Negative Binomial Data
}
\description{
An example two-level dataset, containing an incomplete NB distributed count variable, a continuous individual level predictor, and a continuous group level predictor. The dataset is already in long format, i.e. the groups are stacked upon each other.
}
\usage{data("NB.data")}
\format{
  A data frame with 12604 observations from 50 groups and 7 variables:
  \describe{
    \item{\code{Y}}{incomplete count variale}
    \item{\code{X0}}{column of 1s, corresponding to individual level intercept}
    \item{\code{X1}}{continuous individual level predictor}
    \item{\code{ID}}{Participant identifier; ``ID'' variable}
    \item{\code{GRP}}{Group identifier; variable indicates cluster / group membership}
    \item{\code{Z0}}{column of 1s, corresponding to group level intercept}
    \item{\code{Z1}}{continuous group level predictor}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
simulated data by Kristian Kleinke and Jost Reinecke%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
Kleinke, K., & Reinecke, J. (2013b). \emph{countimp 1.0 -- A multiple imputation package for incomplete count data} [Technical Report]. University of Bielefeld, Faculty of Sociology, available from \url{www.uni-bielefeld.de/soz/kds/pdf/countimp.pdf}.

}
\examples{
data( "NB.data" )
head(NB.data)
}
\keyword{datasets}
