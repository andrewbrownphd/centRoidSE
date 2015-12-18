\name{centRoidSE-package}
\alias{centRoidSE-package}
\docType{package}
\title{
Method to calculate standard errors and effect sizes from multiple reported values.
}
\description{
This package contains one function, \code{"centroidMethod"}, that uses multiple values reported in original research to get consistent, unbiased estimates of values of point estimates and standard errors for meta-analysis.
}
\details{
\tabular{ll}{
Package: \tab centRoidSE\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-12-18\cr
License: \tab MIT\cr
}
}
\author{
Andrew W Brown; John A Dawson

Maintainer: Andrew W Brown<awbrown@uab.edu>
}
\references{
Brown and Dawson. Method for synthesizing multiple reported values to calculate point estimates and standard errors from reported point estimates, 95\% confidence intervals, and p-values. In progress.}

\keyword{ package }
\keyword{ meta-analysis }
\seealso{
\code{\link[centroidMethod:centRoidSE-package]{centroidMethod}}
}
\examples{
centroidMethod(input = c(1.3,1.1,1.7), #Odds ratio of 1.3, 95\% CI: 1.1, 1.7
               decimals = NULL, #let the function determine decimal places
               resolution = 250, #grid of 250 x 250 points
               type = "OR", #Odds ratio
               distribution = "z", #assume z distribution
               df = NULL, #z distribution, so no df
               plot = T, #plot the output
               ci = 95, #95\% CI
               p = NULL #No p-value reported
               ) 
}
