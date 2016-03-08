\name{centroidMethod}
\alias{centroidMethod}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate point estimates and standard errors.
}
\description{
Estimate point estimates and standard errors from multiple data points reported in primary literature. This is particularly useful for backcalculating odds ratios and confidence intervals for meta-analysis. Although simple back-calculation methods exist, they do not always create estimates that are consistent with the reported values. For instance, consider the example of an odds ratio of 1.3, with a 95\% confidence interval of 1.1 to 1.7; the back calculation methods do not forward calculate the correct decimal places (see \code{"Examples"}). The current method creates a potential area of legal values and iteratively tests a series of discrete points for compliance with the rounding constraints; the centroid of the legal values is calculated, returning the point estimate and standard error.
}
\usage{
centroidMethod(input, 
               decimals = NULL, 
               resolution = 250, 
               type = "OR", 
               distribution = "z", 
               df = NULL, 
               plot = FALSE, 
               ci = 95, 
               pInfo = NULL)
}
\arguments{
  \item{input}{
input data as c(PE,LCL,UCL), where PE is the point estimate as OR, log(OR),beta, or mean. Values can be entered as character or numeric. Character entry is particularly useful if there are trailing zeros and the number of decimal places is not specified in \code{"decimals"} (e.g., 1.80 will truncate to 1.8). Data can also be input as c(PE,LCL,UCL,df), where df is the degrees of freedom for a t-distribution (see \code{"df"} and \code{"distribution"})
}
  \item{decimals}{
specifies the number of decimal places for calculating rounding error. \code{"decimals"} can take one of three formats. 1) \code{"NULL"} (default), in which case the function itself will determine the number of decimals from the \code{"input"}. This can be problematic if \code{"input"} is formatted as numeric with trailing zeros, as they will be truncated (e.g., 1.80 will truncate to 1.8, and thus will be assumed to have 1 decimal place). 2) A single number to be used for each of the three numbers. This is useful if the number of decimal places is the same for each of the three \code{"input"} values. 3) A vector for each of the three \code{"input"} values; useful if the number of reported decimal places is different for each of the three \code{"input"} values.
}
  \item{resolution}{
a number that divides the possible standard error and point estimate spaces. The higher the number, the greater the precision of the centroid estimate and the greater the processing time. High resolution numbers will make \code{"plot"} take much longer because it plots each tested point one at a time.
}
  \item{type}{
specifies the type of \code{"input"} to process. Legal values are "OR", "log(OR)", "beta", or "mean". The log(OR), beta, and mean types use the same calculations because they do not need to be transformed, unlike ORs and 95\% CIs that must be log transformed and exponentiated during calculation.
}
  \item{distribution}{
select the distribution to be used to calculate and back-calculate the 95\% confidence intervals. Legal values are "z" and "t". If \code{"distribution"} == 't', \code{"df"} must be specified.
}
  \item{df}{
degrees of freedom (numeric). Required if \code{"distribution"} == 't'. \code{"df"} must be extracted or inferred from the original report.
}
  \item{plot}{
logical for whether to plot a figure that shows the possible standard error and point estimate ranges tested, the plane of values that satisfy all model criteria, and each of the four estimates (SE from the lower and upper confidence intervals, the average SE, and the current centroid method).
}
  \item{ci}{
an integer between 0 and 100 that specifies the confidence interval (e.g., 95\% confidence interval, which is the default, would be entered as 95).
}
  \item{pInfo}{
a vector of three values: the reported p-value, the number of decimal places, and the p-value function. For example, for a study that reported p<0.05, the data would be c(0.05,2,"<"); for a study that reported p=0.020, the data would be c(0.020,3,"==").
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The returned object is a data.frame with two columns: \code{"type"} and SE.
}
\references{
Brown and Dawson. Method for synthesizing multiple reported values to calculate point estimates and standard errors from reported point estimates, 95\% confidence intervals, and p-values. In progress.
}
\author{
Andrew W Brown
}

\examples{

centroidMethod(input = c(1.3,1.1,1.7), #Odds ratio of 1.3, 95\% CI: 1.1, 1.7
               decimals = NULL, #let the function determine decimal places
               resolution = 250, #grid of 250 x 250 points
               type = "OR", #Odds ratio
               distribution = "z", #assume z distribution
               df = NULL, #z distribution, so no df
               plot = TRUE, #plot the output
               ci = 95, #95\% CI
               p = NULL #No p-value reported
               ) 
               
centroidMethod(input = c(1.3,1.1,1.7), #Odds ratio of 1.3, 95\% CI: 1.1, 1.7
               decimals = NULL, #let the function determine decimal places
               resolution = 250, #grid of 250 x 250 points
               type = "OR", #Odds ratio
               distribution = "z", #assume z distribution
               df = NULL, #z distribution, so no df
               plot = TRUE, #plot the output
               ci = 95, #95\% CI
               p = c(0.02,2,"<") #p-value reported as range
               )  
               
centroidMethod(input = c(1.3,1.1,1.7), #Odds ratio of 1.3, 95\% CI: 1.1, 1.7
               decimals = NULL, #let the function determine decimal places
               resolution = 250, #grid of 250 x 250 points
               type = "OR", #Odds ratio
               distribution = "z", #assume z distribution
               df = NULL, #z distribution, so no df
               plot = TRUE, #plot the output
               ci = 95, #95\% CI
               p = c(0.011,3,"==") #p-value reported rounded to 3 decimals
               )       
               
}

\keyword{ ~Meta-analysis }
