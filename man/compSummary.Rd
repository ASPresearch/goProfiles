\name{compSummary}
\alias{compSummary}
\title{This function returns a brief summary of the comparison between two (expanded) profiles.}
\description{
  Function to return a brief summary of the comparison between two (expanded) profiles.}
\usage{
compSummary(l, decs = 6)
}
\arguments{
  \item{l}{A list of comparison results as returned by a call to \texttt{compareGenelists}}
  \item{decs}{Number of decimal places to use in the output}
}
\value{
A data frame with the summarized results of each comparison.
The values contained are: \texttt{Sqr.Eucl.Dist}: The squared euclidean distance, 
\texttt{Standard Err}: The standard error estimate, \texttt{pValue} p value of the test, 
\texttt{low conf.int}Lower value for the desired confidence interval, 
\texttt{up conf.int}Upper value for the desired condfidence interval.
}
\author{Alex Sanchez}
\examples{
data(sampleProfiles)
comparedMF <-compareGOProfiles (pn=expandedWelsh01[['MF']], 
                          qm  = expandedSingh01[['MF']], 
                          pqn0= commonExpandedWelsh01Singh01[['MF']])
print(comparedMF)
print(compSummary(comparedMF))
}
\keyword{htest}

