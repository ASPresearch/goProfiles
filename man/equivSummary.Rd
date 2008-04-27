\name{equivSummary}
\alias{equivSummary}
\title{This function returns a brief summary of the equivalence test between two profiles.}
\description{
  Function to return a brief summary of the equivalence test between two profiles.
  If In its current version it is better that \code{equivalentGOProfiles} is called with option
\code{simplify} set to \code{FALSE} before \code{equivSummary} can be used}
\usage{
equivSummary(l, decs = 6)
}
\arguments{
  \item{l}{A list of comparison results as returned by a call to \texttt{compareGenelists}}
  \item{decs}{Number of decimal places to use in the output}
}
\value{
A data frame with the summarized results of each comparison.
The values contained are: \texttt{Sqr.Eucl.Dist}: The squared euclidean distance, 
\texttt{Standard Err}: The standard error estimate, \texttt{pValue} p value of the equivalence test, 
\texttt{up conf.int}Upper value for the desired condfidence interval.
\texttt{d0}Threshold value for equivalence test.
\texttt{Equivalent?}Numerical value set to 1 if profiles can be considered equivalent and to zero if they cannot.
}
\author{Alex Sanchez}
\examples{
data(sampleProfiles)
comparedMF <-compareGOProfiles (pn=expandedWelsh01[['MF']], 
                          qm  = expandedSingh01[['MF']], 
                          pqn0= commonExpandedWelsh01Singh01[['MF']])
equivMF<- equivalentGOProfiles(comparedMF)
print(equivSummary(equivMF, decs=5))
}
\keyword{htest}

