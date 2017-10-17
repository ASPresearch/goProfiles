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
  \item{l}{A list of comparison results as returned by a call to \code{compareGenelists}}
  \item{decs}{Number of decimal places to use in the output}
}
\value{
A data frame with the summarized results of each comparison.
The values contained are: \code{Sqr.Eucl.Dist}: The squared euclidean distance, 
\code{Standard Err}: The standard error estimate, \code{pValue} p value of the equivalence test, 
\code{up conf.int}Upper value for the desired condfidence interval.
\code{d0}Threshold value for equivalence test.
\code{Equivalent?}Numerical value set to 1 if profiles can be considered equivalent and to zero if they cannot.
}
\author{Alex Sanchez}

\seealso{'equivalentGOProfiles'}

\examples{
# data(prostateIds)
# expandedWelsh <- expandedProfile(welsh01EntrezIDs[1:100], onto="MF",
#                        level=2, orgPackage="org.Hs.eg.db")
# expandedSingh <- expandedProfile(singh01EntrezIDs[1:100], onto="MF",
#                        level=2, orgPackage="org.Hs.eg.db")
#commonGenes <- intersect(welsh01EntrezIDs[1:100], singh01EntrezIDs[1:100])
#commonExpanded <- expandedProfile(commonGenes, onto="MF", level=2, orgPackage="org.Hs.eg.db")

# equivMF <-equivalentGOProfiles (pn=expandedWelsh, 
#                          qm  = expandedSingh, 
#                          pqn0= commonExpanded)
#print(equivSummary(equivMF, decs=5))
}
\keyword{htest}

