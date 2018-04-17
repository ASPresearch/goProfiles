\name{compareProfilesLists}
\alias{compareProfilesLists}
\title{Compares two of expanded profiles}
\description{
This function compares two lists (``sensu R lists'') of expanded profiles by successive 
calls to function \code{compareGOProfiles} following 
the methodology developed by Sanchez, Salicru and Ocan\~a (2007)
}
\usage{
compareProfilesLists(expanded1, expanded2, common.expanded=NULL, relationType,
method = "lcombChisq", ab.approx = "asymptotic", confidence = 0.95, ...)
}
\arguments{
  \item{expanded1}{First expanded profile to be compared}
  \item{expanded2}{Second expanded profile to to be compared}
  \item{common.expanded}{Expanded profile made from the genes appearing in both lists of genes}
  \item{relationType}{Type of relation between gene lists compared through the expanded profiles. It can be INCLUSION, INTERSECTION or DISJOINT}
\item{method}{The approximation method to the sampling distribution under the null hypothesis specifying that the samples pn and qm come from the same population. See the 'Details' section below}
\item{confidence}{The confidence level of the confidence interval in the result}
\item{ab.approx}{The approximation used for computing 'a' and 'b' coefficients (see details)}
  \item{\dots}{Other arguments for the methods 'basicProfile' or 'compareGoProfiles'}
}
\value{
The result of the comparison is a list with a variable number of arguments, depending for which ontologies has been performed the comparison. Each list member is an object of class 'htest' corresponding to the output of the function \code{compareGOProfiles}}
\references{Sanchez-Pla, A., Salicru, M. and Ocana, J.
Statistical methods for the analysis of high-throughput data based on
functional profiles derived from the gene ontology.
Journal of Statistical Planning and Inference, 2007}
\author{Alex Sanchez}
\seealso{\code{\link{compareGeneLists}}, \code{\link{expandedProfile}} }
\examples{
#[NOT RUN] 
#data(ProstateIds)
#expandedWelsh <- expandedProfile(welsh01EntrezIDs[1:100], onto="MF",
#                        level=2, orgPackage="org.Hs.eg.db")
#expandedSingh <- expandedProfile(singh01EntrezIDs[1:100], onto="MF",
#                        level=2, orgPackage="org.Hs.eg.db")
#commonGenes <- intersect(welsh01EntrezIDs[1:100], singh01EntrezIDs[1:100])
#commonExpanded <- expandedProfile(commonGenes, onto="MF", level=2, orgPackage="org.Hs.eg.db")
#comparedMF<- compareProfilesLists (expandedWelsh, expandedSingh, commonExpanded, 
#                                   relationType="COMMON")
#print(comparedMF)
#print(compSummary(comparedMF))
}
\keyword{htest}
