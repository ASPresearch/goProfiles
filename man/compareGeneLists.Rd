\name{compareGeneLists}
\alias{compareGeneLists}
\title{Compares two lists of genes by building (expanded) profiles and comparing them}
\description{
This function wraps all the needed steps to compare two lists of genes following 
the methodology developed by Sanchez, Salicru and Ocan\~a (2007)
}
\usage{
compareGeneLists(genelist1, genelist2, idType = "Entrez", onto = "ANY",
level = 2, orgPackage, 
method = "lcombChisq", ab.approx = "asymptotic", confidence = 0.95, compareFunction="compareGOProfiles", ...)
}
\arguments{
  \item{genelist1}{First gene set to be compared}
  \item{genelist2}{Second gene set to be compared}
  \item{idType}{Type of identifiers for the genes. May be 'Entrez' (default), BiocProbes or GoTermsFrame. See the 'Details' section below}
  \item{onto}{Ontology on which the profile has to be built}
  \item{level}{Level of the ontology at which the profile has to be
    built}
  \item{orgPackage}{Name of a Bioconductor's organism annotations package
    ('org.Xx-eg-db')}
\item{method}{The approximation method to the sampling distribution under the null hypothesis specifying that the samples pn and qm come from the same population. See the 'Details' section below}
\item{confidence}{The confidence level of the confidence interval in the result}
\item{ab.approx}{The approximation used for computing 'a' and 'b' coefficients (see details)}
\item{compareFunction}{Allows to use 'fitGOProfile' (sample vs population) or 'compareGOProfiles' (sample1 vs sample2)}
  \item{\dots}{Other arguments for the methods 'basicProfile' or 'compareGoProfiles'}
}
\value{
The result of the comparison is a list with a variable number of arguments, depending for which ontologies has been performed the comparison. Each list member is an object of class 'htest' corresponding to the output of the function \code{compareGOProfiles}}
\references{Sanchez-Pla, A., Salicru, M. and Ocana, J.
Statistical methods for the analysis of high-throughput data based on
functional profiles derived from the gene ontology.
Journal of Statistical Planning and Inference, 2007}
\author{Alex Sanchez}
\seealso{\code{\link{compareGOProfiles}}, \code{\link{basicProfile}} }
\examples{
data(prostateIds)
prostateCompared<- compareGeneLists (welsh01EntrezIDs[1:500],
singh01EntrezIDs[1:500], level=2, onto='MF', orgPackage="org.Hs.eg.db")
print(prostateCompared)
# print(compSummary(prostateCompared))
}
\keyword{htest}
