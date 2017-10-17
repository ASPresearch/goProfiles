\name{compareGOProfiles}
\alias{compareGOProfiles}
\title{Comparison of lists of genes through their functional profiles}
\description{
Compare two samples of genes in terms of their GO profiles \code{pn} and \code{qm}. Both
samples may share a common subsample of genes, with GO profile \code{pqn0}.
'compareGOProfiles' implements some inferential procedures based on
asymptotic properties of the squared euclidean distance between 
the contracted versions of pn and qm
}
\usage{                                
compareGOProfiles(pn, qm = NULL, pqn0 = NULL, n = ngenes(pn), m = ngenes(qm), n0 = ngenes(pqn0), method = "lcombChisq", 
ab.approx = "asymptotic", confidence = 0.95, nsims = 10000, simplify = T, ...)
}
\arguments{
\item{pn}{an object of class ExpandedGOProfile representing one or more
"sample" expanded GO profiles for a fixed ontology (see the 'Details' section)}
  \item{qm}{an object of class ExpandedGOProfile representing one or more
"sample" expanded GO profiles for a fixed ontology (see the 'Details' section)}
\item{pqn0}{an object of class ExpandedGOProfile representing one or more
"sample" expanded GO profiles for a fixed ontology (see the 'Details' section)}
\item{n}{a numeric vector with the number of genes profiled in each column of pn.
This parameter is included to allow the possibility of exploring the consequences of varying sample sizes, other than the true sample size in pn.}
\item{m}{a numeric vector with the number of genes profiled in each column of qm.}
\item{n0}{a numeric vector with the number of genes profiled in each column of pqn0.}
\item{method}{the approximation method to the sampling distribution under the null hypothesis specifying that the samples pn and qm come from the same population. See the 'Details' section below}
\item{confidence}{the confidence level of the confidence interval in the result}
\item{ab.approx}{the approximation used for computing 'a' and 'b' coefficients (see details)}
\item{nsims}{some inferential methods require a simulation step; the number of simulation replicates is specified with this parameter}
\item{simplify}{should the result be simplified, if possible? See the 'Details' section}
\item{...}{Other arguments needed}
}

\details{
 An object of S3 class 'ExpandedGOProfile' is, essentially, a 'data.frame' object
 with each column representing the relative frequencies in all observed node
 combinations, resulting from profiling a set of genes, for a given and fixed
 ontology. The row.names attribute codifies the node combinations and each
 data.frame column (say, each profile) has an attribute, 'ngenes', indicating the
 number of profiled genes.
 The arguments 'pn', 'qm' and 'pqn0' are compared in a column by column wise, 
 recycling columns, if necessary, in order to perform max(ncol(pn),ncol(qm),ncol(pqn0))
 comparisons (each comparison resulting in an object of class 'GOProfileHtest',
 an specialization of 'htest').  
 In order to be properly compared, these arguments are expanded by row, according
 to their row names. That is, the data arguments can have unequal row numbers. Then,
 they are expanded adding rows with zero frequencies, in order to make them
 comparable.

 In the i-th comparison (i from 1 to max(ncol(pn),ncol(qm),ncol(pqn0))),
 the parameters n, m and n0 are included to allow the possibility of exploring the consequences of varying sample
 sizes, other than the true sample sizes included as an attribute in pn, qm and pqn0.

 When qm = NULL, the genes profiled in pn are compared with a subsample of them, those profiled in pqn0 (compare
 a set of genes with a restricted subset, e.g. those overexpressed under a disease). In this case we take qm=pqn0.
 When pqn0 = NULL, two profiles with no genes in common are compared.

 Let Pn and Qm correspond to the contracted functional profiles (the total counts or relative frequencies
 of hits in each one of the s GO categories being compared) obtained from pn and qm.
 If P stands for the "population" profile originating the sample profile Pn[,j], Q for the profile originating Qm[,j]
 and d(,) for the squared euclidean distance, if P != Q, the distribution of sqrt(nm/(n+m))(d(Pn[,j],Qm[,j]) - d(P,Q))/se(d) is
 approximately standard normal, N(0,1). This provides the basis for the confidence interval in the result field
 icDistance.
 When P=Q, the asymptotic distribution of (nm/(n+m)) d(Pn[,j],Qm[,j]) corresponds to the distribution of a mixture of independent
 chi-square random variables, each one with one degree of freedom. The sampling distribution under H0 P=Q may be directly
 computed from this distribution (approximating it by simulation) (method="lcombChisq") or by a chi-square
 approximation to it, based on two correcting constants a and b (method="chi-square").
 These constants are chosen to equate the first two moments of both distributions (the linear
 combination of chi-square random variables distribution and the approximating chi-square distribution).
 When method="chi-square", the returned test statistic value is the chi-square approximation (n d(pn[,j],qm[,j]) - b) / a. Then,
 the result field 'parameter' is a vector containing the 'a' and 'b' values and the number of degrees of freedom, 'df'.
 Otherwise, the returned test statistic value is (nm/(n+m)) d(Pn[,j],Qm[,j]) and 'parameter' contains the coefficients of the linear
 combination of chi-squares.
}
\value{
A list containing max(ncol(pn),ncol(qm),ncol(pqn0)) objects of class 'GOProfileHtest', directly inheriting from 'htest'
or a single 'GOProfileHtest' object if max(ncol(pn),ncol(qm),ncol(pqn0))==1 and simplify == T.
Each object of class 'GOProfileHtest' has the following fields:
\item{profilePn}{the first contracted profile to compute the squared Euclidean distance}
\item{profileQm}{the second contracted profile to compute the squared Euclidean distance}
\item{statistic}{test statistic; its meaning depends on the value of "method", see the 'Details' section.}
\item{parameter}{parameters of the sample distribution of the test statistic, see the 'Details' section.}
\item{p.value}{associated p-value to test the null hypothesis of profiles equality.}
\item{conf.int}{asymptotic confidence interval for the squared euclidean distance. Its attribute "conf.level" contains its nominal confidence level.}
\item{estimate}{squared euclidean distance between the contracted profiles. Its attribute "se"
contains its standard error estimate.}
\item{method}{a character string indicating the method used to perform the test.}
\item{data.name}{a character string giving the names of the data.}
\item{alternative}{a character string describing the alternative hypothesis (always 'true squared Euclidean distance between the contracted profiles is greater than zero'}
}
\seealso{fitGOProfile, equivalentGOProfiles}
\references{Sanchez-Pla, A., Salicru M. and Ocana, J. Statistical methods for the analysis of highthroughput data based on functional profiles derived from the gene ontology. Journal of
Statistical Planning and Inference, 2007.}
\author{Jordi Ocana}

\examples{
# [NOT RUN COMPLETELY]
data(prostateIds)
expandedWelsh <- expandedProfile(welsh01EntrezIDs[1:100], onto="MF",
                       level=2, orgPackage="org.Hs.eg.db")
expandedSingh <- expandedProfile(singh01EntrezIDs[1:100], onto="MF",
                       level=2, orgPackage="org.Hs.eg.db")
commonGenes <- intersect(welsh01EntrezIDs[1:100], singh01EntrezIDs[1:100])
commonExpanded <- expandedProfile(commonGenes, onto="MF", level=2, orgPackage="org.Hs.eg.db")
# comparedMF <-compareGOProfiles (pn=expandedWelsh, 
#                          qm  = expandedSingh, 
#                          pqn0= commonExpanded)
# print(comparedMF)
# print(compSummary(comparedMF))
#
}
\keyword{htest}
