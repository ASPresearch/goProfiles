\name{fitGOProfile}
\alias{fitGOProfile}
\title{Does a "sample" GO profile 'pn', observed in a sample of 'n' genes, fit a "population" or "model" p0?}
\description{
'fitGOProfile' implements some inferential procedures to solve the preceding question.
These procedures are based on asymptotic properties of the squared euclidean distance
between the contracted versions of pn and p0
}
\usage{
fitGOProfile(pn, p0, n = ngenes(pn), method = "lcombChisq", 
ab.approx = "asymptotic", confidence = 0.95, simplify = T)
}
\arguments{
  \item{pn}{an object of class ExpandedGOProfile representing one or more "sample" expanded GO profiles for a fixed ontology (see the 'Details' section)}
  \item{p0}{  an object of class ExpandedGOProfile representing one or more "population" or "theoretical" expanded GO profiles (see also the 'Details' section)  }
  \item{n}{ a numeric vector with the number of genes profiled in each column of pn. This parameter is included to allow the possibility of exploring the consequences of varying sample sizes, other than the true sample size in pn}
  \item{method}{the approximation method to the sampling distribution under the null hypothesis "p = p0", where p is the 'true' population profile originating each column of pn. See the 'Details' section below }
  \item{ab.approx}{the method used to compute the constants 'a' and 'b' described in the paper. See the 'Details' section}
  \item{confidence}{the confidence level of the confidence interval in the result}
  \item{simplify}{ should the result be simplified, if possible? See the 'Details' section}
}
\details{An object of class 'ExpandedGOProfile' is, essentially, a 'data.frame' object
	 with each column representing the relative frequencies in all observed node
	 combinations, resulting from profiling a set of genes, for a given and fixed
	 ontology. The row.names attribute codifies the node combinations and each
	 data.frame column (say, each profile) has an attribute, 'ngenes', indicating the
	 number of profiled genes. (Actually, the 'ngenes' attribute of each 'p0' column
	 is ignored and is taken as if it were infinite, 'Inf'.)
	 The arguments 'pn' and 'p0' are compared in a column by column wise, 
	 recycling columns, if necessary, in order to perform max(ncol(pn),ncol(p0))
	 comparisons (each comparison resulting in an object of class 'htest').  
	 In order to be properly compared, 'pn' and 'p0' are expanded by row, according
	 to their row names. That is, both arguments can have unequal row numbers. Then,
	 they are expanded adding rows with zero frequencies, in order to make them
	 comparable.

	 In the i-th comparison (i from 1 to max(ncol(pn),ncol(p0))),
	 if p stands for the profile originating the sample profile pn[,i] and d(,)
	 for the squared euclidean distance, 
	 if p != p0[,i], the distribution of
	 sqrt(n)(d(pn[,i],p0[,i]) - d(p,p0[,i]))/se
	 is approximately standard normal, N(0,1). This provides
	 the basis for the confidence interval in the result field conf.int.
	 When p==p0[,i], the asymptotic distribution of
	 n d(pn[,i],p0[,i]) is the distribution of a linear combination of independent
	 chi-square random variables, each one with one degree of freedom.
	 This sampling distribution may be directly computed (approximating it by simulation, method="lcombChisq")
	 or approximated by a chi-square distribution, based on two correcting constants a and b (method="chi-square").
	 These constants are chosen to equate the first two moments of both distributions (the distribution of a linear combination
	 of chi square variables and the approximating chi-square distribution).
	 When method="chi-square", the returned test statistic value is the chi-square approximation (n d(pn,p0) - b) / a. Then,
	 the result field 'parameter' is a vector containing the 'a' and 'b' values and the number of degrees of freedom, 'df'.
	 Otherwise, the returned test statistic value is n d(pn,p0) and 'parameter' contains the coefficients of the linear
	 combination of chi-squares
}
\value{A list containing max(ncol(pn),ncol(p0)) objects of class 'htest', 
or a single 'htest' object if ncol(pn)==1 and ncol(p0)==1 and simplify == T.
Each 'htest' object has the following fields:
\item{statistic}{test statistic; its meaning depends on the value of "method", see the 'Details' section}
\item{parameter}{parameters of the sample distribution of the test statistic, see the 'Details' section}
\item{p.value}{associated p-value to test the null hypothesis "pn[,i] is a random sample taken from p0[,i]"}
\item{conf.int}{asymptotic confidence interval for the squared euclidean distance. Its attribute "conf.level" contains its nominal confidence level}
\item{estimate}{squared euclidean distance between the contracted pn and p0 profiles. Its attribute "se" contains its standard error estimate}
\item{method}{a character string indicating the method used to perform the test}
\item{data.name}{a character string giving the names of the data}
\item{alternative}{a character string describing the alternative hypothesis}
}
\references{Sanchez-Pla, A., Salicru, M. and Ocana, J.
Statistical methods for the analysis of high-throughput data based on
functional profiles derived from the gene ontology.
Journal of Statistical Planning and Inference, 2007.
}
\author{Jordi Ocana}
\seealso{compareGOProfiles}
\examples{
#data(sampleProfiles)
#comparedMF <-fitGOProfile(pn=expandedWelsh01[['MF']], 
#                          p0  = expandedSingh01[['MF']])
#print(comparedMF)
#print(compSummary(comparedMF))
}
\keyword{htest}
