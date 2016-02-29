\name{equivalentGOProfiles}
\alias{equivalentGOProfiles}
\alias{equivalentGOProfiles.GOProfileHtest}
\alias{equivalentGOProfiles.ExpandedGOProfile}
\alias{equivalentGOProfiles.default}
\title{Are two lists of genes equivalent in terms of their Gene Ontology profiles?}
\description{
  Performs an equivalence test based on the squared Euclidean distance
  between the Gene Ontology profiles of two lists of genes.
  Equivalence is declared if the upper limit d.sup of a one-sided confidence 
  interval [0, d.sup] for the distance is lesser than the equivalence limit d0.
}

\usage{
equivalentGOProfiles(goObject, \ldots)
\method{equivalentGOProfiles}{GOProfileHtest}(goObject, equivEpsilon = 0.05, d0 = NULL, confidence = NULL, \ldots)
\method{equivalentGOProfiles}{ExpandedGOProfile}(goObject, qm=NULL, pqn0=NULL,
    n = ngenes(goObject), m = ngenes(qm), n0 = ngenes(pqn0),
    confidence = 0.95,
    equivEpsilon = 0.05, d0 = NULL,
    simplify = FALSE, \ldots)
\method{equivalentGOProfiles}{default}(goObject, \ldots)
}

\arguments{
  \item{goObject}{ an object related to GO profiles or comparisons between them }
  \item{qm}{ an expanded GO profile, i.e. and object of class 'ExpandedGOProfile' }
  \item{pqn0}{ an expanded GO profile, i.e. and object of class 'ExpandedGOProfile' }
  \item{n}{a numeric vector with the number of genes profiled in each column of goObject.
  This parameter is included to allow the possibility of exploring the consequences of varying sample sizes, other than the true sample size in goObject.}
  \item{m}{a numeric vector with the number of genes profiled in each column of qm.}
  \item{n0}{a numeric vector with the number of genes profiled in each column of pqn0.}
  \item{confidence}{ the nominal confidence level of the one-sided confidence interval on the distance }
  \item{d0}{ a positive value specifying the equivalence limit }
  \item{equivEpsilon}{ a positive value used to compute 'd0' if it is not directly available }                             
  \item{simplify}{ should the result be simplified, if possible? See the 'Details' section }
  \item{\dots}{ further arguments, tipically the same than to 'compareGOProfiles' }
}
\details{
  An object of S3 class "ExpandedGOProfile" is, essentially, a "data.frame" object
  with each column representing the relative frequencies in all observed node
  combinations, resulting from profiling a set of genes, for a given and fixed
  ontology. The 'row.names' attribute codifies the node combinations and each
  "data.frame" column (say, each profile) has an attribute, 'ngenes', indicating the
  number of profiled genes. 
  
  In the 'ExpandedGOProfile' interface, the arguments 'goObject', 'qm' and 'pqn0' are compared in a column
  by column wise, recycling columns, if necessary, in order to perform max(ncol(goObject),ncol(qm),ncol(pqn0))
  equivalence tests (each test resulting in an object of class 'htest').  
  In order to be properly tested, these arguments are expanded by row, according
  to their row names. That is, the data arguments can have unequal row numbers. Then,
  they are expanded adding rows with zero frequencies, in order to make them
  comparable. In the i-th comparison (i from 1 to max(ncol(goObject),ncol(qm),ncol(pqn0))),
  the parameters n, m and n0 are included to allow the possibility of exploring the 
  consequences of varying sample sizes, other than the true sample sizes included as 
  an attribute in goObject, qm and pqn0. When qm = NULL, the genes profiled in goObject are compared with a 
  subsample of them, those profiled in pqn0 (is there equivalence between a set of genes and a restricted subset,
  e.g. those overexpressed under a disease, in terms of their profiles?).
  When pqn0 = NULL, an equivalence test between two profiles with no genes in common is performed.
  
  In the 'GOProfileHtest' interface, the one-sided confidence interval for the squared Euclidean distance
  is computed from the distance and its standard error stored in the corresponding fields of the argument
  goObject, itself typically an object of class 'GOProfileHtest' resulting from a call to 'compareGOProfiles' with 
  simplify=T.
 
  In the default interface, the 'goObject' argument is previously converted into an object of class 'ExpandedGOProfile'
  and then this interface is used.
 
  If the argument 'd0' is not provided it is computed as \eqn{d0 <- s * equivEpsilon^2},
  where 's' stands for the number of non empty GO nodes in any of the GO profiles
  being compared.
  
}
\value{
In the 'ExpandedGOProfile' interface, the result is an object of class "list" containg one or more "htest" objects, each of which may come from previous profiles comparisons.
In the other interfaces, the result is a single "htest" object.
Each one of these "htest" objects has the following fields:
\item{statistic}{test statistic, (distance - d0) / se}
\item{parameter}{d0 and the sample sizes (number of genes) n and m}
\item{p.value}{associated p-value to test the null hypothesis of profiles inequivalence}
\item{conf.int}{asymptotic one-sided confidence interval for the squared euclidean distance. Its attribute "conf.level" contains its nominal confidence level.}
\item{estimate}{squared euclidean distance between the contracted profiles. Its attribute "se"
contains its standard error estimate}
\item{data.name}{a character string giving the names of the data}
\item{alternative}{a character string describing the alternative hypothesis (always 'Equivalence or similarity, true squared Euclidean distance between the contracted profiles is less than d0'}
}
\author{Jordi Ocana}

\seealso{'compareGOProfiles'}
\examples{

data(prostateIds)

expandedWelsh <- expandedProfile(welsh01EntrezIDs[1:100], onto="ANY",
                                 level=2, orgPackage="org.Hs.eg.db")
expandedSingh <- expandedProfile(singh01EntrezIDs[1:100], onto="ANY",
                                 level=2, orgPackage="org.Hs.eg.db")
commonGenes <- intersect(welsh01EntrezIDs[1:100], singh01EntrezIDs[1:100])
commonExpanded <- expandedProfile(commonGenes, onto="ANY", level=2, 
                                  orgPackage="org.Hs.eg.db")

### FUnciona si fem: 

equivMF <-equivalentGOProfiles (expandedWelsh[["MF"]], 
                              qm  = expandedSingh[["MF"]], 
                              pqn0= commonExpanded[["MF"]])

equivsList <- lapply(1:length(expandedWelsh), 
                     function (onto){
                          equivalentGOProfiles (expandedWelsh[[onto]], 
                                qm  = expandedSingh[[onto]], 
                                pqn0= commonExpanded[[onto]])
                     }
)

}
\keyword{htest}
