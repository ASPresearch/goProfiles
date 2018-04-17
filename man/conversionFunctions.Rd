\name{conversionFunctions}
\alias{as.GOTerms.frame}
\alias{as.GOTerms.list}
\alias{BioCpack2EntrezIDS}
\alias{BioCpack2Profiles}
\alias{BioCprobes2Entrez}
\alias{GOTermsFrame2GOTermsList}
\title{Functions to transformconvert objects between different types}
\description{
These functions transform data from one classtype into another, or pack
simple processes such as compute the profiles needed for one annotations
package.
}
\usage{
as.GOTerms.frame(myGOTermsList, na.rm=TRUE)
as.GOTerms.list(genelist, probeType, orgPackage=NULL, anotPkg=NULL, 
onto="any", na.rm=TRUE)
BioCpack2EntrezIDS(anotPkg, na.rm=TRUE)
BioCpack2Profiles(anotPkg, orgPackage, level=2, na.rm=TRUE, expanded=FALSE)
BioCprobes2Entrez(probeslist , anotPkg, na.rm=TRUE)
GOTermsFrame2GOTermsList(myGOTermsFrame, evid=FALSE)
}
\arguments{
  \item{myGOTermsList}{GOTermsList object to transform}
  \item{myGOTermsFrame}{GOTermsFrame object to transform}
  \item{genelist}{List of genes (Entrez Ids) to transform}
  \item{evid}{Type of evidence supporting the selected GO Terms}
  \item{na.rm}{Flag indicating if those ids returning NA must be removed from the output}
  \item{probeType}{Type of probes to transform into Entrez Ids}
  \item{probeslist}{List of probes to transform into Entrez Ids}
  \item{orgPackage}{Name of the organism ('org.Xx.eg.db') annotation package}
  \item{anotPkg}{Name of the chip annotation package}
  \item{level}{GO level at which the profile is built}
  \item{onto}{ontology}
  \item{expanded}{Flag to decide if an expanded profile has to be computed} 
}
\value{
Every function returns a transformed object or a list of computed profiles
}
\details{
Not yet available
}
\author{Alex Sanchez}
\examples{
data(CD4Ids)
myGOTermsList <- GOTermsList(CD4LLids[1:5], orgPkg="org.Hs.eg.db")
myGOTermsFrame<- as.GOTerms.frame(myGOTermsList, na.rm=TRUE)
GOTermsFrame2GOTermsList(myGOTermsFrame, evid=FALSE)
}
\keyword{manip}
