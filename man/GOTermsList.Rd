\name{GOTermsList}
\alias{GOTermsList}
\alias{getAncestorsLst}
\alias{getGOLevel}
\title{Functions to create and manage lists of GO terms associated with a vector of 'Entrez' identifiers}
\description{
These functions prepare data to be processed by the 'basicProfile' function.
To create a profile a set of GOterms belonging to one or more ontologies is needed
The terms belonging to each gene must be given separately so that they can be counted.
This function queries the environment 'GOENTREZID2GO' with the vector of Entrez terms
and formats the output into a list whose components -one per Entrez term- contain the most
specific GO identifiers associated with this term.
}
\usage{
GOTermsList(LLids, onto = "any", evid = "any", na.rm = FALSE, orgPkg )
getAncestorsLst(GOtermslist, onto, unique.ancestor=TRUE, na.rm=TRUE, combine=TRUE)
getGOLevel(onto, level)
}
\arguments{
  \item{LLids}{Character vector of \code{Entrez} (formerly \code{Locuslink} identifiers)}
  \item{onto}{ontology to be queried using the genes list}
  \item{evid}{type of evidence supporting the selected GO Terms}
  \item{na.rm}{flag indicating if those ids returning NA must be removed from the output}
  \item{orgPkg}{Organism annotation package ('org.Xx.eg.db') required to obtain the GO terms associated with the Entrez identifiers}
  \item{GOtermslist}{List produced by a call to function \code{GOTermsList}}
  \item{unique.ancestor}{Flag to remove repeated ancestor identifiers}
  \item{combine}{Flag to combine ancestors}
  \item{level}{GO level at which the profile is built}
}
\value{
A list whose components -one per Entrez term- are character vectors with the most
specific GO identifiers associated with this term
}
\details{
During the call to this function there may appear 
two types of NAs.

By one side if a name is not mapped in LocusLink this yields an NA that must
be eliminated because nothing can be found through LL about this name

By another side if a gene is identified in LL but yields NA
it seems to mean that it is not mapped in the GO

This may be eliminated but it may be worth the pity to keep track of them
and to put these terms in an 'Seemingly unnanotated' category.
In the case that its number was very high it migt suggest reviewing the list
or reconsidering the results.
}
\author{Alex Sanchez}
\seealso{\code{getAncestorsLst}}
\examples{
#data(CD4Ids)
#simpleLLids<- as.character(c(2189,5575,5569,11)) #1 is not a Locuslink identifier
#simpleGOlist<- GOTermsList (simpleLLids, orgPkg="org.Hs.eg.db")
#print(simpleGOlist.CC<-GOTermsList (simpleLLids,"CC", orgPkg="org.Hs.eg.db"))
#print(simpleGOlist.IEA<-GOTermsList (simpleLLids,evid="IEA",na.rm=TRUE, orgPkg="org.Hs.eg.db"))
}
\keyword{manip}
