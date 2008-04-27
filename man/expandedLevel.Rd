\name{expandedLevel}
\alias{expandedLevel}
\alias{expandTerm}
\title{Function to create expanded levels which can contain GO Terms at
  different GO levels}
\description{
This function, combined with function \code{expandTerm}, allows to
create mixed levels which can contain terms belonging to different GO
levels.
Specifically one can take one (or several, but one by one) term at a
given GO level and expand it into its children terms using function
expandTerm and then combine them into a new level using this function.
}
\usage{
expandedLevel(LevelTerms, Term2Expand, onto)
expandTerm(GOTerm, onto)
}
\arguments{
  \item{LevelTerms}{Other terms which have not been expanded, and will
    be combined with the expanded ones}
  \item{Term2Expand}{The GO term which will be substituted by its
    children terms}
  \item{GOTerm}{The GO term which will be substituted by its
    children terms}
  \item{onto}{The ontology ('MF','BP','CC'}
}
\value{
The value returned is the vector combining the original terms with the
children of the term that had to be expanded.
}
\author{Alex Sanchez}
\examples{
got<-toTable(GOTERM)[,2:3]
desc<-function(s) got[got[,1]==s,2]
MFLevel2<-getGOLevel("MF",2)
bindingLevel2<-MFLevel2 [2]
bindingLevel3 <- expandTerm(bindingLevel2,"MF")
print(descbindingLevel3<-as.matrix(sapply(bindingLevel3,desc )))
mixedLevel<-c(MFLevel2[-2],bindingLevel3)
print(mixedLevel<-as.matrix(sapply(mixedLevel,desc )))
}
\keyword{manip}
