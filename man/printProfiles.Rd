\name{printProfiles}
\alias{printProfiles}
\title{Print functional profiles}
\description{
Prints basic functional profiles created with the 'basicProfile' instruction.
Allows for several formatting operations such as truncating long labels, 
removing empty categories or choosing between absolute or relative frequencies.
If several profiles have to be printed together they must be first merged using the 'mergeProfiles' function.
}
\usage{
printProfiles(aProf, aTitle = "Functional Profile", anOnto = NULL, percentage = FALSE, 
      Width=25, emptyCats=FALSE)
}
\arguments{
  \item{aProf}{Functional profile to plot}
  \item{aTitle}{Title for the figures}
  \item{anOnto}{Ontology (to appear in the title)}
  \item{percentage}{Plot absolute or relative frequencies (not summing to 100)}
  \item{Width}{Maximum width for the description of GO categories}
  \item{emptyCats}{Set to 'TRUE' if empty categories should appear in the profile}
}
\value{
The printout
}
\examples{
require(goProfiles)
data(prostateIds)
welsh.MF <- basicProfile (welsh01EntrezIDs[1:100], onto="MF", level=2, orgPackage="org.Hs.eg.db") 
singh.MF <- basicProfile (singh01EntrezIDs[1:100], onto="MF", level=2, orgPackage="org.Hs.eg.db") 
printProfiles(welsh.MF,'Functional profiles for Welsh dataset',percentage=TRUE, anOnto='MF')
welsh.singh.MF <-mergeProfilesLists(welsh.MF, singh.MF, profNames=c("Welsh", "Singh"))
printProfiles(welsh.singh.MF, percentage=TRUE, emptyCats=TRUE)
}
\author{Alex Sanchez}
\keyword{hplot}
