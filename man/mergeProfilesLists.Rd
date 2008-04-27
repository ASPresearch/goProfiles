\name{mergeProfilesLists}
\alias{mergeProfilesLists}
\title{Combines two lists of profiles into one}
\description{
Combines two lists of profiles, that is two lists with three components, 
'MF', 'BP', 'CC' into a single one.
}
\usage{
mergeProfilesLists(profilesList1, profilesList2, emptyCats = F, profNames = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{profilesList1}{First list to combine}
  \item{profilesList2}{Second list to combine}
  \item{emptyCats}{Boolean. Set to TRUE if there are empty categories that should be accounted for
  in any of the profiles}
  \item{profNames}{Names for the profiles (optional). If missing they are set to
  'Frequency-1', 'Frequency-2',etc.}
}
\value{
A list of profiles with more than one column each.
}
\author{Alex Sanchez}
\examples{
data(sampleProfiles)
merged<-mergeProfilesLists (basicSingh01, basicWelsh01, profNames=c('Welsh','Singh'))
}
\keyword{utilities}
