\name{plotProfiles}
\alias{plotProfiles}
\title{Plot functional profiles}
\description{
Plots basic functional profiles created with the 'basicProfile' instruction.
If several profiles have to be plot together they must be first merged using the 'mergeProfiles' function.
The labels of the Y-axis of the plots are the descriptions of the GO Terms. If the label is longer than 20 characters
it is truncated and ended by three dots.}
\usage{
plotProfiles(aProf, aTitle = "Functional Profile", anOnto = NULL, percentage = FALSE, 
HORIZVERT = TRUE, legendText = NULL, colores = c("white", "red"), multiplePlots = F, multipleWindows = T, labelWidth=25,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{aProf}{Functional profile to plot}
  \item{aTitle}{Title for the figures}
  \item{anOnto}{Ontology (to appear in the title)}
  \item{percentage}{Plot absolute or relative frequencies (not summing to 100)}
  \item{HORIZVERT}{Plot horizontal or vertical bars}
  \item{legendText}{Text of the legend for the plot}
  \item{colores}{Colors to be used}
  \item{multiplePlots}{Plot all profiles for a given dataset in one figure}
  \item{multipleWindows}{Open a new window after each plot}
  \item{labelWidth}{Width of Y axis labels (Names of GO categories) in the plot}
  \item{...}{Other graphical parameters that should be passed for plotting}
}
\value{
The plot
}
\examples{
require(goProfiles)
data(prostateIds)
welsh.MF <- basicProfile (welsh01EntrezIDs[1:100], onto="MF", level=2, orgPackage="org.Hs.eg.db") 
singh.MF <- basicProfile (singh01EntrezIDs[1:100], onto="MF", level=2, orgPackage="org.Hs.eg.db") 
plotProfiles(welsh.MF,'Functional profiles for Welsh dataset',percentage=TRUE)
welsh.singh.MF <-mergeProfilesLists(welsh.MF, singh.MF, profNames=c("Welsh", "Singh"))
plotProfiles(welsh.singh.MF , percentage=TRUE, multiplePlots=TRUE, labelWidth=30)
}
\author{Alex Sanchez}
\keyword{hplot}
