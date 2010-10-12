\name{expandedProfile}
\alias{expandedProfile}
\title{Builds expanded profiles}
\description{
Expanded profiles are used mainly for comparison of profiles based on the theory developed by Sanchez et al (2007) (see references)}
\usage{
expandedProfile(genelist, idType = "Entrez", onto = "ANY", level = 2,orgPackage=NULL, anotPackage=NULL, 
	multilevels = NULL, ord = TRUE, na.rm = TRUE, percentage = TRUE)
}
\arguments{
  \item{genelist}{List of genes on which the Profile has to be based}
  \item{idType}{Type of identifiers for the genes. Use 'Entrez' preferably}
  \item{onto}{Ontology on which the profile has to be built}
  \item{level}{Level of the ontology at which the profile has to be
    built}
   \item{orgPackage}{Name of a Bioconductor's organism annotations package
    ('org.Xx-eg-db'). This field must be provided if the gene list passed
    to the function is either a character vector of 'Entrez' (NCBI) identifiers or a character vector of probe names}
  \item{anotPackage}{Name of Bioconductor annotations package. This field must be provided if the gene list passed to the function is a character vector of probe names}
  \item{ord}{Set to 'TRUE' if the profile has to appear ordered by the category names}
  \item{multilevels}{If it is not NULL it must be a vector of GO categories that defines the level at where the profile is built}
  \item{na.rm}{Set to 'TRUE' if NAs should be removed}
  \item{percentage}{Set to 'TRUE' if the profile must be built using percentages}
}
\value{
An object of class GOProfile containing an expanded profile}
\details{
The function admits three types of entries: Entrez ('Entrez'), Bioconductor probe set names ('BioCprobes') 
or a special type of data frames ('GOTermsFrames'). If the identifier type are 'BioCprobes' 
then an annotation package name must be provided too.
}
\references{Sanchez-Pla, A., Salicru, M. and Ocana, J.
Statistical methods for the analysis of high-throughput data based on
functional profiles derived from the gene ontology.
Journal of Statistical Planning and Inference,
Volume 137, Issue 12, Pages 3975-3989, 2007.
}
\author{Alex Sanchez}
\seealso{basicProfile}
\examples{
data(CD4Ids)
CD4.Expanded <-expandedProfile(genelist=CD4LLids[1:50], onto='MF', level=2, orgPackage="org.Hs.eg.db")
}
\keyword{htest}
