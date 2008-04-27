\name{basicProfile}
\alias{basicProfile}
\title{Builds basic functional profile}
\description{
Compute basic functional profile for a given list of genes/GO identifiers, a given ontology at a given level of the GO}
\usage{
basicProfile(genelist, idType = "Entrez", onto = "ANY", level = 2,orgPackage=NULL, anotPackage=NULL, 
	ord = TRUE, multilevels = NULL, empty.cats = TRUE, cat.names = TRUE, na.rm = FALSE)
}
\arguments{
  \item{genelist}{List of genes on which the Profile has to be based}
  \item{idType}{Type of identifiers for the genes. May be 'Entrez' (default), BiocProbes or GoTermsFrame (see details below).}
  \item{onto}{Ontology on which the profile has to be built}
  \item{level}{Level of the ontology at which the profile has to be
    built}
  \item{orgPackage}{Name of a Bioconductor's organism annotations package
    ('org.Xx-eg-db'). This field must be provided if the gene list passed
    to the function is either a character vector of 'Entrez' (NCBI) identifiers or a character vector of probe names}
  \item{anotPackage}{Name of Bioconductor's microarray annotations package. This field must be provided if the gene list passed to the function is a character vector of probe names}
  \item{ord}{Set to 'TRUE' if the profile has to appear ordered by the category names}
  \item{multilevels}{If it is not NULL it must be a vector of GO categories that defines the level at where the profile is built}
  \item{empty.cats}{Set to 'TRUE' if empty categories should appear in the profile}
  \item{cat.names}{Set to 'TRUE' if the profile has to contain the names of categories}
  \item{na.rm}{Set to 'TRUE' if NAs should be removed}
}
\value{
An object of class GOProfile (one or more data frames in a list named by the ontologies)
}
\details{
The function admits three types of entries: Entrez ('Entrez'), Bioconductor probe set names ('BioCprobes') 
or a special type of data frames ('GOTermsFrames'). If the identifier type are 'BioCprobes' 
then an annotation package name must be provided too.
}

\references{Sanchez-Pla, A., Salicru, M. and Ocana, J.
Statistical methods for the analysis of high-throughput data based on
functional profiles derived from the gene ontology.
Journal of Statistical Planning and Inference,
Volume 137, Issue 12, Pages 3975-3989, 2007
}
\author{Alex Sanchez}
\seealso{expandedProfile}
\examples{
data(CD4Ids)
CD4.MF.Profiles <-basicProfile(genelist=CD4LLids, onto='MF', level=2, orgPackage="org.Hs.eg.db")
print(CD4.MF.Profiles)
}
\keyword{htest}
