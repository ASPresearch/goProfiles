\name{omimIds}
\alias{omimIds}
\alias{diseaseIds}
\alias{morbidmapIds}
\alias{dominantIds}
\alias{recessiveIds}
\alias{dominantIdsEBI}
\alias{recessiveIdsEBI}
\alias{dominantIdsNCBI}
\alias{recessiveIdsNCBI}
\docType{data}                                                  
\title{Entrez identifiers for disease-related genes in the OMIM database}
\description{Entrez identifiers for several lists of genes related with human disease.
\item \code{diseaseIds} contains the Entrez identifiers corresponding to disease-related genes found in the 
OMIM database. This list has been manually curated by Nuria Lopez-Bigas et al. who kindly provided it to us. 
\item \code{morbidmapIds} contains the Entrez identifiers for all the genes in the morbidmap table. 
This list would correspond to disease-related genes if there had been no manual curation, as in the previous list ('diseaseIds')
\item \code{dominantIds} contains the Entrez identifiers for dominant genes after manual curation by Nuria Lopez-Bigas who has kindly allowed us to include them in the package.
\item \code{recessiveIds} contains the Entrez identifiers for recessive genes after manual curation by Nuria Lopez-Bigas who has kindly allowed us to include them in the package.
\item \code{dominantIdsEBI} contains the Entrez identifiers for dominant genes in the EBI version of the OMIM database recovered using SRS with the term 'dominant' in the KEYWORDS field.
\item \code{recessiveIdsEBI} contains the Entrez identifiers for recessive genes in the EBI version of the OMIM database recovered using SRS with the term 'recessive' in the KEYWORDS field.
\item \code{dominantIdsNCBI} contains the Entrez identifiers for dominant genes in the NCBI version of the OMIM database recovered using ENTREZ with the term 'dominant' in the CLINICAL field.
\item \code{recessiveIdsNCBI} contains the Entrez identifiers for recessive genes in the NCBI version of the OMIM database recovered using ENTREZ with the term 'recessive' in the CLINICAL field.
}
\details{Lopez-Bigas et al. analyzed the distribution of functional categories in genes causing disease in human.
They did several comparisosn which can also be done using \code{goProfiles}. In order to perform these comparisons
we first tried to obtain the same lists of genes using standard database browsers, such as 'SRS', at the European Bioinformatics Institute, or 'Entrez', at the National Center for Biotechnological Information. Curiously both approaches provided very different lists so we asked the authors for their data and they kindly provided them to us. In order to facilitat the use of functions included in \code{goProfiles} we have trimmed the list of recessive and dominant genes so that (i) They become exclussive (no gene belows to both lists) (2) They are both included in the \code{diseaseIds} list. This eliminated 39 genes (out of 639) from the list of recessive genes and  52 genes (out of 414) from the list of dominant genes}
\references{
Lopez-Bigas, N., Blencowe, B.J. and Ouzounis, C.A.,
Highly consistent patterns for inherited human diseases at the molecular level,
Bioinformatics, 2006, 22 (3), 269-277.
}
\usage{data(omimIds)}
\format{Each dataset is a character vector with a different number of elements which (should) correspond to valid Entrez identifiers}
\examples{
data(omimIds)
}
\keyword{datasets}
