\name{prostateIds}
\alias{prostateIds}
\alias{welsh01EntrezIDs}
\alias{welsh05EntrezIDs}
\alias{singh01EntrezIDs}
\alias{singh05EntrezIDs}
\docType{data}
\title{Prostate cancer-related genes}
\description{
Entrez identifiers for genes related with Prostate Cancer selected from two datasets analyzed by Welsh et al. (2001) and Singh et al. (2002) respectively.
The genes have been selected from freely available datasets in the internet using a standard workflow for selecting differentially expressed genes.
The dataset contains 4 character vectors, each corresponding to the entrez identifiers of the genes selected at a 5\% and 1\% significance level 
from the Welsh and Singh dataset respectively.
\item \code{welsh05EntrezIDs} List of genes selected from Welsh et al. study at a 0.05 significance level.
\item \code{welsh01EntrezIDs} List of genes selected from Welsh et al. study at a 0.01 significance level.
\item \code{singh05EntrezIDs} List of genes selected from Singh et al. study at a 0.05 significance level.
\item \code{singh01EntrezIDs} List of genes selected from Singh et al. study at a 0.01 significance level.
}
\usage{data(prostateIds)}
\source{
\item John B. Welsh, Lisa M. Sapinoso, Andrew I. Su, Suzanne G. Kern, Jessica Wang-Rodriguez,
Christopher A. Moskaluk, Jr. Frierson, Henry F., and Garret M. Hampton. Analysis of Gene
Expression Identifies Candidate Markers and Pharmacological Targets in Prostate Cancer.
Cancer Res, 61(16):5974-5978, 2001.
\item Singh, Dinesh and Febbo, Phillip G and Ross, Kenneth and Jackson, Donald
G and Manola, Judith and Ladd, Christine and Tamayo, Pablo and Renshaw,
Andrew A and D'Amico, Anthony V and Richie, Jerome P and Lander, Eric S
and Loda, Massimo and Kantoff, Philip W and Golub, Todd R and Sellers,
William R.Gene expression correlates of clinical prostate cancer behavior.
Cancer Cell, 2002, Mar., 1(2) 203-209, 2002.
}
\format{Each dataset is a character vector with a different number of elements which (should) correspond to valid Entrez identifiers}
\examples{
data(prostateIds)
}
\keyword{datasets}
