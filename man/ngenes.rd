\name{ngenes}
\alias{ngenes}
\alias{ngenes.default}
\alias{ngenes.numeric}
\alias{ngenes.matrix}
\alias{ngenes.ExpandedGOProfile}
\alias{ngenes.BasicGOProfile}
\title{Returns the number of genes that lead to this GO profile (an object of
class ExpandedGOProfile, BasicGOProfile or assimilable to them)}
\description{The information contained in one or more lists of genes may be
summarized by their GO profiles, that is to say, the absolute or relative 
frequencies of annotations or hits in all 
the classes or nodes of a given leven in a given GO ontology, or by the 
corresponding frequencies in a selected set of nodes (possibly belonging to more
than one GO level but not hierarchicaly related). This function returns the
number of genes in each list that were annotated to compute the profiles
}
\usage{
ngenes(pn, i=NULL)
\method{ngenes}{default}(pn, i=NULL)
\method{ngenes}{numeric}(pn, i=NULL)
\method{ngenes}{matrix}(pn, i=NULL)
\method{ngenes}{ExpandedGOProfile}(pn, i=NULL)
\method{ngenes}{BasicGOProfile}(pn, i=NULL)
}
\arguments{
\item{pn}{an object of class ExpandedGOProfile or BasicGOProfile representing 
one or more "sample" expanded GO profiles for a fixed ontology, or a numeric 
vector interpretable as a GO profile (expanded or not), or a
frequency matrix (see the 'Details' section)}
\item{i}{i-th profile in the case of more than one profiles. A vector with
the number of genes of all profiles is returned if this argument is absent}
}

\details{
 Given a list of n genes, and a set of s GO nodes X, Y, Z, ... in a given ontology
 (BP, MF or CC), its  associated (contracted) "basic profile" is the frequencies vector (either
 absolute or relative frequencies) of annotations or hits of the n genes in each node.
 For a given node, say X, this frequency includes all annotations for X alone, for X and Y,
 for X and Z and so on. Thus, as relative frequencies, its sum is not necessarily one,
 or as absolute frequencies their sum is not necessarily n.
 On the other hand, an "expanded profile" corresponds to the frequencies in ALL OBSERVED NODE
 COMBINATIONS. That is, if n genes have been profiled, the expanded profile stands
 for the frequency of all hits EXCLUSIVELY in nodes X, Y, Z, ..., jointly with
 all hits simultaneously in nodes X and Y (and only in X and Y), simultaneously in X and Z,
 in Y and Z, ... , in X and Y and Z (and only in X,Y,Z), and so on. Thus, their sum is one.

 An object of S3 class 'ExpandedGOProfile' is, essentially, a 'data.frame' object
 with each column representing an expanded profile.
 The row.names attribute codifies the node combinations and each
 data.frame column (say, each profile) has an attribute, 'ngenes', indicating the
 number of profiled genes.
}

\value{A vector with the number of genes annotated in one or more GO profiles}
\seealso{BasicGOProfile object, ExpandedGOProfile object}
\author{Jordi Ocana}

\examples{
require("org.Hs.eg.db")
data(prostateIds) 
# To improve speed, use only the first 100 genes:
list1 <- welsh01EntrezIDs[1:100]
prof1 <- expandedProfile(list1, onto="MF", level=2, orgPackage="org.Hs.eg.db", na.rm=TRUE)$MF
length(list1)
# Only a subset of the initial gene list are annotated in the profile
ngenes(prof1)
}
