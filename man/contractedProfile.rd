\name{contractedProfile}
\alias{contractedProfile}
\alias{contractedProfile.ExpandedGOProfile}
\alias{contractedProfile.default}
\title{Converts an expanded GO profile into a basic (contracted) GO profile}
\description{
  Converts an object of class 'ExpandedGOProfile', or assimilable to it,
  in an object of class 'BasicGOProfile'
}

\usage{
contractedProfile(prof, nams = NULL)
\method{contractedProfile}{ExpandedGOProfile}(prof, nams = NULL)
\method{contractedProfile}{default}(prof, nams = NULL)
}

\arguments{
  \item{prof}{ an expanded GO profile, i.e. and object of class 'ExpandedGOProfile', or
  a numeric vector assimilable to an expanded profile, see the "details" section }
  \item{nams}{ optionally, the names of the annotated combinations of GO nodes
  whose frequency is represented in the expanded profile, see the "details" section }
}
\details{
 Given a list of n genes, and a set of s GO nodes X, Y, Z, ... in a given ontology
 (BP, MF or CC), its  associated (contracted) "profile" is the frequencies vector (either
 absolute or relative frequencies) of annotations or hits of the n genes in each node.
 For a given node, say X, this frequency includes all annotations for X alone, for X and Y,
 for X and Z and so on. Thus, as relative frequencies, its sum is not necessarily one,
 or as absolute frequencies their sum is not necessarily n. Basic contracted profiles
 are represented by objects of S3 class 'BasicGOProfile'.
 On the other hand, an "expanded profile" corresponds to the frequencies in ALL OBSERVED NODE
 COMBINATIONS. That is, if n genes have been profiled, the expanded profile stands
 for the frequency of all hits EXCLUSIVELY in nodes X, Y, Z, ..., jointly with
 all hits simultaneously in nodes X and Y (and only in X and Y), simultaneously in X and Z,
 in Y and Z, ... , in X and Y and Z (and only in X,Y,Z), and so on. Thus, their sum is one.
 Expanded profiles are represented by objects of S3 class 'ExpandedGOProfile'.
 The generic function 'contractedProfile' "contracts" an expanded profile, either
 represented by a 'ExpandedGOProfile' object or a numeric vector interpretable as
 an expanded profile, in order to obtain its contracted profile representation.
 
 The \code{rownames} attribute of an 'ExpandedGOProfile' or, equivalently, the
 \code{names} attribute of a vector representing an expanded profile, or the
 \code{nams} argument, must represent the GO nodes combinations separating the
 node names with dots, ".", for example: "X", "Y", "Z", "X.Y", "X.Z", "Y.Z",
 "X.Y.Z" and so on.
}
\value{
An object of class 'BasicGOProfile' the contracted profile representation of the
expanded profile}
\author{Jordi Ocana}

\examples{
data(prostateIds)
expandedWelsh <- expandedProfile(welsh01EntrezIDs[1:100], onto="MF",
                        level=2, orgPackage="org.Hs.eg.db")
reContractedWelsh <- contractedProfile(expandedWelsh[["MF"]])
print(expandedWelsh)
print(reContractedWelsh)
class(reContractedWelsh)
ngenes(reContractedWelsh)
}


