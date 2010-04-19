\name{fisherGOProfiles}
\alias{fisherGOProfiles}
\alias{fisherGOProfiles.numeric}
\alias{fisherGOProfiles.matrix}
\alias{fisherGOProfiles.BasicGOProfile}
\alias{fisherGOProfiles.ExpandedGOProfile}
\title{GO Class-by-class Fisher tests in lists of genes characterized by their functional profiles}
\description{
Given two lists of genes, both characterized by their frequencies of annotations
(or "hits") in the same set of GO nodes (also designated as GO terms or GO classes),
for each node determine if the annotation frequencies depart from what is expected
by chance. The annotation frequencies are specified in the "GO profiles" arguments
\code{pn}, \code{qm} and \code{pn}.
Both samples may share a common subsample of genes, with GO profile
\code{pqn0}. The analysis is based on the Fisher's exact test, as is
implemented by \code{fisher.test} R function, followed by p-value adjustment for
multitesting based on function \code{p.adjust}. Usually, this function will be
called after a significant result on \code{compareGOProfiles} which performs
global (all GO nodes simultaneously) profile comparisons (with better
type I and type II error control), to identify the more rellevant nodes.}
\usage{
fisherGOProfiles(pn, \ldots)
\method{fisherGOProfiles}{numeric}(pn, qm=NULL, pqn0=NULL,
    n = ngenes(pn), m = ngenes(qm), n0 = ngenes(pqn0),
    method = "BH", simplify=T, expanded=F, \ldots)
\method{fisherGOProfiles}{matrix}(pn, n, m, method = "BH", \ldots)
\method{fisherGOProfiles}{BasicGOProfile}(pn, qm=NULL, pqn0=NULL,
    method = "BH", goIds=T, \ldots)
\method{fisherGOProfiles}{ExpandedGOProfile}(pn, qm=NULL, pqn0=NULL,
    method = "BH", simplify=T, \ldots)
}
\arguments{
\item{pn}{an object of class \code{BasicGOProfile} or \code{ExpandedGOProfile}
representing a "sample" GO profile for a fixed ontology, or a numeric vector
interpretable as a GO profile (expanded or not), or a two-dimensional
frequency matrix (see the 'Details' section). This is a required argument}
\item{qm}{similarly, an object representing a "sample" GO profiles for a fixed ontology}
\item{pqn0}{an object representing a "sample" GO profile for a fixed ontology}
\item{n}{the number of genes profiled in pn}
\item{m}{the number of genes profiled in qm}
\item{n0}{the number of genes profiled in pqn0}
\item{method}{the p-values adjusting method for multiple comparisons; the same
possibilities as in standard R function \code{p.adjust}}
\item{expanded}{boolean; are these numeric vectors representing expanded profiles?}
\item{simplify}{should the result be simplified, if possible? See the 'Details' section}
\item{goIds}{if TRUE, each node is represented by its GO identifier}
\item{...}{other arguments (to be passed to \code{p.adjust} or \code{fisher.test} functions)}
}

\details{
 Given a list of \code{n} genes, and a set of \code{s} GO classes or nodes
 X, Y, Z, ... in a given ontology
 (BP, MF or CC), its  associated ("contracted" or "basic") "profile" is the
 absolute frequencies vector of annotations or hits of the \code{n} genes in each
 one of the \code{s} GO nodes.
 For a given node, say X, this frequency includes all annotations for X alone, for X and Y,
 for X and Z and so on. Thus, as relative frequencies, its sum is not necessarily one,
 or as absolute frequencies their sum is not necessarily \code{n}.
 On the other hand, an "expanded profile" corresponds to the relative frequencies
 in ALL NODE COMBINATIONS. That is, if \code{n} genes have been profiled, the
 expanded profile stands
 for the frequency of all hits EXCLUSIVELY in node X, exclusively in node Y,
 exclusively in Z, ..., jointly with
 all hits simultaneously in nodes X and Y (and only in X and Y), simultaneously in X and Z,
 in Y and Z, ... , in X and Y and Z (and only in X,Y,Z), and so on.
 Thus, their sum is one.
 
 Let \code{n}, \code{m} and \code{n0} designate the total number of genes
 profiled in \code{pn}, \code{qm} and \code{pqn0} respectively.
 According to these profiles, n[i], m[i] and n0[i] genes are annotated
 for node 'i', i = 1, \ldots, \code{s}. Note that the sum of all the n[i] not
 necessarily equals \code{n} and so on.
 If not NULL, \code{pqn0} stands for the profile of the \code{n0}
 genes common to the gene lists that gave rise to \code{pn} and \code{qm}.
 \code{fisherGOProfiles} builds a \code{s}x2 absolute frequencies matrix
 \tabular{rcc}{
  GO node 1             \tab N[1,1] \tab N[1,2] \cr
  GO node 2             \tab N[2,1] \tab N[2,2] \cr
  \ldots                \tab \ldots \tab \ldots \cr
  GO node \code{s}      \tab N[2,1] \tab N[s,2] \cr
 }
 with column totals N1 and N2 (not necessarily equal to the column sums)
 and performs a Fisher's exact test over each one of the 2x2 tables
 \tabular{rcc}{
  GO node i             \tab N[i,1]     \tab N[i,2]     \cr
  All nodes except i    \tab N1 - N[i,1]\tab N2 - N[i,2]\cr
 }
 followed by a p-value correction for multiplicity in testing.
 If \code{pqn0} is NULL, then both gene lists do not have any genes in common,
 N[i,1] = n[i] and N[i,2] = m[i], and N1 = n, N2 = m, n0 = 0.
 Otherwhise (if \code{pqn0} is not NULL) N[i,1] = n[i] - n0[i], N1 = n - n0 and
 N[i,2] = n0[i], N2 = n0 if \code{qm} is NULL, or N[i,2] = m[i], N2 = m if \code{qm}
 is not NULL.
 
 In other words, this function provides a general setting for diverse, common
 in practice, situations where a node-by-node analysis is required.
 When \code{pqn0} = NULL, two lists with no genes in common are compared.
 Otherwise, when \code{qm} = NULL, the genes profiled in \code{pn} are compared
 with a subsample of them, those profiled in \code{pqn0} (a set of genes vs a restricted subset,
 e.g. those overexpressed under a disease). Finally, if both arguments \code{qm}
 and \code{pqn0} are not NULL (\code{pn} is always required) two gene lists with
 some genes in common are analised.

 If both \code{qm} and \code{pqn0} are NULL, \code{pn} should correspond to an
 absolute frequencies matrix with \code{s} rows and 2 columns.

 The arguments \code{n}, \code{m} or \code{n0} are only required in case of
 numeric vectors or matrices specifying profiles but lacking  the 'ngenes' attribute.
}

\value{
A list containing max(ncol(pn),ncol(qm),ncol(pqn0)) p-values numeric vectors,
or a single p-values vector if max(ncol(pn),ncol(qm),ncol(pqn0))==1 and simplify == T.
}
\seealso{fitGOProfile, compareGOProfiles, equivalentGOProfiles}
\references{Sanchez-Pla, A., Salicru M. and Ocana, J. Statistical methods for the analysis of highthroughput data based on functional profiles derived from the gene ontology. Journal of
Statistical Planning and Inference, 2007.}
\author{Jordi Ocana}

\examples{

}
