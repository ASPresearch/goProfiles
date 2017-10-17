
compareGOProfiles <- function(pn, qm=NULL, pqn0=NULL, n = ngenes(pn), m = ngenes(qm), n0 = ngenes(pqn0),
    method = "lcombChisq", ab.approx = "asymptotic", confidence = 0.95,
#    d0 = NULL, equivEpsilon = 0.05,   # <----------------------- suprimir de l'ajuda
    nsims = 10000,
    simplify=T, ...){
#
# Compare two samples of genes in terms of their GO profiles pn and qm. Both
# samples may share a common subsample of genes, with GO profile pqn0.
# Are they just random samples taken from the same "population" GO profile?
#
# Description:
# 'compareGOProfiles' implements some inferential procedures to solve the preceding question. These procedures are based on
# asymptotic properties of the squared euclidean distance between the contracted versions of pn and qm
#
# Usage:
#     compareGOProfiles(pn, qm=NULL, pqn0=NULL, n = ngenes(pn), m = ngenes(qm), n0 = ngenes(pqn0),
#           method = "lcombChisq", ab.approx = "asymptotic", confidence = 0.95,
#           d0 = NULL, equivEpsilon = 0.05,
#           nsims = 10000,
#           simplify=T, ...)
#
# Arguments:
# pn:           an object of class ExpandedGOProfile representing one or more
#               "sample" expanded GO profiles for a fixed ontology (see the 'Details' section)
# qm:           an object of class ExpandedGOProfile representing one or more
#               "sample" expanded GO profiles for a fixed ontology (see the 'Details' section)
# pqn0:         an object of class ExpandedGOProfile representing one or more
#               "sample" expanded GO profiles for a fixed ontology (see the 'Details' section)
# n:            a numeric vector with the number of genes profiled in each column of pn.
#               This parameter is included to allow the possibility of exploring
#               the consequences of varying sample sizes, other than the true sample size in pn.
# m:            a numeric vector with the number of genes profiled in each column of qm.
# n0:           a numeric vector with the number of genes profiled in each column of pqn0.
# method:       the approximation method to the sampling distribution under the null hypothesis specifying that
#               the samples pn and qm come from the same population. See the 'Details' section below
# confidence:   the confidence level of the confidence intervals in the result
# d0:           the equivalence limit for the squared Euclidean distance.
#               The population profiles are considered equivalent if their true distance is less than d0
# equivEpsilon: if d0 is not given, it is computed as d0 = s * equivEpsilon^2, where s stands for the number of
#               classes or GO nodes being compared
# nsims:        some inferential methods require a simulation step; the number of simulation replicates is specified with
#               this parameter
# simplify:     should the result be simplified, if possible? See the 'Details' section
#
# Details:
# An object of S3 class 'ExpandedGOProfile' is, essentially, a 'data.frame' object
# with each column representing the relative frequencies in all observed node
# combinations, resulting from profiling a set of genes, for a given and fixed
# ontology. The row.names attribute codifies the node combinations and each
# data.frame column (say, each profile) has an attribute, 'ngenes', indicating the
# number of profiled genes.
# The arguments 'pn', 'qm' and 'pqn0' are compared in a column by column wise,
# recycling columns, if necessary, in order to perform max(ncol(pn),ncol(qm),ncol(pqn0))
# comparisons (each comparison resulting in an object of class 'htest').
# In order to be properly compared, these arguments are expanded by row, according
# to their row names. That is, the data arguments can have unequal row numbers. Then,
# they are expanded adding rows with zero frequencies, in order to make them
# comparable.
#
# In the i-th comparison (i from 1 to max(ncol(pn),ncol(qm),ncol(pqn0))),
# the parameters n, m and n0 are included to allow the possibility of exploring the consequences of varying sample
# sizes, other than the true sample sizes included as an attribute in pn, qm and pqn0.
#
# When qm = NULL, the genes profiled in pn are compared with a subsample of them, those profiled in pqn0 (compare
# a set of genes with a restricted subset, e.g. those overexpressed under a disease). When pqn0 = NULL, two profiles
# with no genes in common are compared.
#
# If p stands for the "population" profile originating the sample profile pn[,j], q for the profile originating qm[,j]
# and d(,) for the squared euclidean distance, if p != q, the distribution of sqrt(n)(d(pn[,j],qm[,j]) - d(p,q))/se is
# approximately standard normal, N(0,1). This provides the basis for the confidence interval in the result field
# icDistance.
# When p=q, the asymptotic distribution of n d(pn[,j],qm[,j]) corresponds to the distribution of a mixture of independent
# chi-square random variables with one degree of freedom. The sampling distribution under H0 p=q may be directly
# computed from this distribution (approximating it by simulation) (method="lcombChisq") or by a chi-square
# approximation to it, based on two correcting constants a and b (method="chi-square").
# These constants are chosen to equate the first two moments of both distributions (the linear
# combination of chi-square random variables distribution and the approximating chi-square distribution).
# When method="chi-square", the returned test statistic value is the chi-square approximation (n d(pn[,j],qm[,j]) - b) / a. Then,
# the result field 'parameter' is a vector containing the 'a' and 'b' values and the number of degrees of freedom, 'df'.
# Otherwise, the returned test statistic value is n d(pn[,j],qm[,j]) and 'parameter' contains the coefficients of the linear
# combination of chi-squares.
#
# The 'print' method for the resulting objects (of class 'GOProfileHtest') is the same than the 'print' method for the
# S3 class 'htest' and prints the results for the difference test H0: p==q vs. H1: p != q.
# Silently, the result of 'compareGOProfiles' contains also some fields related to the equivalence test
#                           H0: d(p,q) >= d0 vs. H1: d(p,q) < d0.
# The equivalence results are displayed by means of function 'equivalentGOProfiles'
#
# Value:
# A list containing max(ncol(pn),ncol(qm),ncol(pqn0)) objects of class 'GOProfileHtest', directly inheriting from 'htest'
# or a single 'GOProfileHtest' object if max(ncol(pn),ncol(qm),ncol(pqn0))==1 and simplify == T.
# Each object of class 'GOProfileHtest' has the following fields:
# profilePn:            the first contracted profile to compute the squared Euclidean distance
# profileQm:            the second contracted profile to compute the squared Euclidean distance
# statistic:            test statistic; its meaning depends on the value of "method", see the 'Details' section.
# parameter:            parameters of the sample distribution of the test statistic, see the 'Details' section.
# p.value:              associated p-value to test the null hypothesis of profiles equality.
# conf.int:             asymptotic confidence interval for the squared euclidean distance. Its attribute "conf.level"
#                       contains its nominal confidence level.
# estimate:             squared euclidean distance between the contracted profiles. Its attribute "se"
#                       contains its standard error estimate.
# method:               a character string indicating the method used to perform the test.
# data.name:            a character string giving the names of the data.
# alternative:          a character string describing the alternative hypothesis (always
#                       "true squared Euclidean distance between the contracted profiles is greater than zero").
# one.sided.conf.int    a one-sided confidence interval for the squared euclidean distance: [0, d.sup]. Equivalence
#                       is declared if d.sup < d0
# d0:                   the equivalence threshold for the squared Euclidean distance
#
# References: citar papers en construcció
#
# See Also:     fitGOProfile, equivalentGOProfiles
#
# Examples:
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# ATENCIO, CAL REVISAR I ACTUALITZAR ELS EXEMPLES
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  pnName <- deparse(substitute(pn))
  qmName <- deparse(substitute(qm))
  pqn0Name <- deparse(substitute(pqn0))
  pn <- as.ExpandedGOProfile(pn)
  if (!is.null(qm)) qm <- as.ExpandedGOProfile(qm)
  if (!is.null(pqn0)) pqn0 <- as.ExpandedGOProfile(pqn0)
  test.jkl <- function(i) {
    j <- i %% ncolPn + 1
    vecPn <- pn[,j]
    names(vecPn) <- rownames(pn)
    if (is.null(qm)) {
      vecQm <- NULL
      qmName.k <- NULL
    }
    else {
      k <- i %% ncolQm + 1
      vecQm <- qm[,k]
      names(vecQm) <- rownames(qm)
      qmName.k <- paste(qmName,"[",k,"]", sep="")
    }
    if (is.null(pqn0)) {
      vecPQn0 <- NULL
      pqn0Name.l <- NULL
    }
    else {
      l <- i %% ncolPQn0 + 1
      vecPQn0 <- pqn0[,l]
      names(vecPQn0) <- rownames(pqn0)
      pqn0Name.l <- paste(pqn0Name,"[",l,"]", sep="")
    }
    result.jkl <- internal.compareGOProf(vecPn, vecQm, vecPQn0, n[j], m[k], n0[l],
        method, ab.approx, confidence, nsims)
    result.jkl$data.name <-
        paste(pnName,"[",j,"] and ", qmName.k, " and ", pqn0Name.l, sep="")
    result.jkl
  }
  maxncol <- max(ncolPn <- ncol(pn),ncolQm <- ncol(qm),ncolPQn0 <- ncol(pqn0))
  result <- lapply(0:(maxncol-1), test.jkl)
  if (simplify && (maxncol == 1)) {
    result[[1]]$data.name <- paste(pnName," and ",qmName," and ",pqn0Name, sep="")
    return(result[[1]])
  }
  else
    return(result)
}