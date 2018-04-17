`fitGOProfile` <-
function(pn, p0, n = ngenes(pn),
    method = "lcombChisq", ab.approx = "asymptotic", confidence = 0.95, simplify=T)
{
#
# Does a "sample" GO profile 'pn', observed in a sample of 'n' genes, fit a
# "population" or "model" p0?
#
# Description:
# 'fitGOProfile' implements some inferential procedures to solve the preceding question.
# These procedures are based on asymptotic properties of the squared euclidean distance
# between the contracted versions of pn and p0
#
# Usage:
#   fitGOProfile(pn, p0, n = ngenes(pn),
#       method = "lcombChisq", ab.approx = "asymptotic", confidence = 0.95, simplify=T)
#
# Arguments:
# pn:           an object of class ExpandedGOProfile representing one or more
#               "sample" expanded GO profiles for a fixed ontology (see the 'Details' section)
# p0:           an object of class ExpandedGOProfile representing one or more
#               "population" or "theoretical" expanded GO profiles (see also the 'Details' section)
# n:            a numeric vector with the number of genes profiled in each column of pn.
#               This parameter is included to allow the possibility of exploring
#               the consequences of varying sample sizes, other than the true sample size in pn.
# method:       the approximation method to the sampling distribution under the null hypothesis "p = p0", where p is the
#               "true" population profile originating each column of pn. See the 'Details' section below
# confidence:   the confidence level of the confidence interval in the result
# simplify:     should the result be simplified, if possible? See the 'Details' section
#
# Details:
# An object of class 'ExpandedGOProfile' is, essentially, a 'data.frame' object
# with each column representing the relative frequencies in all observed node
# combinations, resulting from profiling a set of genes, for a given and fixed
# ontology. The row.names attribute codifies the node combinations and each
# data.frame column (say, each profile) has an attribute, 'ngenes', indicating the
# number of profiled genes. (Actually, the 'ngenes' attribute of each 'p0' column
# is ignored and is taken as if it were infinite, 'Inf'.)
# The arguments 'pn' and 'p0' are compared in a column by column wise, 
# recycling columns, if necessary, in order to perform max(ncol(pn),ncol(p0))
# comparisons (each comparison resulting in an object of class 'htest').  
# In order to be properly compared, 'pn' and 'p0' are expanded by row, according
# to their row names. That is, both arguments can have unequal row numbers. Then,
# they are expanded adding rows with zero frequencies, in order to make them
# comparable.

# In the i-th comparison (i from 1 to max(ncol(pn),ncol(p0))),
# if p stands for the profile originating the sample profile pn[,i] and d(,)
# for the squared euclidean distance, 
# if p != p0[,i], the distribution of
# sqrt(n)(d(pn[,i],p0[,i]) - d(p,p0[,i]))/se
# is approximately standard normal, N(0,1). This provides
# the basis for the confidence interval in the result field conf.int.
# When p==p0[,i], the asymptotic distribution of
# n d(pn[,i],p0[,i]) is the distribution of a linear combination of independent
# chi-square random variables, each one with one degree of freedom.
# This sampling distribution may be directly computed (approximating it by simulation, method="lcombChisq")
# or approximated by a chi-square distribution, based on two correcting constants a and b (method="chi-square").
# These constants are chosen to equate the first two moments of both distributions (the distribution of a linear combination
# of chi square variables and the approximating chi-square distribution).
# When method="chi-square", the returned test statistic value is the chi-square approximation (n d(pn,p0) - b) / a. Then,
# the result field 'parameter' is a vector containing the 'a' and 'b' values and the number of degrees of freedom, 'df'.
# Otherwise, the returned test statistic value is n d(pn,p0) and 'parameter' contains the coefficients of the linear
# combination of chi-squares.
#
# Value:
# A list containing max(ncol(pn),ncol(p0)) objects of class 'htest', 
# or a single 'htest' object if ncol(pn)==1 and ncol(p0)==1 and simplify == T.
# Each 'htest' object has the following fields:
# statistic:            test statistic; its meaning depends on the value of "method", see the 'Details' section.
# parameter:            parameters of the sample distribution of the test statistic, see the 'Details' section.
# p.value:              associated p-value to test the null hypothesis "pn[,i] is a random sample taken from p0[,i]".
# conf.int:             asymptotic confidence interval for the squared euclidean distance. Its attribute "conf.level"
#                       contains its nominal confidence level.
# estimate:             squared euclidean distance between the contracted pn and p0 profiles. Its attribute "se"
#                       contains its standard error estimate.
# method:               a character string indicating the method used to perform the test.
# data.name:            a character string giving the names of the data.
# alternative:          a character string describing the alternative hypothesis (always
#                       "true squared Euclidean distance between the contracted profiles is greater than zero").
#
# References: citar papers en construccio
#
# See Also:     compareGOProfiles, simPnP0, simSeriesPnP0
#
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# ATENCIO, CAL REVISAR I ACTUALITZAR ELS EXEMPLES
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Examples:
# dataPath <- ".\\sampleDatasets\\"
# TCells_vec_BPprofile <- dgetGO(dataPath,"TCells_vec_BPprofile.txt")
# hgu95a_vec_BPprofile <- dgetGO(dataPath,"hgu95a_vec_BPprofile.txt")
#
# TCells_vec_MFprofile <- dgetGO(dataPath,"TCells_vec_MFprofile.txt")
# hgu95a_vec_MFprofile <- dgetGO(dataPath,"hgu95a_vec_MFprofile.txt")
#
# TCells_vec_CCprofile <- dgetGO(dataPath,"TCells_vec_CCProfile.txt")
# hgu95a_vec_CCprofile <- dgetGO(dataPath,"hgu95a_vec_CCProfile.txt")
# 
# internal.fitGOProf(TCells_vec_BPprofile, hgu95a_vec_BPprofile, method="chi-square")
# # Equivalent to method="lcombChisq":
# internal.fitGOProf(TCells_vec_BPprofile, hgu95a_vec_BPprofile)
# The same comparison but assuming a sample size "ngenes" for TCells_vec_BPProfile of 1000 (instead of its 140 true value):
# internal.fitGOProf(TCells_vec_BPprofile, hgu95a_vec_BPprofile, 1000)
# # Then, a significant result!!
#
# # Or equivalently:
# internal.fitGOProf(TCells_vec_BPprofile, hgu95a_vec_BPprofile, n=1000)
# Include the original expanded GO profiles in the result (displaying all GO categories, even if they are empty):
# internal.fitGOProf(TCells_vec_BPprofile, hgu95a_vec_BPprofile, dataInResult=T)
#
# internal.fitGOProf(TCells_vec_MFprofile, hgu95a_vec_MFprofile, method="chi-square")
# internal.fitGOProf(pn=TCells_vec_MFprofile, p0=hgu95a_vec_MFprofile)
# internal.fitGOProf(pn=TCells_vec_MFprofile, p0=hgu95a_vec_MFprofile, n=1000)
#
# internal.fitGOProf(TCells_vec_CCprofile, hgu95a_vec_CCprofile, method="chi-square")
# internal.fitGOProf(pn=TCells_vec_CCprofile, p0=hgu95a_vec_CCprofile)
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#
#
  pnName <- deparse(substitute(pn))
  p0Name <- deparse(substitute(p0))
  pn <- as.ExpandedGOProfile(pn)
  p0 <- as.ExpandedGOProfile(p0)
  test.jk <- function(i) {
    j <- i %% ncolPn + 1
    k <- i %% ncolP0 + 1
    vecPn <- pn[,j]
    names(vecPn) <- rownames(pn)
    vecP0 <- p0[,k]
    names(vecP0) <- rownames(p0)
    result.jk <- internal.fitGOProf(vecPn, vecP0, n[j], method, ab.approx,
      confidence)
    result.jk$data.name <- paste(pnName,"[",j,"] and ",p0Name,"[",k,"]", sep="")
    result.jk
  }
  maxncol <- max(ncolPn <- ncol(pn),ncolP0 <- ncol(p0))
  result <- lapply(0:(maxncol-1), test.jk)
  if (simplify && (maxncol == 1)) {
    result[[1]]$data.name <- paste(pnName," and ",p0Name, sep="")
    return(result[[1]])
  }
  else
    return(result)
}

