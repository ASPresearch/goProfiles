`internal.fitGOProf` <-
function(pn, p0, n = attr(pn,"ngenes"),
    method = "lcombChisq", ab.approx = "asymptotic", confidence = 0.95)
{
# Is the "sample" GO profile pn just a random sample taken from the "population" GO profile p0?
# Description:
# 'internal.fitGOProf' implements some inferential procedures to solve the preceding question. These procedures are based on
# asymptotic properties of the squared euclidean distance between the contracted versions of pn and p0
#
# Usage:
#   internal.fitGOProf(pn, p0, n = attr(pn,"ngenes"), method = "lcombChisq", confidence = 0.95)
#
# Arguments:
# pn:           a numeric vector representing a "sample" expanded GO profile.
#               Its 'names' attribute should codify the node combinations of the
#               profiled genes (in a given ontology) and, additionally, it must
#               have a 'ngenes' attribute specifying the number of profiled genes.
#               The values in 'pn' are interpreted as relative frequencies
# p0:           a vector representing a "population" or "theoretical" expanded GO profile.
#               It must have the same structure than pn ('names' and 'ngenes' attributes,
#               etc) but the value of the attribute 'ngenes' is internally taken as
#               the infinity, 'Inf'
# n:            the number of genes in the sample pn. This parameter is included to allow the possibility of exploring
#               the consequences of varying sample sizes, other than the true sample size in the attribute 'ngenes' of pn.
# method:       the approximation method to the sampling distribution under the null hypothesis "p = p0", where p is the
#               "true" population profile originating pn. See the 'Details' section below
# confidence:   the confidence level of the confidence interval in the result
#
# Details:
# If p stands for the profile originating the sample profile pn and d(,) for the squared euclidean distance, 
# if p != p0, the distribution of sqrt(n)(d(pn,p0) - d(p,p0))/se is approximately standard normal, N(0,1). This provides
# the basis for the confidence interval in the result field conf.int.
# When p=p0, the asymptotic distribution of
# n d(pn,p0) is the distribution of a linear combination of independent chi-square random variables,
# each one with one degree of freedom.
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
# An object of class 'htest', containing the following fields:
# statistic:            test statistic; its meaning depends on the value of "method", see the 'Details' section.
# parameter:            parameters of the sample distribution of the test statistic, see the 'Details' section.
# p.value:              associated p-value to test the null hypothesis "pn is a random sample taken from p0".
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
# Include the original expanded GO profiles in the result (displaying all GO categories, even if the are empty):
# internal.fitGOProf(TCells_vec_BPprofile, hgu95a_vec_BPprofile, dataInResult=T)
#
# internal.fitGOProf(TCells_vec_MFprofile, hgu95a_vec_MFprofile, method="chi-square")
# internal.fitGOProf(pn=TCells_vec_MFprofile, p0=hgu95a_vec_MFprofile)
# internal.fitGOProf(pn=TCells_vec_MFprofile, p0=hgu95a_vec_MFprofile, n=1000)
#
# internal.fitGOProf(TCells_vec_CCprofile, hgu95a_vec_CCprofile, method="chi-square")
# internal.fitGOProf(pn=TCells_vec_CCprofile, p0=hgu95a_vec_CCprofile)
#
#
    dataNames <- paste(deparse(substitute(pn)), " and ", deparse(substitute(p0)), sep="")
    fullNams <- unique(c(names(pn),names(p0)))
    pn <- fullGOProfile(pn, fullNams)
    p0 <- fullGOProfile(p0, fullNams)
    d <- dEuclid2(contrPn <- contractedProfile(pn), contractedProfile(contrP0 <- p0))     
#    dims <- dimsGO(c(names(pn),names(p0)))
#    ncateg <- dims$ncateg
#    simult <- dims$simult
    names(d) <- "sample squared Euclidean distance"
    attr(contrPn,"ngenes") <- n
    attr(contrP0,"ngenes") <- Inf
    if (pmatch(method,"chi-square",nomatch=F)) {
        stat <- chiPnP0Correct(d2=d, p0=p0, n=n, ab.approx=ab.approx)
        names(stat) <- "(n*d2 - b) / a"
        parameters <- c(attr(stat,"a"),attr(stat,"b"),attr(stat,"df"))
        names(parameters) <- c("a","b","df")
        attr(method,"argument") <- method
        method <- "chi-square statistic (affine approximation n*d2 aprox. a*X2+b)"
        pvalue <- pchisq(stat,df=attr(stat,"df"),lower.tail=F)
    }
    else {
        stat <- n * d
        names(stat) <- "n*d2"
        #parameters <- eigen(internal.covGO(p0), symmetric=T, only.values=T)$values
        #names(parameters) <- paste("coef",1:length(parameters),sep="")
        attr(method,"argument") <- method
        method <- "linear combination of chi-squares statistic"
        #pvalue <- pvaluePnP0LinCombChisq(d, p0=p0, n=n, betas=parameters)
        pvalue <- pvaluePnP0LinCombChisq(d, p0=p0, n=n, attrs=T)
        parameters <- attr(pvalue,"betas")
        names(parameters) <- paste("coef",1:length(parameters),sep="")
    }
    attributes(pvalue) <- NULL
    se <- sqrt(internal.varDifStatPP0(pn=pn, p0 = p0) / n)
    names(se) <- "distance standard error"
    attr(d,"se") <- se
    d.precision <- qnorm(p=(1-confidence)/2, lower.tail = F) * se
    icDistance <- c(d - d.precision, d + d.precision)
    attr(icDistance,"conf.level") <- confidence
    names(icDistance) <- NULL
    result <- list(
        statistic = stat, parameter = parameters, p.value = pvalue, conf.int = icDistance, estimate = d,
        method = method, data.name = dataNames,
        alternative = "true squared Euclidean distance between the contracted profiles is greater than zero"
    )
    class(result) <- "htest"
    return(result)
}

