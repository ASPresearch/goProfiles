
internal.compareGOProf <- function(pn, qm, pqn0, n, m, n0,
  method, ab.approx, confidence)
{
# Compare two "sample" GO profiles pn and qm. Are they just random samples taken
# from the same "population" GO profile?
# Description:
# 'internal.compareGOProf' implements some inferential procedures to solve the
# preceding question. These procedures are based on
# asymptotic properties of the squared euclidean distance between the contracted
# versions of pn and qm
#
# Usage:
#   internal.compareGOProf(pn, qm=NULL, pqn0=NULL, n = attr(pn,"ngenes"), m = attr(qm,"ngenes"), n0 = attr(pqn0,"ngenes"),
#       method = "lcombChisq", confidence = 0.95)
#
#
# Arguments:
# pn:           a numeric vector representing a "sample" expanded GO profile.
#               Its 'names' attribute should codify the node combinations of the
#               profiled genes (in a given ontology) and, additionally, it must
#               have a 'ngenes' attribute specifying the number of profiled genes.
#               The values in 'pn' are interpreted as relative frequencies
# qm:           a numeric vector representing a "sample" expanded GO profile.
#               The comments to pn are also valid.
# pqn0:         a numeric vector representing a "sample" expanded GO profile
#               corresponding to the common genes in pn and qm
# n:            the number of genes in the sample pn if it should differ from its 'ngenes' attribute
# m:            the number of genes in the sample qm if it should differ from its 'ngenes' attribute
# n0:           the number of genes in the sample pqn0 if it should differ from its 'ngenes' attribute.
# method:       the approximation method to the sampling distribution under the null hypothesis specifying that
#               the samples pn and qm come from the same population. See the 'Details' section below
# confidence:   the confidence level of the confidence interval in the result
#
# Details:
# The parameters n, m and n0 are included to allow the possibility of exploring the consequences of varying sample
# sizes, other than the true sample sizes included as an attribute in pn, qm and pqn0.
#
# When qm = NULL, the genes profiled in pn are compared with a subsample of them, those profiled in pqn0 (compare
# a set of genes with a restricted subset, e.g. those overexpressed under a disease). When pqn0 = NULL, two profiles
# with no genes in common are compared.
#
# If p stands for the "population" profile originating the sample profile pn, q for the profile originating qm
# and d(,) for the squared euclidean distance, if p != q, the distribution of sqrt(n)(d(pn,qm) - d(p,q))/se is
# approximately standard normal, N(0,1). This provides the basis for the confidence interval in the result field
# icDistance.
# When p=q, the asymptotic distribution of n d(pn,qm) corresponds to the distribution of a mixture of independent
# chi-square random variables with one degree of freedom. The sampling distribution under H0 p=q may be directly
# computed from this distribution (approximating it by simulation) (method="lcombChisq") or by a chi-square
# approximation to it, based on two correcting constants a and b (method="chi-square").
# These constants are chosen to equate the first two moments of both distributions (the linear
# combination of chi-square random variables distribution and the approximating chi-square distribution).
# When method="chi-square", the returned test statistic value is the chi-square approximation (n d(pn,qm) - b) / a. Then,
# the result field 'parameter' is a vector containing the 'a' and 'b' values and the number of degrees of freedom, 'df'.
# Otherwise, the returned test statistic value is n d(pn,qm) and 'parameter' contains the coefficients of the linear
# combination of chi-squares.
#
# Value:
# An object of class 'GOProfileHtest', containing the following fields:
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

# References: citar papers en construccio
#
# See Also:     internal.fitGOProf, simPQIntersect, simSeriesPQIntersect
#
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# ATENCIO, CAL REVISAR I ACTUALITZAR ELS EXEMPLES
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Examples:
# dataPath <- ".\\sampleDatasets\\"  # put here the adequate path
# TCells_vec_BPprofile <- dgetGO(dataPath,"TCells_vec_BPprofile.txt")
# hgu95a_vec_BPprofile <- dgetGO(dataPath,"hgu95a_vec_BPprofile.txt")
#
# TCells_vec_MFprofile <- dgetGO(dataPath,"TCells_vec_MFprofile.txt")
# hgu95a_vec_MFprofile <- dgetGO(dataPath,"hgu95a_vec_MFprofile.txt")
#
# TCells_vec_CCprofile <- dgetGO(dataPath,"TCells_vec_CCProfile.txt")
# hgu95a_vec_CCprofile <- dgetGO(dataPath,"hgu95a_vec_CCProfile.txt")
#
# Compare the full microarray profile 'hgu95a_vec_BPprofile' with a subset of it 'TCells_vec_BPprofile'
# TCellsBPcompare.moments.asymp <- compareGOProfiles(hgu95a_vec_BPprofile, pqn0=TCells_vec_BPprofile, method="chi-square")
#
# TCellsBPcompare.linCombChisq <- compareGOProfiles(hgu95a_vec_BPprofile, pqn0=TCells_vec_BPprofile, method="lcombChisq")
#
# Compare the profiles of two disjoint sets of genes
# notTCells_vec_BPprofile <- (attr(hgu95a_vec_BPprofile,"ngenes") * hgu95a_vec_BPprofile - attr(TCells_vec_BPprofile,"ngenes") * fullGOProfile(TCells_vec_BPprofile,names(hgu95a_vec_BPprofile))) / (attr(hgu95a_vec_BPprofile,"ngenes") - attr(TCells_vec_BPprofile,"ngenes"))
# attr(notTCells_vec_BPprofile,"ngenes") <- attr(hgu95a_vec_BPprofile,"ngenes") - attr(TCells_vec_BPprofile,"ngenes")
#
# TCellsBPcomparePQ.moments.asymp <- compareGOProfiles(notTCells_vec_BPprofile, TCells_vec_BPprofile, method="chi-square")
# # Equivalent to:
# TCellsBPcomparePQ.moments.asymp <- compareGOProfiles(pn=notTCells_vec_BPprofile, qm=TCells_vec_BPprofile, method="chi-square")
# TCellsBPcomparePQ.linCombChisq <- compareGOProfiles(notTCells_vec_BPprofile, TCells_vec_BPprofile, method="lcombChisq")
#
# Compare two (simulated) profiles sharing 100 genes:
# profilP0 <- generate.multinomial(multinomial(n=100,p=hgu95a_vec_BPprofile))
# profilP <- profilP0 + generate.multinomial(multinomial(n=200,p=hgu95a_vec_BPprofile))
# profilQ <- profilP0 + generate.multinomial(multinomial(n=300,p=hgu95a_vec_BPprofile))
# profilP0 <- profilP0 / 100
# profilP <- profilP / 300
# profilQ <- profilQ / 400
# profilP <- as.vector(profilP)
# names(profilP) <- colnames(profilQ)
# profilQ <- as.vector(profilQ)
# names(profilQ) <- names(profilP)
# profilP0 <- as.vector(profilP0)
# names(profilP0) <- names(profilP)
# attr(profilP0, "ngenes") <- 100
# attr(profilP, "ngenes") <- 300
# attr(profilQ, "ngenes") <- 400
#
#
# compareGOProfiles(pn=profilP, qm=profilQ, pqn0=profilP0, method="chi-square")
# compareGOProfiles(pn=profilP, qm=profilQ, pqn0=profilP0, method="lcombChisq")
#
#
    dataNames <- paste(deparse(substitute(pn)), " and ", deparse(substitute(qm)), " and ", deparse(substitute(pqn0)), sep="")
    fullNams <- unique(c(names(pn),names(qm),names(pqn0)))
    pn <- fullGOProfile(pn, fullNams)
    qm <- fullGOProfile(qm, fullNams)
    pqn0 <- fullGOProfile(pqn0, fullNams)
    contrPn <- contractedProfile.default(pn)
#    contrPn0 <- contrQm <- NULL

    if (pmatch(method,"chi-square",nomatch=F)) {
        if (is.null(qm)) {
            # compare a profile with a subset of it
            d <- dEuclid2(contrPn, contrQm <- contractedProfile.default(pqn0))
            #d <- dEuclid2(contrPn, contrPn0 <- contractedProfile.default(pqn0))
            stat <- chiRestrict(d2=d, pCommon=pn, n=n, m=n0, ab.approx=ab.approx)
            m <- n0
        } else if (is.null(pqn0)) {
            # compare two disjoint profiles (no gens in common)
            d <- dEuclid2(contrPn, contrQm <- contractedProfile.default(qm))
            stat <- chiDisjoint(d2=d, pCommon=(n*pn+m*qm)/(n+m), n=n, m=m, ab.approx=ab.approx)
        } else {
            # compare two intersecting profiles (some genes are specific of pn, some are specific of qm, and some common genes are profiled in pqn0)
            d <- dEuclid2(contrPn, contrQm <- contractedProfile.default(qm))
            stat <- chiIntersect(d2=d, pCommon=(n*pn+m*qm-n0*pqn0)/(n+m-n0),
                n=n, m=m, pq=pqn0, n0=n0,
                ab.approx=ab.approx)
            #contrPn0 <- contractedProfile.default(pqn0)
        }
        dgf <- attr(stat,"df")
        parameters <- c(attr(stat,"a"),attr(stat,"b"),dgf)
        names(parameters) <- c("a","b","df")
        attr(method,"argument") <- method
        method <- "chi-square statistic (affine approximation n*d2 aprox. a*X2+b)"
        pvalue <- pchisq(stat,df=attr(stat,"df"),lower.tail=F)
    } else {
        if (is.null(qm)) {
            d <- dEuclid2(contrPn, contrQm <- contractedProfile.default(pqn0))
            pvalue <- pvalueRestrictLinCombChisq(d, pCommon=pn, n=n, m=n0, attrs=T)
            stat <- (n*n0/(n+n0)) * d
            names(stat) <- "(n*n0/(n+n0)) * d2"
            m <- n0
        } else if (is.null(pqn0)) {
            d <- dEuclid2(contrPn, contrQm <- contractedProfile.default(qm))
            pvalue <- pvalueExcludingLinCombChisq(d, pCommon=(n*pn+m*qm)/(n+m), n=n, m=m, attrs=T)
            stat <- (n*m/(n+m)) * d
            names(stat) <- "(n*m/(n+m)) * d2"
        }
        else {
            d <- dEuclid2(contrPn, contrQm <- contractedProfile.default(qm))
            pvalue <- pvalueIntersectLinCombChisq(d, pCommon=(n*pn+m*qm-n0*pqn0)/(n+m-n0), n=n, m=m, pq=pqn0, n0=n0, attrs=T)
            #contrPn0 <- contractedProfile.default(pqn0)
            stat <- (n*m/(n+m)) * d
            names(stat) <- "(n*m/(n+m)) * d2"
        }
        s <- sum((contrPn[,3] + contrQm[,3]) > 0) # classes being compared in the profiles
        betas <- attr(pvalue,"betas")
        parameters <- c(s, betas)
        names(parameters) <- c("number of classes", paste("coef",1:length(betas),sep=""))
        attr(method,"argument") <- method
        method <- "linear combination of chi-squares statistic"
    }

    names(d) <- "sample squared Euclidean distance"
    attributes(pvalue) <- NULL
    se <- sqrt(internal.varDifStatPnQm(pn=pn, qm=qm, pqn0 = pqn0) * ((n+m)/(n*m)))
    names(se) <- "distance standard error"
    attr(d,"se") <- se
    p <- (1-confidence)/2
    d.precision <- qnorm(p=p, lower.tail = F) * se
    icDistance <- c(d - d.precision, d + d.precision)
    names(icDistance) <- c(paste("two-sided lower",p), paste("two-sided upper",p))
    attr(icDistance,"conf.level") <- confidence
    result <- list(
        profilePn = contrPn,
        profileQm = contrQm,
        statistic = stat, parameter = parameters, p.value = pvalue, conf.int = icDistance, estimate = d,
        method = method, data.name = dataNames,
        alternative = "true squared Euclidean distance between the contracted profiles is greater than zero"
    )
    class(result) <- c("GOProfileHtest", "htest")
    return(result)
}
