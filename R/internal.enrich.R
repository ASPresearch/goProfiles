
internal.enrich <- function(pn, qm, pqn0, n, m, n0, method, ...)
{
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
# method:       corresponds to the 'p.adjust' function argument 'method'
# Details:
# The parameters n, m and n0 are included to allow the possibility of exploring the consequences of varying sample
# sizes, other than the true sample sizes included as an attribute in pn, qm and pqn0.
#
# When qm = NULL, the genes profiled in pn are compared with a subsample of them, those profiled in pqn0 (compare
# a set of genes with a restricted subset, e.g. those overexpressed under a disease). When pqn0 = NULL, two profiles
# with no genes in common are compared.
#
    dataNames <- paste(deparse(substitute(pn)), " and ", deparse(substitute(qm)), " and ", deparse(substitute(pqn0)), sep="")
    fullNams <- unique(c(names(pn),names(qm),names(pqn0)))
    pn <- fullGOProfile(pn, fullNams)
    qm <- fullGOProfile(qm, fullNams)
    pqn0 <- fullGOProfile(pqn0, fullNams)
    contrPn <- contractedProfile.default(pn)

    if (is.null(qm)) {
        # compare a profile with a subset of it
        contrPQn0 <- contractedProfile.default(pqn0)
        internal.enrichProfile(contrPn * n - contrPQn0 * n0, contrPQn0 * n0, n - n0, n0, method, ...)
    }
    else if (is.null(pqn0)) {
        # compare two disjoint profiles (no genes in common)
        contrQm <- contractedProfile.default(qm)
        internal.enrichProfile(contrPn * n, contrQm * m, n, m, method, ...)
    }
    else {
        # compare two intersecting profiles (n genes are specific of pn, m are specific of qm,
        # and n0 common genes are profiled in pqn0)
        internal.enrichProfile(contrPn * n - contractedProfile.default(pqn0) * n0, contractedProfile.default(qm) * m, n - n0, m, method, ...)
    }
}
