`pvalueRestrictLinCombChisq` <-
function(d2, pCommon, n=attr(pCommon,"ngenes"), m,
    lambda = n/(n+m), coefStat = m*lambda,
    sigma = internal.covGO(pCommon), betas=eigen((2*lambda-1)*sigma, symmetric=T, only.values=T)$values[1:ncol(sigma)],
    nsims = 10000, chisMat=NULL, attrs=F)
{
    if (is.null(chisMat))
        pvalue <- plcombChisq(coefStat * d2, betas, lower.tail = F, nsims=nsims)
    else
        pvalue <- sum( as.numeric( betas %*% chisMat) > coefStat * d2) / ncol(chisMat)
    if (attrs) attr(pvalue,"betas") <- betas
    return(pvalue)
}

