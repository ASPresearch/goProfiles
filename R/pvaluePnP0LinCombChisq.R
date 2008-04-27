`pvaluePnP0LinCombChisq` <-
function(d2, p0, n=attr(p0,"ngenes"), sigma = internal.covGO(p0), betas = eigen(sigma, symmetric=T, only.values=T)$values, 
    nsims = 10000, chisMat=NULL, attrs=F)
{
    if (is.null(chisMat))
        pvalue <- plcombChisq(n * d2, betas, lower.tail = F, nsims=nsims)
    else
        pvalue <- sum( as.numeric( betas %*% chisMat) > n * d2) / ncol(chisMat)
    if (attrs) attr(pvalue,"betas") <- betas
    return(pvalue)
}

