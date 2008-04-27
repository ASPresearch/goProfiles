`pvalueExcludingLinCombChisq` <-
function(d2, pCommon, n=attr(pCommon,"ngenes"), m,
    lambda = n/(n+m), coefStat = m*lambda,
    sigma = internal.covGO(pCommon), betas=NULL,
    nsims = 10000, chisMat=NULL, attrs=F)
{
    if (is.null(betas)) {
      k <- ncol(sigma)
      betas <- eigen(
        cbind(rbind(diag(k),-diag(k)),rbind(-diag(k),diag(k))) 
        %*%
        cbind(rbind((1-lambda)*sigma,matrix(0,nrow=k,ncol=k)),rbind(matrix(0,nrow=k,ncol=k),lambda*sigma)),
        symmetric=T,
        only.values=T)$values[1:k]
    }
    if (is.null(chisMat))
        pvalue <- plcombChisq(coefStat*d2, betas, lower.tail = F, nsims=nsims)
    else
        pvalue <- sum( as.numeric( betas %*% chisMat) > coefStat * d2) / ncol(chisMat)
    if (attrs) attr(pvalue,"betas") <- betas
    return(pvalue)
}

