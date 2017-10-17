`pvalueIntersectLinCombChisq` <-
function(d2, pCommon, n=attr(pCommon,"ngenes"), m, pq, n0=attr(pq,"ngenes"),
    lambda = n/(n+m), coefStat = m*lambda,

    sigma0=internal.covGO(pCommon), betas=NULL,
    nsims = 10000, chisMat=NULL, attrs=F)
{
    if (is.null(betas)) {
        n1 <- n - n0
        m1 <- m - n0
        k <- ncol(sigma0)
        theta <- n0 / (n + m)
        betas <- eigen(
            cbind(rbind(diag(k),-diag(k)),rbind(-diag(k),diag(k))) 
            %*%
            cbind(rbind((1-lambda)*sigma0, theta*sigma0),
                  rbind(theta*sigma0, lambda*sigma0)),
        symmetric=T,
        only.values=T)$values
    }
    if (is.null(chisMat))
        pvalue <- plcombChisq(coefStat*d2, betas, lower.tail = F, nsims=nsims)
    else
        pvalue <- sum( as.numeric( betas %*% chisMat) > coefStat * d2) / ncol(chisMat)
    if (attrs) attr(pvalue,"betas") <- betas
    return(pvalue)
}

