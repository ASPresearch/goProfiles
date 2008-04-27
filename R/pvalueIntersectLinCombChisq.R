`pvalueIntersectLinCombChisq` <-
function(d2, pCommon, n=attr(pCommon,"ngenes"), m, pq, n0=attr(pq,"ngenes"),
    lambda = n/(n+m), coefStat = m*lambda,
    sigma0=internal.covGO(pq), betas=NULL,
    nsims = 10000, chisMat=NULL, attrs=F)
{
    if (is.null(betas)) {
        n1 <- n - n0
        m1 <- m - n0
        p1 <- (n/n1)*pCommon - (n0/n1)*pq
        p1[p1 < 0] <- 0
        p1 <- p1 / sum(p1)
        sigma.p <- internal.covGO(p1)
        q1 <- (m/m1)*pCommon - (n0/m1)*pq
        q1[q1 < 0] <- 0
        q1 <- q1 / sum(q1)
        sigma.q <- internal.covGO(q1)
        k <- ncol(sigma0)
        theta <- n0 / (n + m)
        betas <- eigen(
            cbind(rbind(diag(k),-diag(k)),rbind(-diag(k),diag(k))) 
            %*%
            cbind(rbind((1-lambda)*((n0/n)*sigma0+(n1/n)*sigma.p), theta*sigma0),
                  rbind(theta*sigma0, lambda*((n0/m)*sigma0+(m1/m)*sigma.q))),
        symmetric=T,
        only.values=T)$values#[1:k]
    }
    if (is.null(chisMat))
        pvalue <- plcombChisq(coefStat*d2, betas, lower.tail = F, nsims=nsims)
    else
        pvalue <- sum( as.numeric( betas %*% chisMat) > coefStat * d2) / ncol(chisMat)
    if (attrs) attr(pvalue,"betas") <- betas
    return(pvalue)
}

