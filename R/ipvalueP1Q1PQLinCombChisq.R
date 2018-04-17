`ipvalueP1Q1PQLinCombChisq` <-
function(i, d2, pCommon, n, m, pq, n0, lambda, coefStat, sigma0, betas)
{
    pvalueIntersectLinCombChisq(d2=d2[i], pCommon=pCommon[,i], n=n, m=m, pq=pq[,i], n0=n0,
        lambda=lambda, coefStat=coefStat, sigma0=sigma0, betas=betas)

}

