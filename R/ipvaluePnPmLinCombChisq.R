`ipvaluePnPmLinCombChisq` <-
function(i, d2, pCommon, n, m, sigma, lambda, coefStat)
{
    pvalueRestrictLinCombChisq(d2=d2[i], pCommon=pCommon[,i], n=n, m=m, lambda=lambda, coefStat=coefStat)
}

