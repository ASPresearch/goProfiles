`pvalueRestrictLinCombChisq` <-
function(d2, pCommon, n=attr(pCommon,"ngenes"), m,
    lambda = n/(n+m), coefStat = m*lambda,
    sigma = internal.covGO(pCommon), 
    betas=eigen((2*lambda-1)*sigma, symmetric=T, only.values=T)$values[1:ncol(sigma)],
    attrs=F)
{
  davies.out <- davies(coefStat*d2, betas)
  if (davies.out$ifault > 0) {
    warning(
      paste("Davies's method returned fault indicator ifault = ", davies.out$ifault, 
            ", see ?davies for more information", sep = ""))
  }
  pvalue <- max(davies.out$Qq, 2.2e-16)
  if (attrs) attr(pvalue,"betas") <- betas
    return(pvalue)
}

