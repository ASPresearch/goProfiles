`pvaluePnP0LinCombChisq` <-
function(d2, p0, n=attr(p0,"ngenes"), sigma = internal.covGO(p0), 
         betas = eigen(sigma, symmetric=T, only.values=T)$values, attrs=F)
{
    # if (is.null(chisMat))
    #     pvalue <- plcombChisq(n * d2, betas, lower.tail = F, nsims=nsims)
    # else
    #     pvalue <- sum( as.numeric( betas %*% chisMat) > n * d2) / ncol(chisMat)
  davies.out <- davies(n * d2, betas)
  if (davies.out$ifault > 0) {
    warning(
      paste("Davies's method returned fault indicator ifault = ", davies.out$ifault, 
            ", see ?davies for more information", sep = ""))
  }
  pvalue <- max(davies.out$Qq, 2.2e-16)
    if (attrs) attr(pvalue,"betas") <- betas
    return(pvalue)
}

