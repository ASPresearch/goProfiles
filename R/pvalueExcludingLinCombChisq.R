`pvalueExcludingLinCombChisq` <-
function(d2, pCommon, n=attr(pCommon,"ngenes"), m,
    lambda = n/(n+m), coefStat = m*lambda,
    sigma = internal.covGO(pCommon), betas = NULL, attrs = F)
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
    # if (is.null(chisMat))
    #     pvalue <- plcombChisq(coefStat*d2, betas, lower.tail = F, nsims=nsims)
    # else
    #     pvalue <- sum( as.numeric( betas %*% chisMat) > coefStat * d2) / ncol(chisMat)
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

