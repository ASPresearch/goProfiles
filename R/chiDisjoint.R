
chiDisjoint <- function(d2, pCommon, n=attr(pCommon,"ngenes"), m,
    lambda = n/(n+m), coefStat = m*lambda,
    sigma = internal.covGO(pCommon), betas = NULL,
    dfr = NULL,
    ab.approx = "asymptotic",
    nsims=10000)
{
    if (is.null(dfr))
        dfr <- sum(contractedProfile(pCommon) > 0)
    if (pmatch(ab.approx,"asymptotic",nomatch=F)) {
        if (is.null(betas)) {
          k <- ncol(sigma)
          betas <- eigen(
            cbind(rbind(diag(k),-diag(k)),rbind(-diag(k),diag(k)))
            %*%
            cbind(rbind((1-lambda)*sigma,matrix(0,nrow=k,ncol=k)),rbind(matrix(0,nrow=k,ncol=k),lambda*sigma)),
            symmetric=T,
            only.values=T)$values#[1:dfr]
        }
        mean.stat <- sum(betas)
        var.stat <- 2 * sum(betas*betas)
    }
    else {
        mnom1 <- multinomial(n=n, p=pCommon)
        mnom2 <- multinomial(n=m, p=pCommon)
        dsamp <- sapply(1:nsims, idEuclid2, rndGO(mnom1, ssize=nsims), rndGO(mnom2, ssize=nsims))
        mean.stat <- coefStat * mean(dsamp)
        var.stat <- (coefStat*coefStat) * var(dsamp)
    }
    a <- sqrt(var.stat/(2*dfr))
    b <- mean.stat - a*dfr
    result <- (coefStat * d2 - b) / a
    attr(result,"a") <- a
    attr(result,"b") <- b
    attr(result,"df") <- dfr
    return (result)
}

