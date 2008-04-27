
chiRestrict <- function(d2, pCommon, n=attr(pCommon,"ngenes"), m,
    coefStat = m*n/(m+n),
    sigma=NULL, betas=NULL,
    dfr = NULL,
    ab.approx = "asymptotic",
    nsims=10000)
{
    if (is.null(dfr))
        dfr <- sum(contractedProfile(pCommon) > 0)
    if (pmatch(ab.approx,"asymptotic",nomatch=F)) {
        if (is.null(sigma)) sigma <- internal.covGO(pCommon)
        if (is.null(betas)) betas <- eigen(((n-m)/(n+m)) * sigma, symmetric=T, only.values=T)$values#[1:dfr]
        mean.stat <- sum(betas)
        var.stat <- 2 * sum(betas*betas)
    }
    else {
        mnom <- multinomial(n=m, p=pCommon)
        pmSampl <- rndGO(mnom, ssize=nsims)
        dsamp <- sapply(1:nsims, idEuclid2, pmSampl, expandRndGO(pmSampl, setN.multinomial(mnom,n=n-m)))
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

