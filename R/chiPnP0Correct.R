`chiPnP0Correct` <-
function(d2, p0, n=attr(p0,"ngenes"),
    sigma = internal.covGO(p0), betas = eigen(sigma, symmetric=T, only.values=T)$values, 
    pp0=NULL, dfr = NULL,
    ab.approx = "asymptotic", 
    nsims=10000)
{
    if (pmatch(ab.approx,"asymptotic",nomatch=F)) {
        mean.stat <- sum(betas)
        var.stat <- 2 * sum(betas*betas)
        if (is.null(dfr)) dfr <- sum(betas > 0)
    }
    else {
        if (is.null(pp0)) pp0 <- contractedProfile(p0)
        if (is.null(dfr)) dfr <- sum(pp0 > 0)
        mnom <- multinomial(n=n, p=p0)
        pnSampl <- rndGO(mnom, ssize=nsims)
        dsamp <- sapply(1:nsims, idEuclid2P0, pnSampl, pp0)
        mean.stat <- n * mean(dsamp)
        var.stat <- (n*n) * var(dsamp)
    }
    a <- sqrt(var.stat/(2*dfr))
    b <- mean.stat - a*dfr
    result <- (n * d2 - b) / a
    attr(result,"a") <- a
    attr(result,"b") <- b
    attr(result,"df") <- dfr
    return (result)
}

