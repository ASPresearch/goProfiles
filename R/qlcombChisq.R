`qlcombChisq` <-
function(p, coefs = c(0.5, 0.5), df = 1, lower.tail = TRUE, nsims = 10000, ncp = 0, ...) {
    nchis <- length(coefs)
    if (lcombChisqParsOK(nchis, df, ncp))
        return( quantile( as.numeric( coefs %*% replicate(nsims, rchisq(n=nchis, df=df, ncp=ncp))), probs = if (lower.tail) p else 1-p, ...))
}

