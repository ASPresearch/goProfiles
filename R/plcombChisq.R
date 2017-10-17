`plcombChisq` <-
function(q, coefs = c(0.5, 0.5), df = 1, lower.tail = TRUE, nsims = 10000, ncp = 0, ...) {
    nchis <- length(coefs)
    if (lcombChisqParsOK(nchis, df, ncp)) {
        p <- sum( as.numeric( coefs %*% replicate(nsims, rchisq(n=nchis, df=df, ncp=ncp))) <= q) / nsims
        if (lower.tail)
            return(p)
        else
            return(1-p)
    }
}

