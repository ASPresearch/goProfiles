`rlcombChisq` <-
function(n, coefs = c(0.5, 0.5), df = 1, ncp = 0) {
    nchis <- length(coefs)
    if (lcombChisqParsOK(nchis, df, ncp))
        return( as.numeric( coefs %*% replicate(n, rchisq(n = nchis, df = df, ncp = ncp))))
}

