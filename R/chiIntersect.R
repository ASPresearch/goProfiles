chiIntersect <- function(d2, pCommon, n=attr(pCommon,"ngenes"), m, pq, n0=attr(pq,"ngenes"),
    lambda = n/(n+m), coefStat = m*lambda,
    sigma0=internal.covGO(pq), betas=NULL,
    dfr = NULL,
    ab.approx = "asymptotic",
    nsims=10000)
{
    if (is.null(dfr))
        dfr <- sum(contractedProfile(pCommon) > 0)
    if (pmatch(ab.approx,"asymptotic",nomatch=F)) {
        if (is.null(betas)) {
            n1 <- n - n0
            m1 <- m - n0
            p1 <- (n/n1)*pCommon - (n0/n1)*pq
            p1[p1 < 0] <- 0
            p1 <- p1 / sum(p1)
            q1 <- (m/m1)*pCommon - (n0/m1)*pq
            q1[q1 < 0] <- 0
            q1 <- q1 / sum(q1)
            sigma.p <- internal.covGO(p1)
            sigma.q <- internal.covGO(q1)
            k <- ncol(sigma0)
            theta <- n0 / (n + m)
            sigma <- cbind(rbind(diag(k),-diag(k)),rbind(-diag(k),diag(k))) %*%
              cbind(rbind((1-lambda)*((n0/n)*sigma0+(n1/n)*sigma.p), theta*sigma0),
                rbind(theta*sigma0, lambda*((n0/m)*sigma0+(m1/m)*sigma.q)))
            mean.stat <- sum(diag(sigma))
            sigma <- sigma * sigma
            var.stat <-  2 * sum(sigma)
        }
      else {
            mean.stat <- sum(betas)
            var.stat <- 2 * sum(betas*betas)
      }
    }
    else {
        n1 <- n - n0
        m1 <- m - n0
        p1 <- (n/n1)*pCommon - (n0/n1)*pq
        p1[p1 < 0] <- 0
        p1 <- p1 / sum(p1)
        q1 <- (m/m1)*pCommon - (n0/m1)*pq
        q1[q1 < 0] <- 0
        q1 <- q1 / sum(q1)
        mnomp <- multinomial(n=n1, p=p1)
        mnomq <- multinomial(n=m1, p=q1)
        mnomp0 <- multinomial(n=n0, p=pq)
        p0Sampl <- as.data.frame(t(generate.multinomial(mnomp0, ssize=nsims)))
        pnSampl <- as.data.frame(t(generate.multinomial(mnomp, ssize=nsims))) + p0Sampl
        pnSampl <- (sapply(pnSampl/n, contractCol, row.names(pnSampl)))
        qmSampl <- as.data.frame(t(generate.multinomial(mnomq, ssize=nsims))) + p0Sampl
        qmSampl <- (sapply(qmSampl/m, contractCol, row.names(qmSampl)))

        dsamp <- sapply(1:nsims, idEuclid2, pnSampl, qmSampl)
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
