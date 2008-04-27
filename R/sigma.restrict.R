`sigma.restrict` <-
function(p1, pq, n1 = attr(p1,"ngenes"), n0 = attr(pq,"ngenes")) {
    sigma0 <- internal.covGO(pq)
    sigma.p <- internal.covGO(p1)
    n <- n0 + n1
    lambda <- n / (n + n0)
    theta <- 1 - lambda
    cbind(rbind(theta*((n0/n)*sigma0 + (n1/n)*sigma.p), theta*sigma0),
          rbind(theta*sigma0, lambda*sigma0))
}

