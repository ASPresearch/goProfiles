`sigma.intersect` <-
function(p1, q1, pq, n1 = attr(p1,"ngenes"), m1 = attr(q1,"ngenes"), n0 = attr(pq,"ngenes")) {
    sigma0 <- internal.covGO(pq)
    sigma.p <- internal.covGO(p1)
    sigma.q <- internal.covGO(q1)
    n <- n0 + n1
    m <- n0 + m1
    lambda <- n / (n + m)
    theta <- n0 / (n + m)
    cbind(rbind((1-lambda)*((n0/n)*sigma0 + (n1/n)*sigma.p), theta*sigma0),
          rbind(theta*sigma0, lambda*((n0/m)*sigma0 + (m1/m)*sigma.q)))

}

