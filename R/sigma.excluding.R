`sigma.excluding` <-
function(p1, q1, n1 = attr(p1,"ngenes"), m1 = attr(q1,"ngenes")) {
    sigma.p <- internal.covGO(p1)
    sigma.q <- internal.covGO(q1)
    zeros <- matrix(0, ncol=ncol(sigma.p), nrow=nrow(sigma.p))
    lambda <- n1 / (n1 + m1)
    cbind(rbind((1-lambda)*sigma.p, zeros), rbind(zeros, lambda*sigma.q))
}

