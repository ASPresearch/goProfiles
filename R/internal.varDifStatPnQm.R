`internal.varDifStatPnQm` <-
function(pn, qm=NULL, pqn0=NULL) {
    fullNams <- unique(c(names(pn),names(qm),names(pqn0)))
    pn <- fullGOProfile(pn, fullNams)
    qm <- fullGOProfile(qm, fullNams)
    pqn0 <- fullGOProfile(pqn0, fullNams)
    n <- attr(pn,"ngenes")
    if (is.null(qm)) {
        # variance of the distance between a profile and a subset of it

        n0 <- attr(pqn0,"ngenes")
        n1 <- n - n0
        sigma <- sigma.restrict(p1=(n*pn-n0*pqn0)/n1, pq=pqn0, n1=n1, n0=n0)
        difPQ <- contractedProfile(pn)[,3]/n - contractedProfile(pqn0)[,3]/n0
    } else {
      if (is.null(pqn0)) {
        # variance of the distance between two disjoint profiles (no genes in common)
        m <- attr(qm,"ngenes")
        sigma <- sigma.excluding(p1=pn, q1=qm, n1=n, m1=m)
        difPQ <- contractedProfile(pn)[,3]/n - contractedProfile(qm)[,3]/m
      } else {
        # variance of the distance between two two intersecting profiles
        # (some genes are specific of pn, some are specific of qm, and some common genes are profiled in pqn0)
        m <- attr(qm,"ngenes")
        n0 <- attr(pqn0,"ngenes")
        n1 <- n - n0
        m1 <- m - n0
        sigma <- sigma.intersect(p1=(n*pn-n0*pqn0)/n1, q1=(m*qm-n0*pqn0)/m1, pq=pqn0, n1=n1, m1=m1, n0=n0)
        difPQ <- contractedProfile(pn)[,3]/n - contractedProfile(qm)[,3]/m
      }
    }
    difPQ <- as.matrix(difPQ)
    as.vector(4 * t(rbind(difPQ,-difPQ)) %*% sigma %*% rbind(difPQ,-difPQ))
}

