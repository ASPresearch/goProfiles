`geometricProfile` <-
function(ncateg, simult=ncateg, n, theta=0.05)
{
# Generates an (artificial) profile according to a geometric model (e.g. for simulation pourposes)
    nams <- fullSetOfNames(ncateg, simult)
    len <- length(nams)
    pGO <- (1 - theta)^(0:(len-1))
    pGO <- pGO / sum(pGO)
    names(pGO) <- nams
    attr(pGO,"ngenes") <- n
    return(pGO)
}

