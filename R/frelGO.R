`frelGO` <-
function(pGO) {
    n <- sum(pGO)
    pGO <- pGO / n
    if (!has.ngenes.attr(pGO))
        attr(pGO,"ngenes") <- n
    return(pGO)
}

