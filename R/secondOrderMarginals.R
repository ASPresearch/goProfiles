`secondOrderMarginals` <-
function(pGO) {
    namesVec <- names(pGO)
    setNamesList <- strsplit(namesVec,"\\.")
    pi.labels <- unique(unlist(setNamesList))
    len.pi. <- length(pi.labels)
    present.in <- matrix(unlist(lapply(as.numeric(pi.labels), belongs.to, setNamesList)), nrow=length(pi.labels), byrow=T)
    # possible bug in R 2.1.0, this does not work:
    #   outer(1:len.pi., 1:len.pi., FUN=sum.and, present.in, pGO)
    # instead use:
# -----------------------------------------------------------------------
    pij. <- matrix(NA, nrow=len.pi., ncol=len.pi.)
    colnames(pij.) <- pi.labels
    rownames(pij.) <- pi.labels
    diag(pij.) <- sapply(1:len.pi., sum.if, present.in, pGO)
    for (i in 1:len.pi.)
        for (j in 1:(i-1))
            pij.[i,j] <- sum.and(i, j, present.in, pGO)
    pij.[upper.tri(pij.,diag=F)] <- t(pij.)[upper.tri(pij.,diag=F)]
    return(pij.)
}

