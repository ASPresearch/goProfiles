`contractedProfile` <-
function(pGO) {
    namesVec <- names(pGO)
    setNamesList <- strsplit(namesVec,"\\.")
    pi.labels <- unique(unlist(setNamesList))
    len.pi. <- length(pi.labels)
    present.in <- matrix(unlist(lapply(as.numeric(pi.labels), belongs.to, setNamesList)), nrow=len.pi., byrow=T)
    pi. <- sapply(1:len.pi., sum.if, present.in, pGO)
    names(pi.) <- pi.labels
    return(pi.)
}

