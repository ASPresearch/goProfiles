contractedProfile.default <- function(pGO, nams = NULL) {
    if (is.null(nams))
      nams <- names(pGO)
    else
      if (length(nams) != length(pGO))
        stop("'length(nams)' and 'length(pGO)' must be equal in method 'contractedProfile.default'")
    setNamesList <- strsplit(nams,"\\.")
    pi.labels <- unique(unlist(setNamesList))
    len.pi. <- length(pi.labels)
    present.in <- matrix(unlist(lapply(as.numeric(pi.labels), belongs.to, setNamesList)), nrow=len.pi., byrow=T)
    pi. <- sapply(1:len.pi., sum.if, present.in, pGO)
    names(pi.) <- pi.labels
    attr(pi.,"ngenes") <- attr(pGO,"ngenes")
    return(pi.)
}