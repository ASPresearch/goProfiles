`covGO.ExpandedGOProfile` <-
function(pn, simplify=T) {
  covMat.i <- function(i) {
    vecPn <- pn[,i]
    names(vecPn) <- rownames(pn)
    internal.covGO(vecPn)
  }
  ncolPn <- ncol(pn)
  result <- lapply(1:ncolPn, covMat.i)
  if (simplify && (ncolPn == 1))
    return(result[[1]])
  else {
    names(result) <- paste("covGO", deparse(substitute(pn)), 1:ncolPn, sep=".")
    return(result)
  }
}

