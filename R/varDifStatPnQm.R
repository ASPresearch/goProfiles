`varDifStatPnQm` <-
function(pn, qm=NULL, pqn0=NULL, simplify=T)
{
  pnName <- deparse(substitute(pn))
  qmName <- deparse(substitute(qm))
  pqn0Name <- deparse(substitute(pqn0))
  pn <- as.ExpandedGOProfile(pn)
  if (!is.null(qm)) qm <- as.ExpandedGOProfile(qm)
  if (!is.null(pqn0)) pqn0 <- as.ExpandedGOProfile(pqn0)
  var.jkl <- function(i) {
    j <- i %% ncolPn + 1
    vecPn <- pn[,j]
    names(vecPn) <- rownames(pn)
    if (is.null(qm)) {
      vecQm <- NULL
      qmName.k <- NULL
    }
    else {
      k <- i %% ncolQm + 1
      vecQm <- qm[,k]
      names(vecQm) <- rownames(qm)
      qmName.k <- paste(qmName,"[",k,"]", sep="")
    }
    if (is.null(pqn0)) {
      vecPQn0 <- NULL
      pqn0Name.l <- NULL
    }
    else {
      l <- i %% ncolPQn0 + 1
      vecPQn0 <- pqn0[,l]
      names(vecPQn0) <- rownames(pqn0)
      pqn0Name.l <- paste(pqn0Name,"[",l,"]", sep="")
    }
    result.jkl <- internal.varDifStatPnQm(vecPn, vecQm, vecPQn0)
    if (!simplify)
      result.jkl$data.name <- 
      paste(pnName,"[",j,"] and ", qmName.k, " and ", pqn0Name.l, sep="")
    result.jkl
  }
  maxncol <- max(ncolPn <- ncol(pn),ncolQm <- ncol(qm),ncolPQn0 <- ncol(pqn0))
  return(sapply(0:(maxncol-1), var.jkl, simplify=simplify))
}

