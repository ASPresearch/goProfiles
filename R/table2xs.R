
table2xs <- function(pn, ...) {
    UseMethod("table2xs")
}

table2xs.BasicGOProfile <- function(pn, qm=NULL, pqn0=NULL, goIds=F)
{
  freq.col <- 3
  namsPn <- rownames(pn)
  n <- ngenes(pn)
  if (!is.null(qm)) {
    namsQm <- rownames(qm)
    m <- ngenes(qm)
  }
  else namsQm <- NULL
  if (!is.null(pqn0)) {
      namsPQn0 <- rownames(pqn0)
      n0 <- ngenes(pqn0)
      if (!is.null(qm)) {
        if (!all((namsPQn0 %in% namsQm)))
          stop("All GO nodes in 'pqn0' should be also present in 'qm'")
      }
      else m <- n0
      if (!all((namsPQn0 %in% namsPn)))
        stop("All GO nodes in 'pqn0' should be also present in 'pn'")
      commn <- namsPn %in% namsPQn0
      pn[commn,freq.col] <- pn[commn,freq.col] - pqn0[namsPn[commn],freq.col]
      n <- n - n0
  }
  else namsPQn0 <- NULL
  nams <- unique(c(namsPn, namsQm, namsPQn0))
  tabFreqs <- matrix(0, nrow = length(nams), ncol=2)
  rownames(tabFreqs) <- nams
  tabFreqs[,1][aux] <- pn[nams[aux <- nams %in% namsPn],freq.col]
  if (!is.null(qm))
    tabFreqs[,2][aux] <- qm[nams[aux <- nams %in% namsQm],freq.col]
  else
    tabFreqs[,2][aux] <- pqn0[nams[aux <- nams %in% namsPQn0],freq.col]
  tabFreqs <- tabFreqs[tabFreqs[,1]|tabFreqs[,2],]
  attr(tabFreqs,"ngenes") <- c(n,m)
  if (goIds) {
    goIdsTable <- unique(c(as.character(pn[rownames(pn),"GOID"]), as.character(qm[rownames(qm),"GOID"])))
    names(goIdsTable) <- unique(c(rownames(pn), rownames(qm)))
    rownames(tabFreqs) <- goIdsTable[rownames(tabFreqs)]
  }
  return(tabFreqs)
}
