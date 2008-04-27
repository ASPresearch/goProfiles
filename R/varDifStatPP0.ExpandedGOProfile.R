`varDifStatPP0.ExpandedGOProfile` <-
function(pn, p0 = NULL, funcProfP0 = NULL, simplify=T) {
  pnName <- deparse(substitute(pn))
  vecP0 <- NULL
  vecFuncProfP0 <- NULL
  if (!is.null(funcProfP0)) {
    p0Name <- deparse(substitute(funcProfP0))
    # ATENCIÓ!!!! FALTA MANEGAR funcProfP0 COM OBJECTE DE CLASSE FunctionalGOProfile:
    # funcProfP0 <- as.FunctionalGOProfile(funcProfP0)
    givenFuncProfP0 <- T
  }
  else if (!is.null(p0)) {
    p0Name <- deparse(substitute(p0))
    p0 <- as.ExpandedGOProfile(p0)
    givenFuncProfP0 <- F
  }
  else stop("NULL reference or population profile")
  var.jk <- function(i) {
    j <- i %% ncolPn + 1
    vecPn <- pn[,j]
    names(vecPn) <- row.names(pn)
    if (givenFuncProfP0) {
      k <- i %% ncolFuncProfP0 + 1
      vecFuncProfP0 <- funcProfP0[,k]
      names(vecFuncProfP0) <- rownames(funcProfP0)
    }
    else {
      k <- i %% ncolP0 + 1
      vecP0 <- p0[,k]
      names(vecP0) <- row.names(p0)
    }
    result.jk <- internal.varDifStatPP0(vecPn, p0=vecP0, funcProfP0=vecFuncProfP0)
    if (!simplify)
      attr(result.jk,"data.name") <- paste(pnName,"[",j,"] and ",p0Name,"[",k,"]", sep="")
    result.jk
  }
  maxncol <- max(
    ncolPn <- ncol(pn),
    max(ncolP0 <- ncol(p0), ncolFuncProfP0 <- ncol(funcProfP0))
  )
  return(sapply(0:(maxncol-1), var.jk, simplify=simplify))
}

