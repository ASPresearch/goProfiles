tolWholeNumber <- .Machine$double.eps^0.5
all.wholenumber <- function(x, tol = tolWholeNumber)  all(abs(x - round(x)) < tol)
is.wholenumber <- function(x, tol = tolWholeNumber)  abs(x - round(x)) < tol

fisherGOProfiles <- function(pn, ...) {
    UseMethod("fisherGOProfiles")
}

fisherGOProfiles.numeric <- function(pn, qm=NULL, pqn0=NULL,
    n = ngenes(pn), m = ngenes(qm), n0 = ngenes(pqn0),
    method = "BH", simplify=T, expanded=F, ...)
{
  if (expanded) {
    attr(pn,"ngenes") <- n
    if (!is.null(qm)) attr(qm,"ngenes") <- m
    if (!is.null(pqn0)) attr(pqn0,"ngenes") <- n0
    fisherGOProfiles.ExpandedGOProfile(
      as.ExpandedGOProfile(pn), qm, pqn0,
      method = method, simplify = simplify, ...
    )
  }
  else {
    namsPn <- names(pn)
    if (!all.wholenumber(pn)) pn <- pn * n
    if (!is.null(qm)) {
      namsQm <- names(qm)
      if (!all.wholenumber(qm)) qm <- qm * m
    }
    else namsQm <- NULL
    if (!is.null(pqn0)) {
      namsPQn0 <- names(pqn0)
      if (!all.wholenumber(pqn0)) pqn0 <- pqn0 * n0
      if (!is.null(qm)) {
        if (!all((namsPQn0 %in% namsQm)))
          stop("All GO nodes in 'pqn0' should be also present in 'qm'")
      }
      else m <- n0
      if (!all((namsPQn0 %in% namsPn)))
        stop("All GO nodes in 'pqn0' should be also present in 'pn'")
      commn <- namsPn %in% namsPQn0
      pn[commn] <- pn[commn] - pqn0[namsPn[commn]]
      n <- n - n0
    }
    else namsPQn0 <- NULL
    nams <- unique(c(namsPn, namsQm, namsPQn0))
    tabFreqs <- matrix(0, nrow = length(nams), ncol=2)
    rownames(tabFreqs) <- nams
    tabFreqs[,1][aux] <- pn[nams[aux <- nams %in% namsPn]]
    if (!is.null(qm))
      tabFreqs[,2][aux] <- qm[nams[aux <- nams %in% namsQm]]
    else
      tabFreqs[,2][aux] <- pqn0[nams[aux <- nams %in% namsPQn0]]
    fisherGOProfiles.matrix(tabFreqs, n, m, method, ...)
  }
}

fisherGOProfiles.matrix <- function(pn, n, m, method = "BH", ...)
{
  if (ncol(pn) != 2) {
    pn <- t(pn)
    if (ncol(pn) != 2) stop("Argument 'pn' should have 2 columns in 'fisherGOProfiles.matrix'")
  }
  if (!all.wholenumber(pn)) stop("Elements of argument 'pn' should be all whole numbers")
  if (missing(m)) {
    if (missing(n)) {
      if (has.ngenes.attr(pn)){
        n <- ngenes(pn, 1)
        m <- ngenes(pn, 2)
      }
      else stop("Argument 'pn' should have a 'ngenes' attribute (a numeric(2)) or arguments 'n', 'm' should be specied")
    }
    else {
      m <- n[2]
      n <- n[1]
    }
  }
  else if (missing(n)) stop("Argument 'n' should be specified if 'm' is given")
  internal.enrichProfile(pn[,1], pn[,2], n, m, method, ...)
}

fisherGOProfiles.BasicGOProfile <- function(pn, qm=NULL, pqn0=NULL,
    method = "BH", goIds=T, ...)
{
  fisherGOProfiles.matrix(table2xs(pn, qm, pqn0, goIds=goIds), method=method)
}

fisherGOProfiles.ExpandedGOProfile <- function(pn, qm=NULL, pqn0=NULL,
    method = "BH", simplify=T, ...)
{
  pnName <- deparse(substitute(pn))
  qmName <- deparse(substitute(qm))
  pqn0Name <- deparse(substitute(pqn0))
  #pn <- as.ExpandedGOProfile(pn)
  n = ngenes(pn)
  m = ngenes(qm)
  n0 = ngenes(pqn0)
  if (!is.null(qm)) qm <- as.ExpandedGOProfile(qm)
  if (!is.null(pqn0)) pqn0 <- as.ExpandedGOProfile(pqn0)
  test.jkl <- function(i) {
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
    result.jkl <- internal.enrich(vecPn, vecQm, vecPQn0, n[j], m[k], n0[l],
        method, ...)
    attr(result.jkl, "data.name") <-
        paste(pnName,"[",j,"] and ", qmName.k, " and ", pqn0Name.l, sep="")
    result.jkl
  }
  maxncol <- max(ncolPn <- ncol(pn),ncolQm <- ncol(qm),ncolPQn0 <- ncol(pqn0))
  result <- lapply(0:(maxncol-1), test.jkl)
  if (simplify && (maxncol == 1)) {
    attr(result[[1]], "data.name") <- paste(pnName," and ",qmName," and ",pqn0Name, sep="")
    return(result[[1]])
  }
  else
    return(result)
}
