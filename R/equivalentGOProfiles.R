equivalentGOProfiles <- function(goObject, ...) {
  UseMethod("equivalentGOProfiles")
}

equivalentGOProfiles.list <- function(goObject, ...) {
  sapply(goObject, equivalentGOProfiles, ...)
}

equivalentGOProfiles.GOProfileHtest <- function(goObject, equivEpsilon = 0.05, d0 = NULL, confidence = NULL, ...)
{
  if (is.null(d0))
    d0 <- sum((goObject$profilePn > 0) | (goObject$profileQm > 0)) * equivEpsilon * equivEpsilon
  if (is.null(confidence))
    confidence <- attr(goObject$conf.int,"conf.level")
  d <- goObject$estimate
  se <- attr(goObject$estimate,"se")
  icDistance.oneSided <- c(0, d.upper <- d + qnorm(confidence) * se)
  names(icDistance.oneSided) <- c("origin", "one-sided upper")

  n <- attr(goObject$profilePn,"ngenes")
  m <- attr(goObject$profileQm,"ngenes")
  stat <- (goObject$estimate - d0) / se
  names(stat) <- "(d - d0) / se(d)"
  attr(stat,"se") <- NULL
  pval <- pnorm(stat)
  names(pval) <- NULL
  params <- c(d0, n, m)
  names(params) <- c("d0", "n", "m")  
  result <- list(
    statistic = stat,
    parameter = params,
    p.value = pval,
    conf.int = icDistance.oneSided,
    estimate = d,
    data.name = deparse(substitute(goObject)),
    alternative = paste("Equivalence or similarity, true squared Euclidean distance between the contracted profiles is less than ",
      d0, sep="")
  )
  class(result) <- "htest"
  return(result)
}

equivalentGOProfiles.ExpandedGOProfile <- function(goObject, qm=NULL, pqn0=NULL,
    n = ngenes(goObject), m = ngenes(qm), n0 = ngenes(pqn0),
    confidence = 0.95,
    equivEpsilon = 0.05, d0 = NULL, 
    simplify = FALSE, ...)
{
  pnName <- deparse(substitute(goObject))
  qmName <- deparse(substitute(qm))
  pqn0Name <- deparse(substitute(pqn0))
  if (!is.null(qm)) qm <- as.ExpandedGOProfile(qm)
  if (!is.null(pqn0)) pqn0 <- as.ExpandedGOProfile(pqn0)
  test.jkl <- function(i) {
    j <- i %% ncolPn + 1
    vecPn <- goObject[,j]
    names(vecPn) <- rownames(goObject)
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
    result.jkl <- internal.equivalentGOProf(vecPn, vecQm, vecPQn0, n[j], m[k], n0[l],
        confidence, d0, equivEpsilon)
    result.jkl$data.name <-
        paste(pnName,"[",j,"] and ", qmName.k, " and ", pqn0Name.l, sep="")
    result.jkl
  }
  maxncol <- max(ncolPn <- ncol(goObject),ncolQm <- ncol(qm),ncolPQn0 <- ncol(pqn0))
  result <- lapply(0:(maxncol-1), test.jkl)
  if (simplify && (maxncol == 1)) {
    result[[1]]$data.name <- paste(pnName," and ",qmName," and ",pqn0Name, sep="")
    return(result[[1]])
  }
  else
    return(result)
}

equivalentGOProfiles.default <- function(goObject, ...)  {
  equivalentGOProfiles( as.ExpandedGOProfile(goObject), ...)
}
