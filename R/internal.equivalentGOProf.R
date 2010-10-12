
internal.equivalentGOProf <- function(pn, qm, pqn0, n, m, n0, confidence, d0, equivEpsilon)
{
    dataNames <- paste(deparse(substitute(pn)), " and ", deparse(substitute(qm)), " and ", deparse(substitute(pqn0)), sep="")
    fullNams <- unique(c(names(pn),names(qm),names(pqn0)))
    pn <- fullGOProfile(pn, fullNams)
    qm <- fullGOProfile(qm, fullNams)
    pqn0 <- fullGOProfile(pqn0, fullNams)
    contrPn <- contractedProfile.default(pn)

    if (is.null(qm)) {
      # equivalence of a profile with a subset of it
      d <- dEuclid2(contrPn, contrQm <- contractedProfile.default(pqn0))
      m <- n0
    }
    else if (is.null(pqn0)) {
      # equivalence of two disjoint profiles (no gens in common)
      d <- dEuclid2(contrPn, contrQm <- contractedProfile.default(qm))
    }
    else {
      # equivalence of two intersecting profiles (some genes are specific of pn, some are specific of qm, and some common genes are profiled in pqn0)
      d <- dEuclid2(contrPn, contrQm <- contractedProfile.default(qm))
    }
    s <- sum((contrPn[,3] + contrQm[,3]) > 0) # classes being compared in the profiles
    if (is.null(d0))
      d0 <- s * equivEpsilon * equivEpsilon

    names(d) <- "sample squared Euclidean distance"
    se <- sqrt(internal.varDifStatPnQm(pn=pn, qm=qm, pqn0 = pqn0) * ((n+m)/(n*m)))
    names(se) <- "distance standard error"
    attr(d,"se") <- se
    icDistance.oneSided <- c(0, d.upper <- d + qnorm(confidence) * se)
    names(icDistance.oneSided) <- c("origin", "one-sided upper")
    stat <- (d - d0) / se
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
      data.name = dataNames,
      alternative = paste("Equivalence or similarity, true squared Euclidean distance between the contracted profiles is less than ",
      d0, sep="")
  )
  class(result) <- "htest"
  return(result)

}
