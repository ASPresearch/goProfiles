
# --------------------------------------------------------------------------------------------------
#        Clustering according to equivalence threshold distance (itself based on dEuclid^2) 
# --------------------------------------------------------------------------------------------------

# ..................................................................................................
#           Building all expanded profiles, for all lists and intersections of lists
# ..................................................................................................
intersectOne.vs.other <- function(iList2, iList, geneLists, onto, ontoLevel, #orgPackage,
                         trace, ...) {
  if (trace) {
    nams <- names(geneLists)
    max.length <- max(str_length(nams))
    cat(str_pad(nams[iList], max.length, "right"), ",",
        str_pad(nams[iList2], max.length, "right"), "|", sep = "")
  }
  listsIntersect <- intersect(geneLists[[iList]], geneLists[[iList2]])
  if (length(listsIntersect) > 0) {
    return(expandedProfile(listsIntersect, onto = onto, level = ontoLevel, ...))
  } else {
    return(NULL)
  }
}

intersectOne.vs.remaining <- function(iList, geneLists, onto, ontoLevel, #orgPackage, 
                             trace, ...) {
  if (trace) {
    cat("\n")
  }
  result <- lapply(1:(iList-1), intersectOne.vs.other,
                   iList = iList, geneLists = geneLists,
                   onto = onto, ontoLevel = ontoLevel, 
                   trace = trace, ...)
  names(result) <- names(geneLists)[1:(iList-1)]
  return(result)
}

allIntersections <- function(geneLists, onto, ontoLevel, 
                           trace = TRUE, ...) {
  lenList <- length(geneLists)
  result <- lapply(2:lenList, intersectOne.vs.remaining, geneLists = geneLists,
                   onto = onto, ontoLevel = ontoLevel, 
                   trace = trace, ...)
  names(result) <- names(geneLists)[-1]
  return(result)
}

marginalProfile <- function(iList, geneLists, onto, ontoLevel, trace, ...) {
  if (trace) {
    cat( "Building profile for list ", names(geneLists)[iList], "\n")
  }
  return(expandedProfile(geneLists[[iList]], onto = onto, level = ontoLevel, ...))
}

# ..................................................................................................
#                     Equivalence test between all lists
# ..................................................................................................
one.vs.other <- function(iList2, iList, 
                         allMarginalProfiles, allIntersectionProfiles, 
                         onto, ontoLevel, #orgPackage,
                         trace, ...) {
  if (trace) {
    # cat(".")
    nams <- names(allMarginalProfiles)
    max.length <- max(str_length(nams))
    cat(str_pad(nams[iList], max.length, "right"), ",",
        str_pad(nams[iList2], max.length, "right"), "|", sep = "")
  }
  getStats(equivalentGOProfiles(allMarginalProfiles[[iList]][[1]], allMarginalProfiles[[iList2]][[1]],
                                allIntersectionProfiles[[iList-1]][[iList2]][[1]],
                            # onto = onto, level = ontoLevel, 
                            # orgPackage = orgPackage, 
                            ...))
}

one.vs.remaining <- function(iList, allMarginalProfiles, allIntersectionProfiles, 
                             onto, ontoLevel, #orgPackage, 
                             trace, ...) {
  if (trace) {
    cat("\n")
  }
  result <- lapply(1:(iList-1), one.vs.other,
                   iList = iList, 
                   allMarginalProfiles = allMarginalProfiles, 
                   allIntersectionProfiles = allIntersectionProfiles,
                   onto = onto, ontoLevel = ontoLevel, 
                   trace = trace, ...)
  names(result) <- names(allMarginalProfiles)[1:(iList-1)]
  return(result)
}

allEquivTests <- function(allMarginalProfiles, allIntersectionProfiles, 
                          onto, ontoLevel, 
                          trace = TRUE, ...) {
  lenList <- length(allMarginalProfiles)
  result <- lapply(2:lenList, one.vs.remaining, 
                   allMarginalProfiles = allMarginalProfiles, 
                   allIntersectionProfiles = allIntersectionProfiles,
                   onto = onto, ontoLevel = ontoLevel, 
                   trace = trace, ...)
  names(result) <- names(allMarginalProfiles)[-1]
  return(result)
}

# ..................................................................................................
#     Euclidean squared distance and its standard error from 'equivalentGOProfiles ' statistic
# ..................................................................................................
getStats <- function(goObject)
{
  d <- goObject[[1]]$estimate
  result <- c(d, attr(d,"se"))
  names(result) <- c("d", "se")
  attr(result, "pn") <- goObject[[1]]$profilePn
  attr(result, "qm") <- goObject[[1]]$profileQm
  return(result)
}


#' For a given level (2, 3, ...) in a GO ontology (BP, MF or CC), compute the equivalence threshold 
#'   distance matrix and generate a dendrogram from it.
#'
#' @param ontoLevel integer (2, 3, ...) level of a GO ontology where the GO profiles are built
#' @param onto character, GO ontology ("BP", "MF" or "CC") under consideration
#' @param geneLists list of character vectors, each vector stands for the gene names in a given gene set
#' @param trace boolean, the full process must be traced? Defaults to TRUE
#' @param onTheFlyDev character, name of the graphical device where to immediately display the resulting
#'   diagram. The appropriate names depend on the operating system. Defaults to \code{NULL} and then
#'   nothing is displayed
#' @param method character, one of the admissible methods in function \code{hclust}. Defaults to "complete"
#' @param jobName character, main plot name, defaults to 
#'   \code{paste("Equivalence cluster", onto, ontoLevel, method, sep = "_")}
#' @param ylab character, label of the vertical axis of the plot, defaults to "Equivalence threshold distance"
#' @param alpha simultaneous nominal significance level for the equivalence tests to be repeteadly performed,
#'   defaults to 0.05
#' @param precis numerical precission in the iterative search of the equivalence threshold distances,
#'   defaults to 0.001
#' @param ... additional arguments to \code{hclust}
#' @return An object of class \code{equivClust}, descending from class \code{hclust}
#'   with some additional attributes:
#'   \describe{
#'     \item{jobName}{The main job name}
#'     \item{sub}{The graphic subtittle}
#'     \item{ylab}{The vertical axis label}
#'     \item{distMat}{The equivalence threshold distance matrix}
#'     \item{allComps}{A list with some information on all the pairwise equivalence tests:
#'     the Euclidean squared distance, its standard error and the corresponding GO profiles}
#'   }
#' @details Do not confuse the threshold distance matrix with the squared distances computed
#' in each equivalence test.
#' @importFrom stringr str_pad
#' @examples
#' \dontrun{
#' data(kidneyGeneLists)
#' clustMF2 <- equivClust(2, "MF", kidneyGeneLists, orgPackage="org.Hs.eg.db")
#' plot(clustMF2)
#' plot(clustMF2, 
#'      main = "Dendrogram (method = complete)", sub = attr(clustMF2, "sub"),
#'      ylab = "Equivalence threshold distance")
#' # With the same data, an UPGMA dendrogram:
#' equivClust(2, "MF", kidneyGeneLists, method = "average",
#' orgPackage="org.Hs.eg.db")
#' }
#' @export
equivClust <- function(ontoLevel, onto, geneLists, 
                               trace = TRUE, onTheFlyDev = NULL, method = "complete", 
                               jobName = paste("Equivalence cluster", onto, ontoLevel, method, sep = "_"), 
                               ylab = "Equivalence threshold distance", 
                               alpha = 0.05, precis = 0.001, ...)
{
  subName <- paste0("Ontology ", onto, " at level ", ontoLevel, sep = "")
  if (trace) {
    cat("\n\n", jobName, subName, "\n")
  }
  s <- length(geneLists)
  h <- s * (s - 1) / 2
  equivDists <- rep(NA, h)
  dIdxs <- unlist(lapply(2:s, function(i) lapply(1:(i-1), function(j) c(i,j))), recursive = FALSE)
  if (trace) {
    cat("\nBuilding marginal profiles:\n\n")
  }
  allMarginalProfiles <- lapply(1:length(geneLists), marginalProfile, 
                                geneLists = geneLists,
                                onto = onto, ontoLevel = ontoLevel, 
                                trace = trace, ...)
  names(allMarginalProfiles) <- names(geneLists)
  if (trace) {
    cat("\nBuilding intersection profiles:\n")
  }
  allIntersectionProfiles <- allIntersections(geneLists, onto = onto, ontoLevel = ontoLevel,
                                              trace = trace, ...)
  if (trace) {
    cat("\n\nPerforming all equivalence tests:\n")
  }
  allComps <- allEquivTests(allMarginalProfiles, allIntersectionProfiles, 
                            onto = onto, ontoLevel = ontoLevel,
                            trace = trace, ...)
  d.se <- matrix(unlist(sapply(allComps, function(dse) unlist(dse))), ncol = h)
  delta <- max(c(1, qnorm(1 - alpha / h)) %*% d.se)
  a <- d.se[1,] / d.se[2,]
  names(a) <- 1:h
  b <- -1 / d.se[2,]
  
  for (ih in h:1) {
    nextStep <- nextEquivDist(a, b, ih, delta, alpha, precis)
    delta <- nextStep$delta
    iDelta <- nextStep$iDelta
    equivDists[iDelta] <- delta
    a[iDelta] <- NA
  }
  distMat <- matrix(0, nrow = s, ncol = s)
  rownames(distMat) <- colnames(distMat) <- names(geneLists)
  
  distMat[upper.tri(distMat)] <- equivDists
  distMat <- as.dist(t(distMat))
  attr(distMat, "dist.method") <- "equivalence"
  
  clust <- hclust(distMat, method = method)
  attr(clust, "jobName") <- jobName
  attr(clust, "ylab") <- ylab
  attr(clust, "sub") <- subName
  attr(clust, "distMat") <- distMat
  attr(clust, "allComps") <- allComps
  
  if (!is.null(onTheFlyDev)) {
    eval(call(onTheFlyDev, width = 20, height = 20))
    # dev.set(dev.list()[onTheFlyDev])
    plot(clust, hang = -1, main = jobName, sub = subName, ylab = ylab)
  }
  
  class(clust) <- c("equivClust", class(clust))
  return(clust)
}

nextEquivDist <- function(a, b, h, delta, alpha, precis) {
  incDelta <- 0.1 * delta
  alphas.holm <- alpha / seq(from = h, to = 1, by = -1)
  repeat {
    p.vals <- ifelse(is.na(a), NA, pnorm(a + b * delta)) 
    p.order <- order(p.vals)
    p.sort <- p.vals[p.order][1:h]
    # if (all(ifelse(is.na(p.sort), NA, p.sort <= alphas.holm), na.rm = TRUE)) {
    if (all(p.sort <= alphas.holm)) {
      if (incDelta < precis) {
        return(list(delta = delta, iDelta = p.order[h]))
      } else {
        delta <- delta - incDelta
      }
    } else {
      incDelta <- 0.5 * incDelta
      delta <- delta + incDelta
    }
  }
}

auxIter <- function(onto, ontoLevels, geneLists, 
                    trace, onTheFlyDev, method, jobName, ylab, alpha, precis,
                    ...) {
  result <- lapply(ontoLevels, equivClust, 
                   onto = onto, geneLists = geneLists,
                   trace = trace, onTheFlyDev = onTheFlyDev, method = method,
                   jobName = jobName, ylab = ylab, alpha = alpha, precis = precis,
                   ...)
  names(result) <- paste("Level", ontoLevels, sep = "")
  return(result)
}

#' For each combination of the specified levels in the choosen GO ontologies, 
#' compute the equivalence threshold distance matrix and generate a dendrogram from it.
#'
#' @param geneLists list of character vectors, each vector stands for the gene names in a given gene set
#' @param ontos character vector, (e.g. c("BP","MF")) indicating the GO ontologies to be analysed
#' @param ontoLevels integer vector (e.g. 2:4) indicating the GO levels in these ontologies 
#'   where the GO profiles are built
#' @param trace boolean, the full process must be traced? Defaults to TRUE
#' @param onTheFlyDev character, name of the graphical device where to immediately display the resulting
#'   diagrams. The appropriate names depend on the operating system. Defaults to \code{NULL} and then
#'   nothing is displayed. Otherwise, successive graphical windows are opened and the successive diagrams
#'   are displayed in them
#' @param method character, one of the admissible methods in function \code{hclust}. Defaults to "complete"
#' @param jobName character, main plot name, defaults to "Equivalence clustering" 
#' @param ylab character, label of the vertical axis of the plot, defaults to "Equivalence threshold distance"
#' @param alpha simultaneous nominal significance level for the equivalence tests to be repeteadly performed,
#'   defaults to 0.05
#' @param precis numerical precission in the iterative search of the equivalence threshold distances,
#'   defaults to 0.001
#' @param ... additional arguments to \code{hclust}
#' 
# CAL MODIFICAR (AQUESTA ?S LA D'equivClust):

#' @return An object of class \code{iterEquivCluster}. It is a list of \code{length(ontos)},
#'   one element for each ontology under study. Each element of this list is itself a list of
#'   \code{length(ontoLevels)} with elements of class \code{equivClust}, standing for the cluster equivalence
#'   analysis performed for each ontology and level analysed
#' @examples
#' \dontrun{
#' data(kidneyGeneLists)
#' kidneyGeneLists
#' genListsClusters <- iterEquivClust(kidneyGeneLists, ontoLevels = 2:3,
#'                                    jobName = "Kidney Gene Lists_Equivalence Clustering (complete)",
#'                                    ylab = "Equivalence threshold distance",
#'                                    orgPackage="org.Hs.eg.db", method = "complete")
#' genListsClusters[["BP"]][["Level 3"]]
#' class(genListsClusters[["BP"]][["Level 3"]])
#' }

#' @export
iterEquivClust <- function(geneLists, ontos = c("BP","MF","CC"), ontoLevels = c(2,3), 
                           trace = TRUE, onTheFlyDev = NULL, 
                           method = "complete", 
                           jobName = "Equivalence clustering", ylab = "Equivalence threshold distance",
                           alpha = 0.05, precis = 0.001,
                           ...)
{
  # if (!is.null(onTheFlyDev)) {
  #   eval(call(onTheFlyDev, 20,20))
  # }
  result <- lapply(ontos, auxIter, ontoLevels = ontoLevels, geneLists = geneLists,
                   trace = trace, onTheFlyDev = onTheFlyDev, method = method, 
                   jobName = jobName, ylab = ylab, alpha = alpha, precis = precis,
                   ...)
  names(result) <- ontos
  class(result) <- c("iterEquivClust", class(result))
  attr(result, "jobName") <- jobName
  attr(result, "ylab") <- ylab
  return(result)
}

