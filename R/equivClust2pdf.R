# --------------------------------------------------------------------------------------------------
# Pdf generation of clusters
# --------------------------------------------------------------------------------------------------

auxPdf <- function(x, jobName, ylab) {
  lapply(x, function(clust) {
    plot(clust, hang = -1, main = jobName, sub = attr(clust, "sub"), ylab = ylab)
  })
}

#' Save the graphical representation of objects of class \code{equivClust} or \code{iterEquivClust}
#' as pdf files.
#'
#' @param x an object of class \code{equivClust} or \code{iterEquivClust} 
#' @param ... additional arguments to function \code{pdf}

#' @return NULL

#' @examples
#' data(clustKidneyMF2)
#' equivClust2pdf(clustKidneyMF2)
#' # And then open file "Equivalence cluster_MF_2_complete.pdf"...
#' equivClust2pdf(clustKidneyMF2, 
#'                jobName = "Method 'complete' dendrogram for level 2 of GO ontology MF")
#' # And then open file "Method 'complete' dendrogram for level 2 of GO ontology MF.pdf"...

#' @export
equivClust2pdf <- function(x, ...) {
  UseMethod("equivClust2pdf")
}

#'
#' @param jobName character, main plot title and file name (it should be correct as a file name!)
#' @param ylab character, label of the plot vertical axis
#' @method equivClust2pdf equivClust
#' @describeIn equivClust2pdf \code{equivClust2pdf} method for class \code{equivClust}
#' @return NULL

#' @export
equivClust2pdf.equivClust <- function(x, jobName, ylab, ...) {
  on.exit({
    dev.off(dev.list()["pdf"])
  })
  if (missing(jobName))
    jobName <- attr(x, "jobName")
  if (missing(ylab))
    ylab <- attr(x, "ylab")
  
  pdf(file = paste0(jobName, ".pdf", sep = ""), ...)
  plot(x, hang = -1, main = jobName, sub = attr(x, "sub"), ylab = ylab)
  return()
}

#' @method equivClust2pdf iterEquivClust
#' @return NULL
#' @describeIn equivClust2pdf \code{equivClust2pdf} method for class \code{iterEquivClust}
#' @export
equivClust2pdf.iterEquivClust <- function(x, jobName, ylab, ...) {
  on.exit({
    dev.off(dev.list()["pdf"])
  })
  if (missing(jobName))
    jobName <- attr(x, "jobName")
  if (missing(ylab))
    ylab <- attr(x, "ylab")
  
  pdf(file = paste0(jobName, ".pdf", sep = ""), ...)
  lapply(x, auxPdf, jobName = jobName, ylab = ylab)
}