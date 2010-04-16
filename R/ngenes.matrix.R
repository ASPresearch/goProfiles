`ngenes.matrix` <-
 function(pn, i=NULL) {
  if (is.null(i)) {
    if (has.ngenes.attr(pn))
      return(attr(pn,"ngenes"))
    else if (has.numGenes.attr(pn))
      return(attr(pn,"numGenes"))
    else
      return(NULL)
  }
  else
    return(ngenes(pn)[i])
}

