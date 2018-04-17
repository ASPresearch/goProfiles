`ngenes.numeric` <-
function(pn, i=NULL) {
  if (has.ngenes.attr(pn))
    return(attr(pn,"ngenes"))
  else if (has.numGenes.attr(pn))
    return(attr(pn,"numGenes"))
  else
    return(NULL)
}

