`ngenes.i` <-
 function(vecPn) {
  if (has.ngenes.attr(vecPn))
    return(attr(vecPn,"ngenes"))
  else if (has.numGenes.attr(vecPn))
    return(attr(vecPn,"numGenes"))
  else
    return(NULL)
}

