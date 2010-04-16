`ngenes.BasicGOProfile` <-
 function(pn, i=NULL) {
  if (is.null(i))
    return(attr(pn, "numGenes"))
  else
    return(attr(pn, "numGenes")[i])
}

