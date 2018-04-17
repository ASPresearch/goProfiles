`ngenes.ExpandedGOProfile` <-
function(pn, i=NULL) {
  if (is.null(i))
    return(sapply(pn, ngenes.i))
  else
    return(ngenes.i(pn[,i]))
}

