`ngenes.default` <-
function(pn, i=NULL) {
  if (is.null(pn)) return(NULL)
  if (is.null(i))
    return(sapply(pn, ngenes.i))
  else
    return(ngenes.i(pn[,i]))
}

