`BioCprobes2Entrez` <-
function (probeslist , anotPkg, na.rm=TRUE)
{
stopifnot(require(anotPkg, character.only = TRUE, quietly=TRUE, warn.conflicts=FALSE))
myenvir<-eval(parse(text = paste(anotPkg,"ENTREZID",sep="")))
myENTREZIDs <- as.character(unlist(mget(probeslist,envir=myenvir, ifnotfound=NA)))
if(na.rm){
  return(myENTREZIDs[!(is.na(myENTREZIDs))])
}else{
  return(myENTREZIDs)}
}

