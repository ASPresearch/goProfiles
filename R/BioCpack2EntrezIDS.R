`BioCpack2EntrezIDS` <-
function (anotPkg, na.rm=FALSE)
{
  stopifnot(require(anotPkg, character.only=TRUE, quietly=TRUE, warn.conflicts=FALSE))
  envENTREZID<-paste(anotPkg, "ENTREZID",sep="")
  envENTREZID <- eval(parse(text = envENTREZID))
  gLL<-as.list(envENTREZID)
  gLL<-as.character(unlist(gLL))
  if (na.rm)
    gLL<-gLL[!is.na(gLL)]
  return(gLL)
}

