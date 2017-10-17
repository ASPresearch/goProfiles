`as.GOTerms.frame` <-
function(myGOTermsList, na.rm=TRUE){
  unlisted<-unlist(myGOTermsList)
  if (na.rm) 
    unlisted<-unlisted[!is.na(unlisted)]
  else
    names(unlisted[is.na(unlisted)])<- 
      sapply(names(unlisted[is.na(unlisted)]),function (ll) paste(ll,"NA","NA",sep="."))
  m <- matrix(unlist(strsplit(names(unlisted),"\\.")), ncol=2,byrow=T)
  m21<-sapply (m[,2],function (s) substr(s,1,2))
  m22<-sapply (m[,2],function (s) substr(s,4,7))
  GOTermsFrame<-data.frame(m[,1], m21, m22, unlisted)
  names(GOTermsFrame)<-c("GeneID","Ontology","Evidence","GOID")
  return(GOTermsFrame)
}

