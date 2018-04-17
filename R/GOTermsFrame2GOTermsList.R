`GOTermsFrame2GOTermsList` <-
function (myGOTermsFrame, evid=FALSE)
{
  x<-tapply(as.character(myGOTermsFrame$GOID),as.character(myGOTermsFrame$GeneID),"[")
  for (i in 1:length(x))
  {
      selected<-myGOTermsFrame$GeneID==names(x[i])
      names(x[[i]])<- as.character(myGOTermsFrame$Ontology[selected])
      if (evid) 
       names(x[[i]])<-paste(names(x[[i]]),as.character(myGOTermsFrame$Evidence[selected]),
          sep="-")
  }
  return(x)
}

