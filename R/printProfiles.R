`printProfiles` <-
function(aProf, aTitle="Functional Profile", anOnto=NULL, percentage=FALSE, Width=25, emptyCats=FALSE) {
  cat(aTitle, paste(rep("=", nchar(aTitle)), collapse=""), sep="\n")
  if (is.data.frame(aProf)){
      if(!is.null(anOnto)) print(paste(anOnto, "ontology"))
        printOne ( aProf=aProf, aTitle=aTitle, percentage=percentage,
                  Width=Width, emptyCats=emptyCats)
  }else{
        for (i in 1:length(aProf)){
          print(paste(names(aProf)[i], "ontology"))
          printOne ( aProf=aProf[[i]], aTitle=aTitle, percentage=percentage,
                  Width=Width, emptyCats=emptyCats) }
  }    
}

`printOne` <-
function (aProf, aTitle="Functional Profile", percentage=FALSE, Width=25, emptyCats=TRUE){
    aProf$Description <-shortName(as.character(aProf$Description), aWidth=Width)
    if (percentage){
        if (ncol(aProf)==3){                                                              
            numGenes<-attr(aProf,"numGenes")
            if (!is.null(numGenes)& (numGenes> 0))
                freq<-round((aProf[,3]/numGenes*100),1)
            aProf$Frequency <- freq
        }else{
            numGenes1<-attr(aProf,"numGenes1")
            if (!is.null(numGenes1)& (numGenes1>0))
                aProf[,3]<-round((aProf[,3]/numGenes1*100),1)
            numGenes2<-attr(aProf,"numGenes2")
            if (!is.null(numGenes2)& (numGenes2> 0))
                aProf[,4]<-round((aProf[,4]/numGenes2*100),1)
          }
        }
    if(!emptyCats){
          if (ncol(aProf)==3)
            nonNullRows <-  which(sapply(aProf[,3],sum)>0)
          else    
             nonNullRows <-  which(apply(aProf[,3:ncol(aProf)],1,sum)>0)
     aProf<-aProf[nonNullRows,]
     }
    print(aProf)
 }
