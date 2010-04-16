`mergeProfiles` <-
function(profile1, profile2, emptyCats=F, profNames=NULL){
### This function combines two individual data-frame-like profiles
### Check 1: Are both profiles of the same type
  if ((is.data.frame(profile1) & (! (is.data.frame(profile2)))) |
   (is.data.frame(profile2) & (! (is.data.frame(profile1)))))
  {stop("ERROR!. Profile es must be of same type")}
### Check 2: Do the profiles contain the same information in rows?
  iguals<-TRUE
  for(i in 1:nrow(profile1)) 
    if (profile1$Description[i]!=profile2$Description[i]) 
      iguals<-FALSE
### If checks are OK
    if (!(iguals)){
          stop("ERROR: Profiles have different row names")
    }else{
        merged<-cbind(profile1,profile2$Frequency)
        attr(merged,"numGenes1")<-ngen1 <- attr(profile1,"numGenes")
        attr(merged,"numGenes2")<-ngen2 <- attr(profile2,"numGenes")
        attr(merged,"numGenes")<-c(ngen1, ngen2)
        attr(merged,"numNAs") <- c(attr(profile1,"numNAs"),attr(profile2,"numNAs")) 
        if (!emptyCats){
          sumFreqs<-profile1$Frequency+profile2$Frequency
          merged<-merged[sumFreqs!=0,]
        }
        if (is.null(profNames))
            names(merged)<-c(names(profile1)[1:2],"Frequency-1","Frequency-2")
        else
            names(merged)<-c(names(profile1)[1:2],profNames[1], profNames[2])
    }
  class(merged) <- class(profile1)
  return(merged)
}

