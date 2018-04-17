`mergeProfilesLists` <-
function (profilesList1, profilesList2, emptyCats=F, profNames=NULL){
### This function combines two lists of profiles
### Check 1: Are they of the same type?
if (length(profilesList1)!=length(profilesList2))
      {stop("ERROR!. Profiles Lists es must be of same length")
  }else{
      merged<-list()
      for(i in 1:length(profilesList1)){
        merged[[i]]<-mergeProfiles(profilesList1[[i]],profilesList2[[i]], emptyCats=emptyCats, profNames=profNames)
        names(merged)[i] <-names(profilesList1)[i]}
  }
  return(merged)
}

