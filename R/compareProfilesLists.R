`compareProfilesLists` <-
function (expanded1, expanded2=NULL, common.expanded=NULL, relationType,
	  method = "lcombChisq", ab.approx = "asymptotic", confidence = 0.95,...)
{ compared <-list() 
  for(i in 1:length(expanded1)){
    if (relationType=="DISJOINT"){
        compared[[i]]<-compareGOProfiles (pn=expanded1[[i]], qm=expanded2[[i]],
                          method=method, ab.aprox=ab.approx, confidence=confidence,...)
      }else{
        if (relationType=="INCLUSION"){
              compared[[i]]<-compareGOProfiles (pn=expanded1[[i]], 
                              pqn0=common.expanded[[i]],
                              method=method, ab.aprox=ab.approx, confidence=confidence,...)
        }else{
                compared[[i]]<-compareGOProfiles (pn=expanded1[[i]], 
                          qm=expanded2[[i]], pqn0=common.expanded[[i]],
                          method=method, ab.aprox=ab.approx, confidence=confidence,...)}
                }
        names(compared)[i]<-attr(expanded1[[i]],"onto")
        }
return(compared)
}  

