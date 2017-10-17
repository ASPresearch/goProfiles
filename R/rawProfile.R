`rawProfile` <-
function (ancestorsLst, ontoLevelsLst, zeros=FALSE)
{
  
# New version :  changed so that zeros are included in the raw profile
#                zeros will be removed at an upper level (at 'oneProfile' function)
  
  rawProf <- sapply(ontoLevelsLst, function(x) sum(ancestorsLst==x))
  
#  rawProf <- table(ancestorsLst)
#   rawProf <- rawProf[names(rawProf) %in% ontoLevelsLst]
#   if(zeros==TRUE){
#        noElements <-ontoLevelsLst[is.na(match(ontoLevelsLst,ancestorsLst))]
#        emptyCats<-vector(mode="numeric", length=length(noElements))
#        names(emptyCats) <- noElements
#        rawProf<-c(rawProf, emptyCats)
#        rawProf<-rawProf[order(names(rawProf))]}
   return(rawProf)
}

