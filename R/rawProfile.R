`rawProfile` <-
function (ancestorsLst, ontoLevelsLst, zeros=FALSE)
{
   rawProf <- table(ancestorsLst)
   rawProf <- rawProf[names(rawProf) %in% ontoLevelsLst]
   if(zeros==TRUE){
        noElements <-ontoLevelsLst[is.na(match(ontoLevelsLst,ancestorsLst))]
        emptyCats<-vector(mode="numeric", length=length(noElements))
        names(emptyCats) <- noElements
        rawProf<-c(rawProf, emptyCats)
        rawProf<-rawProf[order(names(rawProf))]}
   return(rawProf)
}

