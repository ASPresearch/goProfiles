`niceProfile` <-
function(profTable){
    len<-length(names(profTable))
    desc<-as.character(rep(NA,len))
    got<-toTable(GOTERM)[,2:3]
    for (i in 1:len){
        # godsc<-getGOdesc(names(profTable)[i],"ANY")
        # getGOdesc will be deprecated. Omit its use.
        goterm<-names(profTable)[i]
        desc[i]<-ifelse((!is.na(goterm)),got[got[,1]==goterm,2],"NULL")
     }
    niceProf<-data.frame("Description"= desc, "GOID"= names(profTable), 
        "Frequency"= as.vector(profTable))
    return(niceProf)
}

