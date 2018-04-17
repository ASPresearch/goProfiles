`equivSummary` <- 
function (l,decs=6) {
  shortSummary <-function (x, decs) {
    equiv<- ifelse(x["conf.int"][[1]][2] < x["parameter"][[1]][1],1,0)
    results<-c(x["estimate"][[1]], attr(x["estimate"][[1]], "se"), x["p.value"][[1]],  x["conf.int"][[1]][2],x["parameter"][[1]][1],equiv)
    confLevel <- paste(attr(x$conf.int,"conf.level"),"CI", sep="")
    upper <-paste(confLevel, "up",sep=".")
    names(results)<-c("Sqr.Euc.Dist", "StdErr", "pValue",  upper, "d0", "Equivalent? (1=yes,0=not)")
    return(round(results,decs))
    }
 if (!is.list(l[[1]])){
    result<-shortSummary(l,decs)
  }else{
    if (length(l)>1){
    result<-data.frame(tmp=rep(NA,6))
    for (i in 1:length(l)){
        result <- cbind(result, shortSummary(l[[i]],decs))
        colnames(result)[i+1] <- names(l)[i]  
    }
    rownames(result)<-names(shortSummary(l[[i]],decs))
    result<-result[,-1]
  }else{
    result<-shortSummary(l[[1]],decs)
  }}
return(result)
}
