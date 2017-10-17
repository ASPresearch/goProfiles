`toNumVec` <-
function(x){
    charVec<-unlist(strsplit(as.character(x[1]),"\\."))
    numVal<-0
    Len <- length(charVec)                                          
    for (i in 1:Len){
        numVal<-numVal+10^(Len-i)*as.numeric(charVec[i])}
    #return(list(charVec,numVal))
    return(numVal)
        }

