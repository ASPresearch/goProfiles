`dgetGO` <-
function(path="", fileName){
### POSSIBLEMENT A SUPRIMIR
        gObject <- dget(paste(path,fileName,sep=""))
        goVector <- gObject[!is.na(names(gObject))]
        goSum <- sum(goVector)
        goVector <- goVector / goSum
        n <- attr(gObject,"n")
        if (is.null(n))
                attr(goVector,"ngenes") <- goSum
        else
                attr(goVector,"ngenes") <- n
        return(goVector)
}

