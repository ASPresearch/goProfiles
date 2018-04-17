`shortName` <- 
function (aString, aWidth){
    numCh <- nchar(aString)
    aString2 <- substr(aString, 1, aWidth) # Limit names of terms to 'aWidth' characters
    aString3 <- paste(aString2, ifelse(numCh > (aWidth-5), "...", ""), sep="")
    return(aString3)
}

