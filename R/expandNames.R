`expandNames` <-
function(icol, namesVec, ncateg) {
        colName <- namesVec[icol]
        colIndexes <- unlist(strsplit(colName,"\\."))
        indexCol <- as.integer(colIndexes[length(colIndexes)])
        if (indexCol < ncateg)
                result <- paste(colName,as.character((indexCol+1):ncateg),sep=".") 
        else
                result <- NULL
        result
}

