`expandCol` <-
function(icol, mat, ncateg) {
        colName <- colnames(mat)[icol]
        colIndexes <- unlist(strsplit(colName,"\\."))
        indexCol <- as.integer(colIndexes[length(colIndexes)])
        if (indexCol < ncateg) {
                remRows <- (nr <- nrow(mat)) - indexCol
                result <- rbind(matrix(mat[1:indexCol,icol],nrow=indexCol,ncol=remRows), diag(remRows))
                colnames(result) <- paste(colName,as.character((indexCol+1):nr),sep=".")
        }
        else
                result <- NULL
        result
}

