`completeMat` <-
function(i1, i2, mat, ncateg) {
        resList <- lapply(i1:i2, expandCol, mat, ncateg)
        colNames <- unlist(lapply(resList, colnames))
        result <- matrix(unlist(resList),nrow=ncateg)
        colnames(result) <- colNames
        return(result)
}

