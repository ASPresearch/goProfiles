`completeNames` <-
function(i1, i2, namesVec, ncateg) {
        unlist(lapply(i1:i2, expandNames, namesVec, ncateg))
}

