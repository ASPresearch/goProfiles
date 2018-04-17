`contrMatrix` <-
function(ncateg, simult = ncateg) {
        if (ncateg < simult) stop("")
        result <- cbind(diag(ncateg), matrix(NA, ncol=lenGO(ncateg, simult) - ncateg, nrow=ncateg))
        i1 <- 1
        i2 <- ncateg
        kcols <- ncateg
        colnames(result) <- as.character(1:ncol(result))
        for (k in 2:simult) {
                kcols <- kcols * (ncateg - k + 1) / k
                j1 <- i2 + 1
                j2 <- i2 + kcols
                result[,j1:j2] <- partRes <- completeMat(i1,i2,result, ncateg)
                colnames(result)[j1:j2] <- colnames(partRes)
                i1 <- j1
                i2 <- j2
        }
        result
}

