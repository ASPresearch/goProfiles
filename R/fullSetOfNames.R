`fullSetOfNames` <-
function(ncateg, simult=ncateg) {
        if (ncateg < simult) stop("")
        result <- as.character(1:lenGO(ncateg, simult))
        i1 <- 1
        i2 <- ncateg
        kcols <- ncateg
        for (k in 2:simult) {
                kcols <- kcols * (ncateg - k + 1) / k
                j1 <- i2 + 1
                j2 <- i2 + kcols
                result[j1:j2] <- completeNames(i1,i2,result, ncateg)
                i1 <- j1
                i2 <- j2
        }
        result
}

