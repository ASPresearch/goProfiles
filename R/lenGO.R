`lenGO` <-
function(ncateg = 8, simult = ncateg) {
        ncols <- ncateg
        kcols <- ncols
        for (k in 2:simult) {
                kcols <- kcols * (ncateg - k + 1) / k
                ncols <- ncols + kcols 
        }
        return(ncols)
}

