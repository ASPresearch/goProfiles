`dimsGO` <-
function(categNamesGO) {
        categsList <- strsplit(categNamesGO,"\\.")
        ncateg <- max(as.integer(unlist(categsList)))
        simult <- max(sapply(categsList, length))
        return(list(ncateg=ncateg, simult=simult))
}

