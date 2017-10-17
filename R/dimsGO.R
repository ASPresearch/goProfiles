dimsGO <- function(prof) {
    UseMethod("dimsGO")
}

dimsGO.ExpandedGOProfile <- function(prof) {
    dimsGO.character(rownames(prof))
}

dimsGO.character <- function(categNamesGO) {
        categsList <- strsplit(categNamesGO,"\\.")
        ncateg <- max(as.integer(unlist(categsList)))
        simult <- max(sapply(categsList, length))
        return(list(ncateg=ncateg, simult=simult))
}

