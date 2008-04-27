`fullGOProfile` <-
function(pGO, colNames) {
    if (!is.null(pGO)) {
        n <- attr(pGO,"ngenes")
        pGO <- pGO[pmatch(colNames,names(pGO))]
        naPositions <- is.na(pGO)
        names(pGO)[naPositions] <- colNames[naPositions]
        pGO[naPositions] <- 0
        attr(pGO,"ngenes") <- n
    }
    return(pGO)
}

