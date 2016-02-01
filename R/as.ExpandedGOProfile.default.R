`as.ExpandedGOProfile.default` <-
function(x, expandedCatNames=NULL) {
    x <- as.matrix(x)
    result <- as.data.frame(x)
    if (is.null(expandedCatNames))
        rowNames <- row.names(result)
    else
        rowNames <- expandedCatNames
    result <- data.frame(lapply(result, frelGO))
    row.names(result) <- rowNames
    class(result)<- c("ExpandedGOProfile", class(result))
    return(result)
}

