`as.ExpandedGOProfile.ExpandedGOProfile` <-
function(x, expandedCatNames=NULL) {
    if (!is.null(expandedCatNames))
        row.names(x) <- expandedCatNames
    return(x)
}

