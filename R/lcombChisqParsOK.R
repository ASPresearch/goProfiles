`lcombChisqParsOK` <-
function(nchis, df, ncp) {
    if (length(df) > 1)
        if (nchis != length(df)) {
            message("The mixture coefficients vector and the df vector must have the same length")
            return(FALSE)
        }
    if (length(ncp) > 1)
        if (nchis != length(ncp)) {
            message("The mixture coefficients vector and the ncp vector must have the same length")
            return(FALSE)
        }
    return(TRUE)
}

