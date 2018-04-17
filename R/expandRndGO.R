`expandRndGO` <-
function(prevRndGO, mnom) 
{
        ssize <- ncol(prevRndGO)
        n <- attr(prevRndGO,"ngenes")
        x <- as.data.frame(t(generate.multinomial(mnom, ssize)))
        result <- (prevRndGO*n + sapply(x, contractCol, row.names(x))) / (n2 <- n + mnom$n)
        attr(result, "ngenes") <- n2
        return(result)
}

