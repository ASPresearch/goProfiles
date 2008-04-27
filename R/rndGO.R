`rndGO` <-
function(mnom, ssize=10) 
{
        x <- as.data.frame(t(generate.multinomial(mnom, ssize)/mnom$n))
        result <- sapply(x, contractCol, row.names(x))
        attr(result, "ngenes") <- mnom$n
        return(result)
}

