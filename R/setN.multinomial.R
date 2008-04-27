`setN.multinomial` <-
function(distri, n)
{
# changes the n parameter of a multinomial distribution
        if (n < 0) stop("Negative size (n) binomial parameter")
        distri$n <- n
        return(distri)
}

