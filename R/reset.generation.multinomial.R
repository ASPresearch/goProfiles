`reset.generation.multinomial` <-
function(distri)
{
# Resets internal parameters to generate a multinomial distribution 
# by a conditional algorithm
#
        siz <- distri$n
        if (siz < 0) stop("Negative multinomial size parameter")
        prob <- distri$p
        if (any(prob < 0)) stop("Negative multinomial probability parameters")
        sprob <- sum(prob)
        if (sprob != 1.0) {
                if (sprob == 0.0) stop("All multinomial probability parameters are zero")
                prob <- prob / sprob
        }
        lprob <- length(prob)
        if(siz > 0)
                if(lprob > 2) {
                        condProb <- rep(1, lprob - 2)
                        condIndexes <- 2:(lprob - 1)
                        condProb <- condProb - cumsum(prob)[1:(lprob - 2)]
                        condProb <- ifelse(condProb, prob[condIndexes]/condProb, 1)
                        condProb <- pmin(condProb,1)
                        distri$condProb <- condProb
                        names(distri$condProb) <- names(prob)[condIndexes]
                }
                else distri$condProb <- NULL
        return(distri)
}

