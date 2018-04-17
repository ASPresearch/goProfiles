`generate.multinomial` <-
function(distri, ssize = 1)
{
# Generates a multinomial distribution by a conditional algorithm
#
        lprob <- length(distri$p)
        X <- matrix(distri$n, nrow = ssize, ncol = lprob, byrow = T)
        colnames(X) <- names(distri$p)
        if(distri$n > 0)
                if(lprob > 1) {
                        X[,1] <- rbinom(ssize,  distri$n, distri$p[1])
                        if(lprob > 2) {
                                condProb <- matrix(distri$condProb,
                                        nrow = ssize, ncol = length(distri$condProb),
                                        byrow = T)
                                size <- rep(distri$n, ssize)
                                for(i in 2:(lprob - 1)) {
                                  size <- size - X[,i-1]
                                  #X[,i] <- ifelse(size, rbinom(ssize, size, condProb[,i-1]), 0)
                                  X[,i] <- rbinom(ssize, size, condProb[,i-1])
                                }
                        }
                        X[,lprob] <- size - X[,lprob-1]
                }
        return(X)
}

