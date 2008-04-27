`internal.covGO` <-
function(pGO, pij. = secondOrderMarginals(pGO)){
        pi. <- diag(pij.)
        return(pij. - pi. %o% pi.)
}

