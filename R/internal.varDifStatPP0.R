`internal.varDifStatPP0` <-
function(pn, p0 = NULL, funcProfP0 = NULL){
        allNams <- unique(c(names(pn), names(funcProfP0), names(p0)))
        pn <- fullGOProfile(pn, allNams)
        pij. <- secondOrderMarginals(pn)
        pi. <- diag(pij.)
        if (is.null(funcProfP0))
            funcProfP0 <- contractedProfile(fullGOProfile(p0, allNams))
        else {
            funcProfP0 <- funcProfP0[pmatch(names(pi.),names(funcProfP0))]
            naPositions <- is.na(funcProfP0)
            names(funcProfP0)[naPositions] <- names(pi.)[naPositions]
            funcProfP0[naPositions] <- 0
        }
        return(as.numeric(4 * t(pi.-funcProfP0) %*% internal.covGO(pij. = pij.) %*% (pi.-funcProfP0)))
}

