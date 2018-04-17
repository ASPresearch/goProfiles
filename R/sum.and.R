`sum.and` <-
function(i,j, logicMat, pGO) {
    return(sum(pGO[logicMat[i,]&logicMat[j,]]))
}

