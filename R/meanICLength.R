`meanICLength` <-
function(simObj, nsims = simObj$simDesign$nsims, precis = 0.95) {
    2 * qnorm(p=(1-precis)/2, lower.tail=F) *
    c(
         mean(simObj$se),
         sd(simObj$se)/sqrt(nsims)
    )
}

