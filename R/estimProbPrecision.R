`estimProbPrecision` <-
function(pEstim, n, alpha = 0.025){
# 1 - 2*alpha confidence limits for the estimated (from n Bernouilli trials) probability pEstim
        return(qnorm(1-alpha) * sqrt(pEstim*(1-pEstim)/n))
}

