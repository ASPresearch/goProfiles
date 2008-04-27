`normIntLength` <-
function(mean.se, confid=0.95) 2 * qnorm(p=(1-confid)/2, lower.tail=F) * mean.se

