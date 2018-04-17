dEuclid2 <- function(x, y) {
  UseMethod("dEuclid2")
}

dEuclid2.default <- function(x, y) {
  return(sum((x-y)^2))
}

dEuclid2.BasicGOProfile <- function(x, y) {
  return(sum((x[,3]/ngenes(x) - y[,3]/ngenes(y))^2))
}
