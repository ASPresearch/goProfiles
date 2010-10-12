
internal.enrichProfile <- function(prof, qrof, n, m, method, ...) {
  s <- length(prof)
  pvals <- numeric(s)
  for (i in 1:s){
    freqs_i <- trunc(matrix(c(prof[i], n - prof[i], qrof[i], m - qrof[i]), nrow=2))
    pvals[i] <- fisher.test(freqs_i)$p.value
  }
  names(pvals) <- names(prof)
  adjPvals <- p.adjust(pvals, method=method, ...)
  attr(adjPvals, "unadjusted") <- pvals
  return(adjPvals)
}
