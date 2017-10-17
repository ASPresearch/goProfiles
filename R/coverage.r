`coverage` <-
function(icoverage, t.stat, signifLimit) {
    ifelse(abs(t.stat) <= signifLimit[icoverage], 1, 0)
}

