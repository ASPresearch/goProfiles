`contrast` <-
function(icontrast, stat, signifLimit) {
        ifelse(stat >= signifLimit[icontrast], 1, 0)
}

