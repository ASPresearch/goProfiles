`multinomial` <-
function(n, p) {
        mnomial <- list(n = n, p = p, condProb = NULL)
        reset.generation.multinomial(mnomial)
}

