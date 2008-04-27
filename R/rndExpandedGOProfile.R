`rndExpandedGOProfile` <-
function(p, n, generator = multinomial(n=n, p=p), ssize=1) {
    as.ExpandedGOProfile(t(generate.multinomial(generator, ssize=ssize)))
}

