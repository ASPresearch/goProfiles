`expandedLevel` <-
function (LevelTerms, Term2Expand, onto)
{ expandedterm<-expandTerm(Term2Expand,onto)
  expandedterm<-c(LevelTerms[-Term2Expand],expandedterm)
  return(expandedterm)
}

