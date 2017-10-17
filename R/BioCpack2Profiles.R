`BioCpack2Profiles` <-
function (anotPkg, orgPackage, level=2, na.rm=TRUE, expanded=FALSE){
    stopifnot(require(anotPkg, character.only=TRUE, quietly=TRUE, warn.conflicts=FALSE))
    stopifnot(require(orgPackage, character.only=TRUE, quietly=TRUE, warn.conflicts=FALSE))

    env<-paste(anotPkg, "ENTREZID",sep="")
    env <- eval(parse(text = env))

    gLL<-as.character(unlist(as.list(env)))
    if(na.rm)
      gLL<-gLL[!is.na(gLL)]

    all.Profiles <-basicProfile(genelist=gLL, onto="ANY",
        level=level, empty.cats=TRUE, cat.names=TRUE, na.rm=na.rm, orgPackage=orgPackage)
    all.expanded.Profiles <-NULL
    if (expanded)
      all.expanded.Profiles <- expandedProfile (gLL, onto="ANY",level=level,
        na.rm=na.rm, orgPackage=orgPackage)
    return (list (profiles=all.Profiles, expanded.profiles=all.expanded.Profiles))
 }

