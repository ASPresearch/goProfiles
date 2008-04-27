`haveEntrezID` <-
function (probeslist, anotPkg)
{
 stopifnot(require(anotPkg, character.only=TRUE, quietly=TRUE))
 myenvirENTREZID<-eval(parse(text = paste(anotPkg,"ENTREZID",sep="")))
 return(sapply(mget(probeslist, myenvirENTREZID, ifnotfound=NA),
                 function(x) if (length(x) == 1 && is.na(x)) FALSE else TRUE))
}

