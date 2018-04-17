`haveGO` <-
function (nams, probeType="Entrez", orgPkg, anotPkg=NULL)
{if (probeType=="Entrez"){
	stopifnot(require(orgPkg, character.only = TRUE, quietly=TRUE, warn.conflicts=FALSE))
	envirName <-paste(substr(orgPkg,1,nchar(orgPkg)-3),"GO",sep="")
	myGOenvir<-eval(parse(text = envirName))
	mylist<-mget(as.character(nams),myGOenvir, ifnotfound=NA)
	returnValue <- sapply(mylist,
                 function(x) if (length(x) == 1 && is.na(x)) FALSE else TRUE)
  }else{
    if (!is.null(anotPkg)){
         stopifnot(require(anotPkg, character.only=TRUE, quietly=TRUE, warn.conflicts=FALSE))
         myenvirGO<-eval(parse(text = paste(anotPkg,"GO",sep="")))
         returnValue <- sapply(mget(nams, myenvirGO,ifnotfound=NA),
                 function(x) if (length(x) == 1 && is.na(x)) FALSE else TRUE)
     }else{
         returnValue <- NA}
        }
 returnValue
}

