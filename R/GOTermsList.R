`GOTermsList` <-
function (LLids, onto="any", evid="any", na.rm=TRUE, orgPkg){

stopifnot(require(orgPkg, character.only = TRUE, quietly=TRUE, warn.conflicts=FALSE))
#
# Old version used environments from GO package to be deprecated
## mylist<-mget(as.character(LLids),GOENTREZID2GO,ifnotfound=NA)
#
# New version: uses mget wrappers from AnnotationDbi
envirName <-paste(substr(orgPkg,1,nchar(orgPkg)-3),"GO",sep="")
myGOenvir<-eval(parse(text = envirName))
mylist<-AnnotationDbi::mget(as.character(LLids),myGOenvir, ifnotfound=NA)
### Alternatively to mget...--> use toTable/sql
### got <- toTable(myGOenvir)

  GOIDList <-mylist
	for (i in 1:length(mylist))
	{	if(!is.na(mylist[[i]][1])){
			GOIDList[[i]] <-character()
			for (j in 1:length(mylist[[i]])){
				GOIDList[[i]][j]<- mylist[[i]][[j]][["GOID"]]
				names(GOIDList[[i]])[j] <-paste (
					mylist[[i]][[j]][["Ontology"]],
					mylist[[i]][[j]][["Evidence"]],
					sep="-")
				}
		} else{
			 GOIDList[[i]]<-NA}
	}
	names(GOIDList)<-names(mylist)

# Filter out terms according to ontology and evidence
  
if (!onto=="any")
  GOIDList <-sapply(GOIDList, 
    function(l) l[ sapply(names(l), function(x) substr(x,1,2))==onto])
if (!evid=="any")
  GOIDList <-sapply(GOIDList, 
    function(l) l[ sapply(names(l), function(x) substr(x,4,nchar(x)))==evid])

for(i in 1:length(GOIDList)){
  if (length(GOIDList[[i]])==0) GOIDList[[i]]<-NA
  if (is.null(GOIDList[[i]][[1]])) GOIDList[i]<-NA
}

 NAs<-sapply(GOIDList, function(l) is.na(l[[1]]))
 num.NAs <- sum(NAs)

if (na.rm){
  GOIDList<-GOIDList[!NAs]
}

attr(GOIDList,"numNAs") <- num.NAs
attr(GOIDList, "numGenes") <- length(LLids)-num.NAs
return(GOIDList)
}

