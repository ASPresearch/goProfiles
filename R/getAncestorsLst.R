`getAncestorsLst` <-
function (GOtermslist, onto, unique.ancestor=TRUE, na.rm=TRUE,
                            combine=TRUE){
 vector.is.empty <- function(x) return(length(x) ==0 )
 if (!onto %in% (c("MF", "BP", "CC"))){
      on.exit(cat ("Invalid ontology code. Must be 'MF', 'BP' or 'CC'"))
      return(0)
 }else{
   envirName <-paste("GO",onto,"ANCESTOR",sep="")
   envi<-eval(parse(text = envirName))
 }
### Check (and force)that GOterms lists is formed only be terms of "onto" ontology
 GOtermslist <-sapply(GOtermslist, 
    function(l) l[ sapply(names(l), function(x) substr(x,1,2))==onto],
    simplify=FALSE)
### Create and fill ancestors list
 numAncestors <- length(GOtermslist)
 AncestorsLst <- vector("list", numAncestors)
 for (i in 1:numAncestors){
        x <- GOtermslist[[i]]                                     
        x <- x[!is.na(x)]
        if ((length(x)>0) && (!vector.is.empty(x))) {
            AncestorsLst[[i]]<-unlist(AnnotationDbi::mget(as.character(x),envi,ifnotfound=NA))
            if (combine) AncestorsLst[[i]] <- c(AncestorsLst[[i]], x)
            if (unique.ancestor){
              AncestorsLst[[i]]<-unique(AncestorsLst[[i]]) [-1]} # remove "all"
            names(AncestorsLst)[i]<-names(GOtermslist)[i]       
        }else{ 
            AncestorsLst[[i]]<-NA}
  }
# Falta eliminar  els NA si na.rm=TRUE
# Aqui hi ha un problema 
removeNAs<-function (li){
  li2<-list()
  for (i in 1:length(li))
    if (!is.na(li[i]))
       li2<-c(li2,li[i])
  return (li2)}
if (na.rm) 
  AncestorsLst <-removeNAs(AncestorsLst)
# Si un cop trets tots els NAs no queda res el posem a NULL
if (sum(!is.na(AncestorsLst)) == 0) AncestorsLst <- NULL
 return (AncestorsLst) 
}

