`basicProfile` <-
function (genelist, idType="Entrez", onto="ANY", level=2, orgPackage=NULL, anotPackage =NULL, 
	ord=TRUE, multilevels=NULL, empty.cats = TRUE, cat.names = TRUE, na.rm=TRUE) 

{
oneProfile<-function(GOTermsList, onto, level=2, multilevels=NULL, 
                     empty.cats = FALSE, cat.names = FALSE){
  funcProfile<-NULL
  ancestorsList <-getAncestorsLst(GOTermsList,onto) 

  if (!is.null(ancestorsList)){
      ancestors<-unlist(ancestorsList)
      if (is.null(multilevels))
        ontoLevel<- getGOLevel (onto,level)
      else
        ontoLevel<- multilevels # getGOLevel (onto,level)
      if (is.null(ontoLevel))
        on.exit(cat("No list of terms available for this ontology and level"))
      else
        funcProfile<-rawProfile(ancestors, ontoLevel, empty.cats)
        # recall that empty.cats is not a valid parameter for rawProfile anymore
      if (cat.names){
            funcProfile<-niceProfile(funcProfile)
            if (ord)
              funcProfile<-funcProfile[order(funcProfile$Description),]
            }
     }else 
            {on.exit(cat("No ancestors found for this GOTermsList"))}
 
  if (!is.null(funcProfile)){
    if(!empty.cats){
      if(cat.names)
        funcProfile <-funcProfile[funcProfile$Frequency!=0,] # niceProfile yields a data.frame
      else
        funcProfile <-funcProfile[funcProfile!=0] # rawProfile yields a vector
    }
    attr(funcProfile,"numGenes")<-length(ancestorsList)
                                        # use length(ancestorsList) instead of length(GOTermsList) to account for the cases where ancestors are missing in "onto"
                                        # attr(GOTermsList,"numGenes")
    attr(funcProfile,"numNAs")<-attr(GOTermsList,"numNAs")
    attr(funcProfile,"ontology")<-onto}
  class(funcProfile) <-c("BasicGOProfile", class(funcProfile))
return(funcProfile)
}
    if (idType %in% c("Entrez", "BioCprobes", "GOTermsFrame"))
        GOList <- as.GOTerms.list (genelist, idType, orgPackage, anotPackage,
                                   na.rm=na.rm)
    else
         {on.exit(cat("Gene identifiers are not understood by the program"))
        myprofile<-NULL}
    if (!((onto=="ANY")||(onto=="MF")||(onto=="BP")||(onto=="CC")))
        {on.exit(cat("You must enter MF, BP or CC as ontology"))
        myprofile<-NULL}
    else{
      if (onto=="ANY"){
          ontonames<-c("MF","BP","CC")
      }else{
          ontonames<-onto}    
      names(ontonames)<-ontonames
      myprofile<-lapply (ontonames, function (ontology){
                  oneProfile(GOTermsList=GOList,onto=ontology,
                              level=level, multilevels=multilevels,
                              empty.cats=empty.cats,cat.names=cat.names) }
                         )
      }
return (myprofile)
}

