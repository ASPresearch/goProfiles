`expandTerm` <-
function (GOTerm, onto){
  if (!((onto=="MF")||(onto=="BP")||(onto=="CC"))){
        on.exit(cat("You must enter MF, BP or CC as ontology"))
        expandedTerm<-NULL} # Wouldn't it be better to use NA instead?
 else{
   switch(onto,
      "MF" = {children <- GOMFCHILDREN},
      "BP" = {children <- GOBPCHILDREN},
      "CC" = {children <- GOCCCHILDREN})
       expandedTerm <-unique(unlist(AnnotationDbi::mget(GOTerm,children,ifnotfound=NA)))
      }
 return(expandedTerm)
}

