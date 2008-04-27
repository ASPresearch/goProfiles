`getGOLevel` <-
function (onto, level){
 if (!((onto=="MF")||(onto=="BP")||(onto=="CC")))
        {on.exit(cat("You must enter MF, BP or CC as ontology"))
        GOLevel<-NULL} 
 else{
   switch(onto,
      "MF" = {topNode <- "GO:0003674"; children <- GOMFCHILDREN},
      "BP" = {topNode <- "GO:0008150"; children <- GOBPCHILDREN},
      "CC" = {topNode <- "GO:0005575"; children <- GOCCCHILDREN})
   #GOLevel <-unique(get(topNode,children))
   GOLevel<-topNode
   if ((level <=1)||(level >5)){
    on.exit(cat("You must enter a level greater than 1 and smaller than 6"))
    GOLevel<-NULL}
   else{
      for (i in (1:(level-1)))
          GOLevel <-unique(unlist(AnnotationDbi::mget(GOLevel,children,ifnotfound=NA)))
      GOLevel<-GOLevel[!is.na(GOLevel)]
      }
   }
return(GOLevel)
}

