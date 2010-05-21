`getGOLevel` <-
function (onto, level){
{
    GOChildren <- function (ontoTab, terms) {
    children <- ontoTab[ontoTab[,2] %in% terms,1]
    }

    if (!((onto == "MF") || (onto == "BP") || (onto == "CC"))) {
        on.exit(cat("You must enter MF, BP or CC as ontology"))
        GOLevel <- NULL
    }
    else {
        switch(onto, MF = {
            topNode <- "GO:0003674"
            chTab<-toTable(GOMFCHILDREN)
        }, BP = {
            topNode <- "GO:0008150"
            chTab<-toTable(GOBPCHILDREN)
        }, CC = {
            topNode <- "GO:0005575"
            chTab<-toTable(GOCCCHILDREN)
        })

        GOLevel <- topNode
        if (level < 1) {
            on.exit(cat("You must enter a level greater than 1"))
            GOLevel <- NULL
        }
        else {
          if (level==1) {
             GOLevel <- topNode
          }
          else{
            GOLevel<-topNode
            for (i in (1:(level - 1))){
              # childLevel <- chTab[chTab[,2] %in% GOLevel,1]
              childLevel <- GOChildren(chTab, GOLevel)
              GOLevel <- childLevel
              # GOLevel <- GOLevel[!is.na(GOLevel)] # No se si cal?
            }
          }
        }
          
    return(unique(GOLevel))
     }
  }
}


