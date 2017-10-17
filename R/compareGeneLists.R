`compareGeneLists` <-
function (genelist1, genelist2, idType="Entrez", onto="ANY", level=2, orgPackage,
	  method = "lcombChisq", ab.approx = "asymptotic", confidence = 0.95,
	  compareFunction="compareGOProfiles",...)
{
  ### Step 1: Compute expanded profiles
  
  if (!(idType %in% c("Entrez", "BioCprobes", "GOTermsFrame"))){
    compared<-NULL
    on.exit(cat("Gene identifiers are not understood by the program"))
    }
  if (!((onto=="ANY")||(onto=="MF")||(onto=="BP")||(onto=="CC"))){
    compared<-NULL
    on.exit(cat("You must enter MF, BP or CC as ontology"))
    }
  expanded1 <- expandedProfile (genelist1, idType=idType, onto=onto, level=level, orgPackage=orgPackage, ...) 
  expanded2 <- expandedProfile (genelist2, idType=idType, onto=onto, level=level, orgPackage=orgPackage, ...) 

  ### Step 2: Determine how is the relation  between genelist1 and genelist2
  commonGenes <-intersect(genelist1, genelist2)
  if (length(commonGenes)==0){
    relationType <-"DISJOINT"
  }else{
    if (length(commonGenes)==min(length(genelist1), length(genelist2))){
      relationType <-"INCLUSION"
    }else{
      relationType <-"INTERSECTION"}}
  
  ### Step 3: Compare profiles
  
  # compared <-list()
  numElements <- length(expanded1)
  compared <- vector("list", numElements)
  for(i in 1:numElements){
    if (relationType=="DISJOINT"){
        compared[[i]]<-compareGOProfiles (pn=expanded1[[i]], qm=expanded2[[i]],
                          method=method, ab.approx=ab.approx, confidence=confidence,...)
      }else{
        if (relationType=="INCLUSION"){
            if (length(genelist1) <=length(genelist2)){
              populationProfile <-expanded2; sampleProfile <-expanded1
            }else{
              populationProfile <-expanded1; sampleProfile <-expanded2}
              if (pmatch(compareFunction, "compareGOProfiles")){
            		compared[[i]]<-compareGOProfiles (pn=populationProfile[[i]], 
                              pqn0=sampleProfile[[i]], method=method, ab.approx=ab.approx, 
                              confidence=confidence,...)
              }else{compared[[i]]<-fitGOProfile (p0=populationProfile[[i]], 
                              	pn=sampleProfile[[i]],
                              	method=method, ab.approx=ab.approx, confidence=confidence,...)
              }
        }else{
                common.expanded <-expandedProfile(commonGenes, idType=idType, onto=onto, level=level, orgPackage=orgPackage,...)
                compared[[i]]<-compareGOProfiles (pn=expanded1[[i]], 
                          qm=expanded2[[i]], pqn0=common.expanded[[i]],
                          method=method, ab.approx=ab.approx, confidence=confidence,...)}
                }
        names(compared)[i]<-attr(expanded1[[i]],"onto")
        }
return(compared)
}  

