`expandedProfile` <-
function(genelist, idType="Entrez", onto="ANY", level=2,orgPackage=NULL, anotPackage=NULL,    
  multilevels=NULL, ord=TRUE, na.rm=TRUE, percentage=TRUE){


oneProfile<-function(GOTermsList, onto="ANY", level=2, multilevels=NULL,
                    ord=TRUE, na.rm=TRUE, percentage=TRUE){   
    ancestorsList<-getAncestorsLst(GOTermsList,onto)
    if (!is.null(ancestorsList)){
      if (is.null(multilevels))
        ontoLevel<- getGOLevel (onto,level)
      else
        ontoLevel<- multilevels #getGOLevel (onto,level)
    if (is.null(ontoLevel))
     {errorMsg<-paste("No list of terms available for ",
     "Ontology: ",onto," and level: ",level,sep="")
     on.exit(cat(errorMsg))
        expProf<-NULL}
    else{
        numGenes <- length(ancestorsList) 
        numProfiles <- 0
        # exp.profile <- data.frame(frec=rep (0,length(ontoLevel)))
        exp.profile <- matrix(nrow=length(ontoLevel),ncol=numGenes)
        nams<-character(0)
        for (i in 1:length(ancestorsList))
        {   my.profile<-rawProfile(ancestorsList[[i]],ontoLevel, TRUE)
            # empty.cats is ot used anymore in rawProfile so it should work fine
            if (!is.null(my.profile)){
                exp.profile[,i]<- my.profile
                nams<-c(nams,GOTermsList[i])
                numProfiles<-numProfiles+1
            }
        }
        exp.profile<-exp.profile[,1:numProfiles, drop=FALSE]
        colnames(exp.profile)<-nams #c("0",nams)
        rownames(exp.profile)<-sort(ontoLevel)

        res<-apply(exp.profile, 2,
            function(x){paste(which(x!=0),sep="",collapse=".")})

        expProf<-table(res)
        names(expProf)<-sapply(names(expProf),function(s)
          ifelse(s=="","0",s))
        if (ord){
            o<-sapply(names(expProf),toNumVec)
            expProf<-expProf[order(o)]}
      }
    }
    if(!is.null(expProf)){
      if (na.rm) expProf<-expProf[names(expProf)!="0"]
      if (percentage)
          expProf<- as.ExpandedGOProfile(expProf)
      else
          expProf<- as.data.frame(expProf)
      attr(expProf,"numGenes")<-  numGenes # = length(ancestorsList) 
                                        # use length(ancestorsList) instead of length(GOTermsList)
                                        # to account for the cases where ancestors are missing in "onto"
                                        # length(GOTermsList)
      attr(expProf,"ontology")<-onto}
return(expProf)
}
### Function body starts here 
    myProfile <-NULL
    if (idType %in% c("Entrez", "BioCprobes", "GOTermsFrame"))
        GOList <- as.GOTerms.list (genelist,idType,orgPackage, anotPackage)
    else
        on.exit(cat("Gene identifiers are not understood by the program"))
    if (!((onto=="ANY")||(onto=="MF")||(onto=="BP")||(onto=="CC"))){
      on.exit(cat("You must enter MF, BP, CC or ANY as ontology"))
    }else{
      if (onto=="ANY")
          ontonames<-c("MF","BP","CC")
      else
          ontonames<-onto
      names(ontonames)<-ontonames
      myprofile<-lapply (ontonames, function (ontology){
                  oneProfile(GOTermsList=GOList,onto=ontology,
                              level=level, multilevels=multilevels,
                              ord=ord, na.rm=na.rm, percentage=percentage) }
                        )
      }
return (myprofile)
}

