`as.GOTerms.list` <-
function (genelist, probeType, orgPackage=NULL, anotPkg=NULL, onto="any", na.rm=FALSE)
{
  switch (probeType,
  "Entrez" =
    {
      if(is.null(orgPackage)){
        stop("'orgPackage' cannot be null. An organism annotation package (org.Xx.eg.db) must be provided")
      }else{
        x<-GOTermsList (genelist, onto=onto, na.rm=na.rm, orgPkg=orgPackage)
      }
    },
  "BioCprobes" =
    {
      if(is.null(anotPkg)){
        stop("'anotPkg' cannot be null. A microarray annotation package (CHIPTYPE.db) must be provided")
      }else{
        myLLinkIDs <- BioCprobes2Entrez (probeslist = genelist, anotPkg =anotPkg,
        na.rm=na.rm)}
      if(is.null(orgPackage)){
        stop("'orgPackage' cannot be null. An organism annotation package (org.Xx.eg.db) must be provided")
      }else{
        x <- GOTermsList (myLLinkIDs, onto=onto, na.rm=na.rm, orgPkg=orgPackage)}
    },
  "GOTermsFrame" =
    {
    x<-GOTermsFrame2GOTermsList (genelist)
    })
  return(x)
}

