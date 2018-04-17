`plotProfiles` <-
function(aProf, aTitle="Functional Profile", anOnto=NULL, percentage=FALSE,
                    HORIZVERT=TRUE, legendText=NULL, colores=c("white", "red"),
                    multiplePlots=F, multipleWindows =T, labelWidth=25,...){
###
### Check if there is one or two columns in the profile
  if (is.data.frame(aProf)){
    numProfs <-ncol(aProf)-2
    numCats <- nrow(aProf)
  }else{
    numProfs <-ncol(aProf[[1]])-2
    numCats <- nrow(aProf[[1]])}
  if (numProfs==1){
    plot1Prof ( aProf=aProf, aTitle=aTitle, anOnto=anOnto, 
                percentage = percentage, HORIZVERT=HORIZVERT, 
                legendText=legendText, colores=rainbow(numCats),
                multiplePlots=multiplePlots, labelWidth=labelWidth,...)
  }else{
      plot2Prof ( aProf=aProf, aTitle=aTitle, anOnto=anOnto, 
                percentage = percentage, HORIZVERT=HORIZVERT, 
                legendText=legendText, colores=colores,
                multiplePlots=multiplePlots, labelWidth=labelWidth,...)}
}

`plot1Prof` <-
function(aProf, aTitle="Functional Profile", anOnto=NULL, percentage= FALSE,
                  HORIZVERT=TRUE, legendText, colores=rainbow(16),
                  multiplePlots=F, multipleWindows =F, labelWidth=25,...) {
  if (is.data.frame(aProf)){
    opt<-par (mar=c(4,12,4,4),las = 2, xpd=TRUE)
    plotOne ( aProf=aProf, aTitle=aTitle, anOnto=anOnto, percentage=percentage,
              HORIZVERT=TRUE, legendText=legendText, colores=colores, labelWidth=labelWidth,...)
    par(opt)
  }else{
      if (multiplePlots & (!multipleWindows))
        opt<-par (mar=c(4,12,4,4), xpd=TRUE, mfrow=c(3,1))
      else
        opt<-par (mar=c(4,12,4,4), xpd=TRUE)
        for (i in 1:length(aProf)){
        if ((multipleWindows) & length(aProf)>1 & i >1 ) dev.new() 
        plotOne ( aProf[[i]], aTitle=aTitle, anOnto=names(aProf[i]),
                  percentage=percentage, HORIZVERT=TRUE, 
                  legendText=legendText, colores=colores, labelWidth=labelWidth,...)}
    par(opt)  
   # mtext(aTitle,line=3,adj=0,cex=1.2)
  }    
}

`plotOne` <-
function (aProf, aTitle="Functional Profiles", anOnto=NULL, percentage=FALSE,
                    HORIZVERT=TRUE, legendText, colores,  labelWidth=25,...){
    freq <-aProf$Frequency
    desc <-as.character(aProf$Description)
    if (percentage){
        numGenes<-attr(aProf,"numGenes")
        if (!is.null(numGenes)& !(numGenes==0))
            freq<-round((freq/numGenes*100),1)
        bp<-barplot(freq, horiz=HORIZVERT,  beside=T, 
            legend.text =legendText,
            col= colores, xlim=c(0,100),...)}
    else
        bp<-barplot(freq, horiz=HORIZVERT, beside=T, 
            legend.text =legendText, col= colores,...)
    text(freq,round(bp,1),freq,pos=4,cex=0.8) 
    axis(2, at = bp, labels = shortName(desc, labelWidth), cex.axis = 0.9, las=2)
    if (is.null(anOnto))
        {title(aTitle)}
    else
        {title(main=paste(aTitle,". ", anOnto, " ontology",sep=""))} 
 }

`plot2Prof` <-
function(aProf, aTitle="Functional Profile", anOnto=NULL, percentage=FALSE,
                    HORIZVERT=TRUE, legendText=NULL, colores=c("white", "red"),
                    multiplePlots=F, multipleWindows =F, labelWidth=25,...){
  if (is.data.frame(aProf)){
    opt<-par (mar=c(4,12,4,4), xpd=TRUE,cex.axis=0.01)
    plotTwo ( aProf=aProf, aTitle=aTitle, anOnto=anOnto, 
              percentage = percentage, HORIZVERT=HORIZVERT, 
              legendText=legendText, colores=colores, labelWidth=labelWidth, ...)
    par(opt)
  }else{
    if (multiplePlots& (!multipleWindows))
      opt<-par (mar=c(4,12,4,4), xpd=TRUE, mfrow=c(3,1))
    else
      opt<-par (mar=c(4,12,4,4), xpd=TRUE)
    for (i in 1:length(aProf)){
      if ((multipleWindows) & length(aProf)>1 & i >1 ) dev.new()
      plotTwo (aProf[[i]], aTitle=aTitle, anOnto=names(aProf[i]),
            percentage=percentage, HORIZVERT=HORIZVERT, 
            legendText=legendText, colores=colores,labelWidth=labelWidth,...)}
    par(opt)  
    #mtext(aTitle,line=3,adj=0,cex=1)
  } 
}

`plotTwo` <-
function(aProf, aTitle="Functional Profiles", anOnto=NULL, percentage=FALSE,
                    HORIZVERT=TRUE, legendText, colores, labelWidth=25,...){
    freq <-t(as.matrix(aProf[,3:4]))
    desc <-as.character(aProf$Description)
    opt<-par (mar=c(4,12,4,4), xpd=TRUE,cex.axis=0.01) 
    if (percentage){
        numGenes1<-attr(aProf,"numGenes1")
        numGenes2<-attr(aProf,"numGenes2")
        if (!is.null(numGenes1)& !(numGenes1==0) &
            !is.null(numGenes2)& !(numGenes2==0))
            freq[1,]<-round((freq[1,]/numGenes1*100),2)
            freq[2,]<-round((freq[2,]/numGenes2*100),2)
        bp<-barplot(freq, horiz=HORIZVERT, beside=T, 
            legend.text =legendText,
            col= colores, 
            xlim=c(0,100),...)}
    else
         bp<-barplot(freq, horiz=HORIZVERT, beside=T, 
            legend.text =legendText,
            col= colores, ...)
    text(freq, round(bp,1),freq,pos=4,cex=0.6) 
    axis(1,cex.axis=0.8, labels=seq(0,100,by=20),at=seq(0,100,by=20))
    axis(2,at=(bp[1,]+bp[2,])/2, labels=shortName(desc, labelWidth), cex.axis = 0.6, las=2)
    if (is.null(anOnto))
        {title(aTitle)}
    else
        {title(main=paste(aTitle,". ", anOnto, " ontology",sep=""))}
   par(opt)
}
