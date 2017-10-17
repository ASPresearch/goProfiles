### R code from vignette source 'goProfiles.Rnw'

###################################################
### code chunk number 1: requiredFiles1
###################################################
options(width=60, warn=0, digits=5)
require(goProfiles)


###################################################
### code chunk number 2: prostateIds
###################################################
require(goProfiles)
data(prostateIds)


###################################################
### code chunk number 3: basicProfilesMF
###################################################
welsh.MF <- basicProfile (welsh01EntrezIDs[1:100], onto="MF", level=2, orgPackage="org.Hs.eg.db") 
singh.MF <- basicProfile (singh01EntrezIDs[1:100], onto="MF", level=2, orgPackage="org.Hs.eg.db") 
welsh.singh.MF <-mergeProfilesLists(welsh.MF, singh.MF, profNames=c("Welsh", "Singh"))
printProfiles(welsh.singh.MF, percentage=TRUE)


###################################################
### code chunk number 4: plotProfileMF
###################################################
plotProfiles (welsh.MF, aTitle="Welsh (2001). Prostate cancer data")


###################################################
### code chunk number 5: basicProfilesANY
###################################################
welsh <- basicProfile (welsh01EntrezIDs[1:100], onto="ANY", level=2, orgPackage="org.Hs.eg.db") 


###################################################
### code chunk number 6: comparevisual
###################################################
plotProfiles (welsh.singh.MF, percentage=T,aTitle="Welsh vs Singh", legend=T) 


###################################################
### code chunk number 7: compareDiff1
###################################################
compared.welsh.singh.01.MF <- compareGeneLists (welsh01EntrezIDs[1:100], singh01EntrezIDs[1:100], onto="MF", level=2, orgPackage="org.Hs.eg.db")
print(compSummary(compared.welsh.singh.01.MF))


###################################################
### code chunk number 8: fisherGOProf1
###################################################
list1 <- welsh01EntrezIDs[1:100]
list2 <- singh01EntrezIDs[1:100]
commProf <- basicProfile(intersect(list1, list2), onto="MF", level=2, orgPackage="org.Hs.eg.db")$MF
fisherGOProfiles(welsh.MF$MF, singh.MF$MF, commProf, method="holm")


###################################################
### code chunk number 9: compareEquiv1
###################################################
data(prostateIds)
expandedWelsh <- expandedProfile(welsh01EntrezIDs[1:100], onto="MF",
                        level=2, orgPackage="org.Hs.eg.db")
expandedSingh <- expandedProfile(singh01EntrezIDs[1:100], onto="MF",
                        level=2, orgPackage="org.Hs.eg.db")
commonGenes <- intersect(welsh01EntrezIDs[1:100], singh01EntrezIDs[1:100])
commonExpanded <- expandedProfile(commonGenes, onto="MF", level=2, orgPackage="org.Hs.eg.db")

equivMF <-equivalentGOProfiles (expandedWelsh[['MF']], 
                          qm  = expandedSingh[['MF']], 
                          pqn0= commonExpanded[['MF']])
print(equivSummary(equivMF, decs=5))



