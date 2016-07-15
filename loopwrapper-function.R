library(FD)
library(plyr)
setwd("C:/Users/Clarrien/Documents/traitR")
allTraits <-read.csv("C:/Users/Clarrien/Documents/traitR/allTraits_Final_colsswitched.csv",  strip.white = TRUE,    na.strings = c("NA", ""))
allTraits<-read.csv("~/shinyR/allTraits.csv",na.strings = c("NA", ""))
dim(allTraits)
allTraits.noGnor<-subset(allTraits, GNOR=="una")
allTraits<-allTraits.noGnor
dim(allTraits)

setwd("C:/Users/Clarrien/Documents/traitR/multiaxestraitR")
source("~/TraitR_github/traitR/traitR_standalone_multiple_pcs_works.R")

allTraits$May24.N.leaves <- as.numeric(allTraits$May24.N.leaves)
allTraits$O.June7.Stem.colour<-ordered(allTraits$O.June7.Stem.colour)
allTraits$O.June7.Hairiness<-ordered(allTraits$O.June7.Hairiness)
allTraits$O.July5.Stem.colour<-ordered(allTraits$O.July5.Stem.colour)
allTraits$O.July5.Hairiness<-ordered(allTraits$O.July5.Hairiness)
str(allTraits)

May.PCA.AST <- colnames(allTraits[c(3:4)])
June.PCoA.AST <- colnames(allTraits[c(8:14)])
July.PCoA.Multi <- colnames(allTraits[c(17:23)])
Late.PCA.test <- colnames(allTraits[c(30, 32, 38, 39, 40,41,42,46:52)])#PcA

#RHOP
Late.PCA.RHOP_noBranch <- colnames(allTraits[c(30, 32, 38, 40,41,42,46:52)])#remove branching

#GNOR
Gnor.June.PCoA.multi3.noF<- colnames(allTraits[c(8,9,11,14)])
Gnor.July.PCoA.multi3.noF <- colnames(allTraits[c(17,20:22)])
Gnor.phen.numeric.multi3<- colnames(allTraits[c(30, 32, 46:52)])
Gnor.phen.numeric.multi3.noForce.no542<- colnames(allTraits[c(32, 46:52)])

#Eur reduced
June.reduced<- colnames(allTraits[c(9:12)])
Late.sameasJune2 <-colnames(allTraits[c(30,32,38,47)])
Late.PCA.Eur.Reduced3<- colnames(allTraits[c(32,38,47, 52)])#PcA
Late.PCA.Eur.RedEur2<- colnames(allTraits[c(32,38,47)])#PcA

#same traits only
June.few<-colnames(allTraits[c(9:12)])
Late.few<-colnames(allTraits[c(30, 32, 38, 46)])

#PCoA
list1 <- list(Gnor.June.PCoA.multi3.noF=Gnor.June.PCoA.multi3.noF, Gnor.July.PCoA.multi3.noF=Gnor.July.PCoA.multi3.noF)
list1 <- list(June.PCoA.multi3F=June.PCoA.multi3F)
list1 <- list(June.few=June.few)
list1 <- list(June.reduced=June.reduced)
list1 <- list(June.PCoA.AST=June.PCoA.AST)
#PCA
list1 <- list(Late.PCA.multi3=Late.PCA.multi3)
list1 <- list(Gnor.phen.numeric.multi3.noForce=Gnor.phen.numeric.multi3.noForce)
list1 <- list(May.PCA.multi.all=May.PCA.multi.all)
list1 <- list(Late.PCA.Eur.Reduced2=Late.PCA.Eur.Reduced2)
list1 <- list(Late.sameasJune2=Late.sameasJune2)
list1 <- list(May.PCA.AST=May.PCA.AST)
list1 <- list(Late.PCA.test=Late.PCA.test)
list2 <- colnames(allTraits)[61:69]
list2 <- colnames(allTraits)[60]
list2 <- colnames(allTraits)[72]
list2 <- colnames(allTraits)[65]



#run first to generate column names
traitRvalues.pcoa<-traitR_multPCs.fun(df = allTraits,speciescol = "Species",attackedcol = "GNOR",ancspecies = "A",attackedgroup = "a",traitgrp = colnames(allTraits)[3:5],test="pca",niter = 10,n=3)

out <- matrix("NA", nrow = length(list2) * length(list1), ncol = 36)
colnames(out) <-
  c(
    "traitgrp",
    "attackedspecies",
    colnames(traitRvalues.pcoa$pvals),
    colnames(traitRvalues.pcoa$distances),
    rownames(traitRvalues.pcoa$samp),"Theta1","dotprod1","Theta2","dotprod2",paste("Cent1",colnames(traitRvalues.pcoa$centroid1)[2:4]),paste("Cent2",colnames(traitRvalues.pcoa$centroid1)[2:4]),paste("Cent1A",colnames(traitRvalues.pcoa$centroid2)[2:4]),paste("Cent2A",colnames(traitRvalues.pcoa$centroid1)[2:4])
  )

dir.create("pscores")
dir.create("summ")
#PCoA


loopwrapper.fun()

