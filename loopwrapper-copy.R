library(FD)
library(plyr)
allTraits<-read.csv("~/shinyR/TraitR_pcoa/allTraits.csv")
source("shinyR/TraitR_pcoa/traitR_standalone_pca_pcoa.R")
head(allTraits)
dim(allTraits)
June7<-colnames(allTraits[c(8:14)])
July<-colnames(allTraits[c(17:23)])
late<-colnames(allTraits[c(27:45)])
GNORlab<-colnames(allTraits[c(30,32,45:52)])
lab<-colnames(allTraits[c(30,32,38,39,46:52)]) #nophenology
phen<-colnames(allTraits[c(30,32,38,39,42:52)])

list1<-list(June7=June7,July=July,late=late,lab=lab,phen=phen)
list2<-colnames(allTraits)[64:71]

out<-matrix("NA",nrow=length(list2)*length(list1),ncol=20)
colnames(out)<-c("traitgrp","attackedspecies",colnames(traitRvalues.pcoa$pvals),
                 colnames(traitRvalues.pcoa$distances),
                 rownames(traitRvalues.pcoa$samp))

raitR.fun(df = allTraits,speciescol =
            "Species",attackedcol = "GNOR",ancspecies = "A",attackedgroup =
            "a",traitgrp = phen,test="pca",niter = 10)

z<-1
for(i in 1:length(list2))
  for (j in 1:length(list1)){
traitRvalues.pcoa<-traitR.fun(df = allTraits,speciescol =
                                "Species",attackedcol = list2[i],ancspecies = "A",attackedgroup =
                                "a",traitgrp = list1[[j]],test="pcoa",niter = 10)
out[z,1]<-paste(names(list1[i]))
out[z,2]<-paste(list2[j])
out[z,3:20]  <-as.numeric(cbind(traitRvalues.pcoa$pvals,traitRvalues.pcoa$distances,t(as.matrix(traitRvalues.pcoa$samp))))
z<-z+1
}

head(out)

###### or #####


testl<-lapply(list2,function(x) lapply(list1,function(y){traitR.fun(df = allTraits,speciescol = "Species",attackedcol = x,ancspecies = "A",attackedgroup ="a",traitgrp = y,test="pca",niter = 10)}))
tests<-sapply(list2,function(x) lapply(list1,function(y){traitR.fun(df = allTraits,speciescol = "Species",attackedcol = x,ancspecies = "A",attackedgroup ="a",traitgrp = y,test="pca",niter = 10)}))

x<-flatten(unnest(testl))
y<-unnest(testl)
str(y)
data.frame(y)

x
flatten(x)
library(plyr)
DF <- as.data.frame(
  do.call(rbind.fill.matrix,
          lapply(tests, function(l) {
            res <- unlist(l)
            t(res)
          })
  )
)
DF
renquote <- function(l) if (is.list(l)) lapply(l, renquote) else enquote(l)

x<-(lapply(unlist(renquote(testl)), eval))
y<-(lapply(unlist(renquote(x)), eval))
dim(as.matrix(y))
dim(t(DF))
DF$concert <- names(concertsList)
names(DF) <- gsub("bands.","",names(DF))
library(devtools)
install.packages("pkgutils", repos="http://R-Forge.R-project.org")
library(pkgutils)
flatten(flatten(testl))
str(test)
names(test)
test[[1]]
