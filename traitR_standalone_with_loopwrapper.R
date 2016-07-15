
#function for the histograms for the shiny version
histfun<-function(bsobj,pvalobj,dataobj,coln){
  hist(bsobj[,coln], xlab = paste(names(bsobj[coln])), main= paste(names(pvalobj[coln])))
  abline(v=dataobj[,coln], col="red")
}

#functions for pca/pcoa. 
#Extract the scores and the loadings
#
pcascores <- function(df, speciesCol,attackedCol,traitGrp,n)
{      
  if(is.null(n)) n<-length(traitGrp)-1
  pcnona <- na.omit(subset(df,select=c(speciesCol,attackedCol,traitGrp)))
  pcdf<-subset(pcnona, select=traitGrp)
  pcComp <-princomp(pcdf, scores=TRUE,cor=T)
  Species.Attacked<-interaction(pcnona[[speciesCol]],pcnona[[attackedCol]])
  pcsc<-data.frame(cbind(pcnona[,c(speciesCol,attackedCol)],Species.Attacked,pcComp$scores[,1:n]))
  colnames(pcsc)[4:dim(pcsc)[2]]<-paste("Comp.", 1:n, sep="")
  return(list(pcsc=pcsc,load=pcComp$loadings))
}

pcoascores <- function(df, speciesCol,attackedCol,traitGrp,correction="none",n)
{      
  if(is.null(n)) n<-length(traitGrp)-1
  pcnona <- na.omit(subset(df,select=c(speciesCol,attackedCol,traitGrp)))
  pcdf<-subset(pcnona, select=traitGrp)
  pcg<-gowdis(pcdf)
  pcComp <-pcoa(pcg, correction=correction)
  Species.Attacked<-interaction(pcnona[[speciesCol]],pcnona[[attackedCol]])
  pcsc<-data.frame(cbind(pcnona[,c(speciesCol,attackedCol)],Species.Attacked,pcComp$vectors[,1:n]))
  colnames(pcsc)[4:dim(pcsc)[2]]<-paste("Comp.", 1:n, sep="")
  return(list(pcsc=pcsc,load=NULL))
}



#extract the centroids.
#change this so it does centroids for all PC axes.
#dplyr will do summarise if numeric
centroids<-function(df, speciescol, attackedcol){
  df[,"Species.Attacked"]<-interaction(df[[speciescol]],df[[attackedcol]])
  PC_species<-data.frame(ddply(df, speciescol,numcolwise(mean)))
  PC_nov_att<-data.frame(ddply(df, "Species.Attacked", numcolwise(mean)))
  y<-list(PC_species,PC_nov_att)
  return(y)
  }

#work out the euclidian distances
#dist has probelms if there are categorical variables in the df
dist.fun<-function(df_pca, df_centroid,what,group){
  Cent.O<-na.omit(df_centroid[df_centroid[[what]]==group,])
  pcapoints<-df_pca[df_pca[[what]]==group,]
  mat<-rbind("centroid"=Cent.O[,grep("Comp",colnames(Cent.O))],pcapoints[,grep("Comp",colnames(pcapoints))])
  dist.mat<-dist(mat,method = "euclidian")
  dist.mat<-as.matrix(dist.mat)
  rownames(dist.mat) <-mat[,1]
  x<-as.data.frame(dist.mat)[,1]
  return(sum(x))
}

#work out angle between centroids
angle <- function(Orig,Cord1,Cord2){
  x <-Cord1-Orig
  y<-Cord2-Orig
  dot.prod <- x%*%t(y)
  norm.x <- norm(x,type="2")
  norm.y <- norm(y,type="2")
  theta <- acos(dot.prod / (norm.x * norm.y))
  return(list(theta=as.numeric(theta)*180/pi,dotprod=dot.prod))
}


  
traitR_multPCs.fun<-function(df,speciescol,attackedcol,ancspecies,attackedgroup,traitgrp,niter,test="pca",correction=NULL,n=NULL)
{
  if (!require("pacman")) install.packages("pacman")
  pacman::p_load(FD,ape, plyr)
  #run the pca/pcoa
  df2<- na.omit(subset(df,select=c(speciescol,attackedcol)))
  novelsp<-levels(df2[[speciescol]])[levels(df2[[speciescol]])!=ancspecies]
  novelsp<-novelsp[novelsp!=""]
  if(test=="pca"){
  pca_all<-pcascores(df = df,speciesCol = speciescol,attackedCol = attackedcol,traitGrp = traitgrp,n=n)}
  if(test=="pcoa"){
    pca_all<-pcoascores(df = df,speciesCol = speciescol,attackedCol = attackedcol,traitGrp = traitgrp,correction=correction,n=n)}
  pca<-pca_all$pcsc
  #calculate the centroids and distances
  cen<-centroids(df = pca,speciescol = speciescol,attackedcol = attackedcol)
  anc.G<-dist.fun(df_pca = pca,df_centroid = cen[[2]],what="Species.Attacked",group = paste(ancspecies,".",attackedgroup,sep=""))
  novel.G<-dist.fun(df_pca = pca,df_centroid = cen[[2]],what="Species.Attacked",group = paste(novelsp,".",attackedgroup,sep=""))
  anc.all<-dist.fun(df_pca = pca,df_centroid = cen[[1]],what="Species",group = ancspecies)
  novel.all<-dist.fun(df_pca = pca,df_centroid = cen[[1]],what="Species",group = novelsp)
  cen_real<-cen
  traitR.mat<-matrix(data=NA,nrow=1,ncol=9)
  colnames(traitR.mat)<-c("Loc-anc-ancAtt","Loc-nov-novAtt","Loc-ancAtt-novAtt", "Size-anc-att","Size-nov-att","PC1_A","PC2_A","PC1_G","PC2_G")
 
  # distances output to traitR.mat
  anc.all.cen<-cen[[1]][cen[[1]][speciescol]==ancspecies][2:dim(cen[[1]])[2]]
  anc.att.cen<-cen[[2]][cen[[2]]["Species.Attacked"]==paste(ancspecies,".",attackedgroup,sep="")][2:dim(cen[[2]])[2]]
  nov.all.cen<-cen[[1]][cen[[1]][speciescol]==novelsp][2:dim(cen[[1]])[2]]
  nov.att.cen<-cen[[2]][cen[[2]]["Species.Attacked"]==paste(novelsp,".",attackedgroup,sep="")][2:dim(cen[[2]])[2]]
  traitR.mat[1,1]<-dist(rbind(anc.all.cen,anc.att.cen))
  traitR.mat[1,2]<-dist(rbind(nov.all.cen,nov.att.cen))
  traitR.mat[1,3]<-dist(rbind(anc.att.cen,nov.att.cen))
  traitR.mat[1,4:5]<-c(anc.G,novel.G)
  #traitR.mat[6:7]<-unlist((cen[[1]][1,]-cen[[2]][1,]))
  #traitR.mat[8:9]<-unlist(cen[[1]][2,]-cen[[2]][3,])
  
  traitR.permutations<-matrix(data=NA,nrow=niter,ncol=5)
  colnames(traitR.permutations)<-c("Loc-anc-ancAtt","Loc-nov-novAtt","Loc-ancAtt-novAtt", "Size-anc-att","Size-nov-att")
 
  pb <- txtProgressBar(min = 0, max = niter, style = 3) 
for (i in 1:niter){
  setTxtProgressBar(pb, i)    
  ancestralNewdf <- pca[pca[speciescol]==ancspecies,]
  ancestralNewdf$"Newattack"<-sample(ancestralNewdf[[attackedcol]], replace=F, size=dim(ancestralNewdf)[1])
  novelNewdf <- pca[pca[speciescol]==novelsp,]
  novelNewdf$"Newattack"<-sample(novelNewdf[[attackedcol]], replace=F, size=dim(novelNewdf)[1])
  newpca<-rbind(ancestralNewdf,novelNewdf)
  newpca["Species.Attacked"]<-interaction(newpca[[speciescol]],newpca[["Newattack"]])
 
  cen<-centroids(df = newpca,speciescol = speciescol,attackedcol = "Newattack")
  anc.G<-dist.fun(df_pca = newpca,df_centroid = cen[[2]],what="Species.Attacked",group = paste(ancspecies,".",attackedgroup,sep=""))
  novel.G<-dist.fun(df_pca = newpca,df_centroid = cen[[2]],what="Species.Attacked",group = paste(novelsp,".",attackedgroup,sep=""))
  anc.all<-dist.fun(df_pca = newpca,df_centroid = cen[[1]],what="Species",group = ancspecies)
  novel.all<-dist.fun(df_pca = newpca,df_centroid = cen[[1]],what="Species",group = novelsp)


  # distances output to traitR.mat
  anc.all.cen<-cen[[1]][cen[[1]][speciescol]==ancspecies][2:dim(cen[[1]])[2]]
  anc.att.cen<-cen[[2]][cen[[2]]["Species.Attacked"]==paste(ancspecies,".",attackedgroup,sep="")][2:dim(cen[[2]])[2]]
  nov.all.cen<-cen[[1]][cen[[1]][speciescol]==novelsp][2:dim(cen[[1]])[2]]
  nov.att.cen<-cen[[2]][cen[[2]]["Species.Attacked"]==paste(novelsp,".",attackedgroup,sep="")][2:dim(cen[[2]])[2]]
   traitR.permutations[i,1]<-dist(rbind(anc.all.cen,anc.att.cen))
  traitR.permutations[i,2]<-dist(rbind(nov.all.cen,nov.att.cen))
  traitR.permutations[i,3]<-dist(rbind(anc.att.cen,nov.att.cen))
  traitR.permutations[i,4:5]<-c(anc.G,novel.G)

} 

traitR.pvals<-matrix(data=NA,nrow=1,ncol=5)
  colnames(traitR.pvals)<-c("Location-Ancestral-AncestralAttacked","Location-Novel-NovelAttacked","Location-AncestralAttacked-NovelAttacked", "Spread-AncestralAttacked","Spread-NovelAttacked")
  traitR.pvals[1,1]<-(sum(traitR.permutations[,1] >= traitR.mat[1,1])+1) /(niter+1)
  traitR.pvals[1,2]<-(sum(traitR.permutations[,2] >= traitR.mat[1,2])+1) /(niter+1)
  test_big<-(sum(traitR.permutations[,3] >= traitR.mat[1,3])+1)/ (niter+1)
  test_small<-(sum(traitR.permutations[,3] <= traitR.mat[1,3])+1)/ (niter+1)
  traitR.pvals[1,3]<-2*min(test_big,test_small)
  traitR.pvals[1,4]<-(sum(traitR.permutations[,4] <= traitR.mat[1,4])+1) / (niter+1)
  traitR.pvals[1,5]<-(sum(traitR.permutations[,5] <= traitR.mat[1,5])+1) / (niter+1)
  
  Theta1Origin<-as.matrix(cen_real[[1]][1,2:dim(cen_real[[1]])[2]])
  Theta1Cord1<-as.matrix(cen_real[[1]][2,2:dim(cen_real[[1]])[2]])
  Theta1Cord2<-as.matrix(cen_real[[2]][1,2:dim(cen_real[[2]])[2]])
  
  Theta2Origin<-as.matrix(cen_real[[1]][2,2:dim(cen_real[[1]])[2]])
  Theta2Cord1<-as.matrix(cen_real[[1]][1,2:dim(cen_real[[1]])[2]])
  Theta2Cord2<-as.matrix(cen_real[[2]][2,2:dim(cen_real[[2]])[2]])
  
  theta1<-angle(Theta1Origin,Theta1Cord1,Theta1Cord2)
  theta2<-angle(Theta2Origin,Theta2Cord1,Theta2Cord2)
  
  
  print(list(distances=traitR.mat,pvals=traitR.pvals,samp=table(pca_all$pcsc["Species.Attacked"])))
  return(list(distances=traitR.mat,permutations=traitR.permutations,pvals=traitR.pvals, centroid1=cen_real[[1]],centroid2=cen_real[[2]],pcscores=pca_all$pcsc,loadings=pca_all$load,samp=table(pca_all$pcs["Species.Attacked"]),theta1=theta1,theta2=theta2))

  close(pb)
}


loopwrapper.fun<-function(      df = allTraits,
                                speciescol =
                                  "Species",
                                attackedcol = list2,
                                ancspecies = "A",
                                attackedgroup =
                                  "a",
                                traitgrp = list1,
                                test = "pcoa", niter = 10, correction = "none",n=3, filename="TraitRoutput"){
  #run initial to get column names
  traitRvalues.pcoa <- traitR_multPCs.fun(
    df = allTraits,
    speciescol =
      "Species",
    attackedcol = list2[1],
    ancspecies = "A",
    attackedgroup =
      "a",
    traitgrp = list1[[1]],
    test = "pcoa", niter = 10, correction = "none",n=2)
  
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
  z <- 1
  for (i in 1:length(list2))
    for (j in 1:length(list1)) {
      traitRvalues.pcoa <- traitR_multPCs.fun(
        df = allTraits,
        speciescol =
          "Species",
        attackedcol = list2[i],
        ancspecies = "A",
        attackedgroup =
          "a",
        traitgrp = list1[[j]],
        test = "pcoa", niter = 10, correction = "none",n=2)
      
      out[z, 1] <- paste(names(list1[j]))
      out[z, 2] <- paste(list2[i])
      out[z, 3:20]  <-
        as.numeric(cbind(
          traitRvalues.pcoa$pvals,
          traitRvalues.pcoa$distances,
          t(as.matrix(traitRvalues.pcoa$samp))
        ))
      out[z,21:22]<-as.numeric(traitRvalues.pcoa$theta1)
      out[z,23:24]<-as.numeric(traitRvalues.pcoa$theta2)
      if(dim(traitRvalues.pcoa$centroid1)[2]>3){
        out[z,25:27]<-as.numeric(traitRvalues.pcoa$centroid1[1,2:4])
        out[z,28:30]<-as.numeric(traitRvalues.pcoa$centroid1[2,2:4])
        out[z,31:33]<-as.numeric(traitRvalues.pcoa$centroid2[1,2:4])
        out[z,34:36]<-as.numeric(traitRvalues.pcoa$centroid2[2,2:4])
      }
      else
      {
        out[z,25:26]<-as.numeric(traitRvalues.pcoa$centroid1[1,2:3])
        out[z,28:29]<-as.numeric(traitRvalues.pcoa$centroid1[2,2:3])
        out[z,31:32]<-as.numeric(traitRvalues.pcoa$centroid2[1,2:3])
        out[z,34:35]<-as.numeric(traitRvalues.pcoa$centroid2[2,2:3])
      }
      write.csv(
        traitRvalues.pcoa$pcscores,
        file = paste0("pscores/", names(list1[j]), "_", list2[i], ".csv"),
        quote = F,
        row.names = F
      )
      sink("summaryfiles.txt",append=T)
      print(paste0(names(list1[j]), "_", list2[i], ".txt"))
      print(traitRvalues.pcoa$summ)
      sink()
      z <<- z + 1
    }
  
  write.csv(out,paste0(filename,".csv"),quote=F,row.names=F)
  return(out)
}
