histfun<-function(bsobj,pvalobj,dataobj,coln){
  hist(bsobj[,coln], xlab = paste(names(bsobj[coln])), main= paste(names(pvalobj[coln])))
  abline(v=dataobj[,coln], col="red")
}

centroids<-function(df, speciescol, attackedcol){
  PC_species<-data.frame(ddply(df, speciescol,summarize, meanPC1 = mean(Comp.1), meanPC2=mean(Comp.2)))
  PC_nov_att<-data.frame(ddply(df, c(speciescol,attackedcol), summarize, meanPC1 = mean(Comp.1), meanPC2=mean(Comp.2)) )
  PC_nov_att[,"Species.Attacked"]<-interaction(PC_nov_att[,1],PC_nov_att[,2])
  y<-list(PC_species,PC_nov_att)
  return(y)
  }

dist.fun<-function(df_pca, df_centroid,what,group){
  Cent.O<-na.omit(df_centroid[df_centroid[[what]]==group,c("meanPC1","meanPC2")])
  pcapoints<-df_pca[df_pca[[what]]==group,c("Comp.1","Comp.2")]
  names(Cent.O)<-names(pcapoints)
  mat<-rbind("centroid"=Cent.O,pcapoints)
  dist.mat<-dist(mat,method = "euclidian")
  dist.mat<-as.matrix(dist.mat)
  rownames(dist.mat) <-mat[,1]
  x<-as.data.frame(dist.mat)[,1]
  return(sum(x))
}

pcascores <- function(df, speciesCol,attackedCol,traitGrp)
{      
  pcnona <- na.omit(subset(df,select=c(speciesCol,attackedCol,traitGrp)))
  pcdf<-subset(pcnona, select=traitGrp)
  pcComp <-princomp(pcdf, scores=TRUE,cor=T)
  Species.Attacked<-interaction(pcnona[[speciesCol]],pcnona[[attackedCol]])
  pcsc<-data.frame(cbind(pcnona,pcComp$scores[,1:2],Species.Attacked))
  return(list(pcsc=pcsc,load=pcComp$loadings))
}

pcoascores <- function(df, speciesCol,attackedCol,traitGrp,correction)
  {      
  pcnona <- na.omit(subset(df,select=c(speciesCol,attackedCol,traitGrp)))
  pcdf<-subset(pcnona, select=traitGrp)
  pcg<-gowdis(pcdf)
  pcComp <-pcoa(pcg, correction=correction)
  pcComp2 <-pcComp$vectors[,1:2]
  colnames(pcComp2)<-c("Comp.1","Comp.2")
  Species.Attacked<-interaction(pcnona[[speciesCol]],pcnona[[attackedCol]])
    pcsc<-data.frame(cbind(pcnona,pcComp2,Species.Attacked))
        return(list(pcsc=pcsc,load=NULL))
      }
    



  
traitR.fun<-function(df,speciescol,attackedcol,ancspecies,attackedgroup,traitgrp,niter,test="pca",correction=NULL)
{
  df2<- na.omit(subset(df,select=c(speciescol,attackedcol)))
  novelsp<-levels(df2[[speciescol]])[levels(df2[[speciescol]])!=ancspecies]
  novelsp<-novelsp[novelsp!=""]
  unattacked<-levels(df2[[attackedcol]])[levels(df2[[attackedcol]])!=attackedgroup]
  if(test=="pca"){
  pca_all<-pcascores(df = df,speciesCol = speciescol,attackedCol = attackedcol,traitGrp = traitgrp)}
  if(test=="pcoa"){
    pca_all<-pcoascores(df = df,speciesCol = speciescol,attackedCol = attackedcol,traitGrp = traitgrp,correction=correction)}
  pca<-pca_all$pcsc
  cen<-centroids(df = pca,speciescol = speciescol,attackedcol = attackedcol)
  anc.G<-dist.fun(df_pca = pca,df_centroid = cen[[2]],what="Species.Attacked",group = paste(ancspecies,".",attackedgroup,sep=""))
  novel.G<-dist.fun(df_pca = pca,df_centroid = cen[[2]],what="Species.Attacked",group = paste(novelsp,".",attackedgroup,sep=""))
  anc.all<-dist.fun(df_pca = pca,df_centroid = cen[[1]],what="Species",group = ancspecies)
  novel.all<-dist.fun(df_pca = pca,df_centroid = cen[[1]],what="Species",group = novelsp)
  cen_real<-cen
  traitR.mat<-matrix(data=NA,nrow=1,ncol=9)
  colnames(traitR.mat)<-c("Loc-anc-ancAtt","Loc-nov-novAtt","Loc-ancAtt-novAtt", "Size-anc-att","Size-nov-att","PC1_A","PC2_A","PC1_G","PC2_G")
  #colnames(traitR.mat<-c("Loc-anc-ancAtt","Loc-nov-novAtt","Loc-ancAtt-novAtt", "Size-anc-att","Size-nov-att"))
  # calculate distances between centroids
  anc.all.cen<-cen[[1]][cen[[1]][speciescol]==ancspecies][2:3]
  anc.att.cen<-cen[[2]][cen[[2]]["Species.Attacked"]==paste(ancspecies,".",attackedgroup,sep="")][3:4]
  nov.all.cen<-cen[[1]][cen[[1]][speciescol]==novelsp][2:3]
  nov.att.cen<-cen[[2]][cen[[2]]["Species.Attacked"]==paste(novelsp,".",attackedgroup,sep="")][3:4]
  traitR.mat[1,1]<-dist(rbind(anc.all.cen,anc.att.cen))
  traitR.mat[1,2]<-dist(rbind(nov.all.cen,nov.att.cen))
  traitR.mat[1,3]<-dist(rbind(anc.att.cen,nov.att.cen))
  traitR.mat[1,4:5]<-c(anc.G,novel.G)
  traitR.mat[6:7]<-unlist((cen[[1]][1,2:3]-cen[[2]][1,3:4]))
  traitR.mat[8:9]<-unlist(cen[[1]][2,2:3]-cen[[2]][3,3:4])
  
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


  # calculate distances between centroids
  anc.all.cen<-cen[[1]][cen[[1]][speciescol]==ancspecies][2:3]
  anc.att.cen<-cen[[2]][cen[[2]]["Species.Attacked"]==paste(ancspecies,".",attackedgroup,sep="")][3:4]
  nov.all.cen<-cen[[1]][cen[[1]][speciescol]==novelsp][2:3]
  nov.att.cen<-cen[[2]][cen[[2]]["Species.Attacked"]==paste(novelsp,".",attackedgroup,sep="")][3:4]
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
  print(list(distances=traitR.mat,pvals=traitR.pvals,samp=table(pca_all$pcsc["Species.Attacked"])))
  return(list(distances=traitR.mat,permutations=traitR.permutations,pvals=traitR.pvals, centroid1=cen_real[[1]],centroid2=cen_real[[2]],pcscores=pca_all$pcsc,loadings=pca_all$load,samp=table(pca_all$pcs["Species.Attacked"])))

  close(pb)
}

