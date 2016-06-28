traitR.df <- function(df,
                      speciescol,
                      attackedcol,
                      traitgrp,
                      subsetcol,
                      subsetgrp,
                      traitfac,
                      traitordfac) {
  
  if(!is.null(subsetcol)){
    df2<-df[df[[subsetcol]]==subsetgrp,c(speciescol,attackedcol,traitgrp)]}
  else{
    df2<-df[,c(speciescol,attackedcol,traitgrp)]} 
 if(!is.null(traitfac)|!is.null(traitordfac)){
  if(!is.null(traitfac))
 df2[traitfac]<-data.frame(lapply(df2[traitfac],as.factor))

 if(!is.null(traitordfac)){
   df2[traitordfac]<-data.frame(lapply(df2[traitordfac],as.ordered))
  }
 }
  return(df2)}
  

pcscores <-
  function(pcdf,
           speciesCol,
           attackedCol,traitfac=NULL,traitordfac=NULL,correction="none")
  {
    if(!is.null(traitfac)| !is.null(traitordfac)){
      pcdf<-na.omit(pcdf)
      pcdf2<-pcdf[,3:dim(pcdf)[2]]
      pcg <- gowdis(pcdf2)
      pcComp <<- pcoa(pcg, correction = correction)$vectors[, 1:2]
      colnames(pcComp) <- c("Comp.1", "Comp.2")
    }
    else{
      pcdf<-na.omit(pcdf)
      pcdf2<-pcdf[,3:dim(pcdf)[2]]
      pcComp <-
        princomp(
          pcdf2,
          scores = TRUE,
          cor = T
        )$scores[,1:2]}
    Species.Attacked <-
      interaction(pcdf[[speciesCol]], pcdf[[attackedCol]])
    pcsc <- data.frame(cbind(pcdf, pcComp, Species.Attacked))
    return(pcsc)
  }



centroids <- function(df, speciescol, attackedcol) {
  PC_species <-
    data.frame(ddply(
      df,
      speciescol,
      summarize,
      meanPC1 = mean(Comp.1),
      meanPC2 = mean(Comp.2)
    ))
  PC_nov_att <-
    data.frame(ddply(
      df,
      c(speciescol, attackedcol),
      summarize,
      meanPC1 = mean(Comp.1),
      meanPC2 = mean(Comp.2)
    ))
  PC_nov_att[, "Species.Attacked"] <-
    interaction(PC_nov_att[, 1], PC_nov_att[, 2])
  y <- list(PC_species, PC_nov_att)
  return(y)
}

dist.fun <- function(df_pca, df_centroid, what, group) {
  Cent.O <-
    na.omit(df_centroid[df_centroid[[what]] == group, c("meanPC1", "meanPC2")])
  pcapoints <-
    df_pca[df_pca[[what]] == group, c("Comp.1", "Comp.2")]
  names(Cent.O) <- names(pcapoints)
  mat <- rbind("centroid" = Cent.O, pcapoints)
  dist.mat <- dist(mat, method = "euclidian")
  dist.mat <- as.matrix(dist.mat)
  rownames(dist.mat) <- mat[, 1]
  x <- as.data.frame(dist.mat)[, 1]
  return(sum(x))
}

traitR.perm <-
  function(pcadf, 
           speciescol,
           attackedcol,
           ancspecies,
           attackedgroup,
           niter)
  {
    novelsp <-
      levels(pcadf[[speciescol]])[levels(pcadf[[speciescol]]) != ancspecies]
    novelsp <- novelsp[novelsp != ""]
    unattacked <-
      levels(pcadf[[attackedcol]])[levels(pcadf[[attackedcol]]) != attackedgroup]
    cen <-
      centroids(df = pcadf,
                speciescol = speciescol,
                attackedcol = attackedcol)
    anc.G <-
      dist.fun(
        df_pca = pcadf,
        df_centroid = cen[[2]],
        what = "Species.Attacked",
        group = paste(ancspecies, ".", attackedgroup, sep = "")
      )
    novel.G <-
      dist.fun(
        df_pca = pcadf,
        df_centroid = cen[[2]],
        what = "Species.Attacked",
        group = paste(novelsp, ".", attackedgroup, sep
                      =
                        "")
      )
    anc.all <-
      dist.fun(
        df_pca = pcadf,
        df_centroid = cen[[1]],
        what = "Species",
        group = ancspecies
      )
    novel.all <-
      dist.fun(
        df_pca = pcadf,
        df_centroid = cen[[1]],
        what = "Species",
        group = novelsp
      )
    cen_real <- cen
    traitR.mat <- matrix(data = NA,
                         nrow = 1,
                         ncol = 5)
    colnames(traitR.mat) <-
      c(
        "Loc-anc-ancAtt",
        "Loc-nov-novAtt",
        "Loc-ancAtt-novAtt",
        "Size-anc-att",
        "Size-nov-att"
      )
    
    # calculate distances between centroids
    anc.all.cen <- cen[[1]][cen[[1]][speciescol] == ancspecies][2:3]
    anc.att.cen <-
      cen[[2]][cen[[2]]["Species.Attacked"] == paste(ancspecies, ".", attackedgroup, sep =
                                                       "")][3:4]
    nov.all.cen <- cen[[1]][cen[[1]][speciescol] == novelsp][2:3]
    nov.att.cen <-
      cen[[2]][cen[[2]]["Species.Attacked"] == paste(novelsp, ".", attackedgroup, sep =
                                                       "")][3:4]
    traitR.mat[1, 1] <- dist(rbind(anc.all.cen, anc.att.cen))
    traitR.mat[1, 2] <- dist(rbind(nov.all.cen, nov.att.cen))
    traitR.mat[1, 3] <- dist(rbind(anc.att.cen, nov.att.cen))
    traitR.mat[1, 4:5] <- c(anc.G, novel.G)
    
    traitR.permutations <- matrix(data = NA,
                                  nrow = niter,
                                  ncol = 5)
    colnames(traitR.permutations) <-
      c(
        "Loc-anc-ancAtt",
        "Loc-nov-novAtt",
        "Loc-ancAtt-novAtt",
        "Size-anc-att",
        "Size-nov-att"
      )
    
    withProgress(message = 'Running permutation', value = 0, {
      for (i in 1:niter) {
        incProgress(1 / niter, detail = paste("#:", i))
        ancestralNewdf <- pcadf[pcadf[speciescol] == ancspecies,]
        ancestralNewdf$"Newattack" <-
          sample(ancestralNewdf[[attackedcol]],
                 replace = F,
                 size = dim(ancestralNewdf)[1])
        novelNewdf <- pcadf[pcadf[speciescol] == novelsp,]
        novelNewdf$"Newattack" <-
          sample(novelNewdf[[attackedcol]],
                 replace = F,
                 size = dim(novelNewdf)[1])
        newpca <- rbind(ancestralNewdf, novelNewdf)
        newpca["Species.Attacked"] <-
          interaction(newpca[[speciescol]], newpca[["Newattack"]])
        
        cen <-
          centroids(df = newpca,
                    speciescol = speciescol,
                    attackedcol = "Newattack")
        anc.G <-
          dist.fun(
            df_pca = newpca,
            df_centroid = cen[[2]],
            what = "Species.Attacked",
            group = paste(ancspecies, ".", attackedgroup, sep = "")
          )
        novel.G <-
          dist.fun(
            df_pca = newpca,
            df_centroid = cen[[2]],
            what = "Species.Attacked",
            group = paste(novelsp, ".", attackedgroup, sep = "")
          )
        anc.all <-
          dist.fun(
            df_pca = newpca,
            df_centroid = cen[[1]],
            what = "Species",
            group = ancspecies
          )
        novel.all <-
          dist.fun(
            df_pca = newpca,
            df_centroid = cen[[1]],
            what = "Species",
            group = novelsp
          )
        
        
        # calculate distances between centroids
        anc.all.cen <-
          cen[[1]][cen[[1]][speciescol] == ancspecies][2:3]
        anc.att.cen <-
          cen[[2]][cen[[2]]["Species.Attacked"] == paste(ancspecies, ".", attackedgroup, sep =
                                                           "")][3:4]
        nov.all.cen <-
          cen[[1]][cen[[1]][speciescol] == novelsp][2:3]
        nov.att.cen <-
          cen[[2]][cen[[2]]["Species.Attacked"] == paste(novelsp, ".", attackedgroup, sep =
                                                           "")][3:4]
        traitR.permutations[i, 1] <-
          dist(rbind(anc.all.cen, anc.att.cen))
        traitR.permutations[i, 2] <-
          dist(rbind(nov.all.cen, nov.att.cen))
        traitR.permutations[i, 3] <-
          dist(rbind(anc.att.cen, nov.att.cen))
        traitR.permutations[i, 4:5] <- c(anc.G, novel.G)
        
      }
    })
    traitR.pvals <- matrix(data = NA,
                           nrow = 1,
                           ncol = 5)
    colnames(traitR.pvals) <-
      c(
        "Location-Ancestral-AncestralAttacked",
        "Location-Novel-NovelAttacked",
        "Location-AncestralAttacked-NovelAttacked",
        "Spread-AncestralAttacked",
        "Spread-NovelAttacked"
      )
    traitR.pvals[1, 1] <-
      (sum(traitR.permutations[, 1] >= traitR.mat[1, 1]) + 1) / (niter + 1)
    traitR.pvals[1, 2] <-
      (sum(traitR.permutations[, 2] >= traitR.mat[1, 2]) + 1) / (niter + 1)
    test_big <-
      (sum(traitR.permutations[, 3] >= traitR.mat[1, 3]) + 1) / (niter + 1)
    test_small <-
      (sum(traitR.permutations[, 3] <= traitR.mat[1, 3]) + 1) / (niter + 1)
    traitR.pvals[1, 3] <- 2 * min(test_big, test_small)
    traitR.pvals[1, 4] <-
      (sum(traitR.permutations[, 4] <= traitR.mat[1, 4]) + 1) / (niter + 1)
    traitR.pvals[1, 5] <-
      (sum(traitR.permutations[, 5] <= traitR.mat[1, 5]) + 1) / (niter + 1)
    return(list(
      traitR.mat,
      traitR.permutations,
      traitR.pvals,
      cen_real[[1]],
      cen_real[[2]],
      pcadf,
      table(pcadf$Species.Attacked)
    ))
    
  }
