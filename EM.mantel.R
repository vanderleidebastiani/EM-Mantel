## R function for estimating phylogenetic signal in traits using EM-Mantel test.

# Debastiani, V.J., Duarte, L.d. Evolutionary Models and Phylogenetic Signal Assessment via Mantel Test. Evol Biol 44, 135–143 (2017). 
# https://doi.org/10.1007/s11692-016-9396-1


## Required packages

require(vegan)
require(FD)

## Description

# Calculate the phylogenetic signal based in Mantel test, incorporating the Brownian motion evolutionary model.

## Arguments

# tree = Phylogenetic tree, as phylo object.
# traits = Matrix or data frame containing the traits data, with traits as columns and species as rows. Traits can be numeric, ordered, or factor, as the required by the gowdis function. (Symmetric or asymmetric binary variables should be numeric and only contain 0 and 1. character variables will be converted to factor).
# runs = Number of permutations in assessing significance.
# euclidean = Logical argument (TRUE or FALSE) to specify if use transformation of  euclidean properties in pairwise trait dissimilarities (Default euclidean = TRUE).
# sqrtPhylo = Logical argument (TRUE or FALSE) to specify if use square root transformation of phylogenetic distance (Default sqrtPhylo = FALSE).
# checkdata = Logical argument (TRUE or FALSE) to check if species sequence in the trait data follows the same order as in phylogenetic tree (Default checkdata = TRUE).
# ... = Parameters for gowdis function.

## Value

# perm.NULL = A vector of permuted Mantel statistic under the null hypothesis that the distances between both matrices are not related.
# perm.BM = A vector of permuted Mantel statistic under the null hypothesis that the traits evolve under a Brownian motion evolutionary model.
# r.Mantel = The Mantel observed statistic.
# p.NULL = The p value from no phylogenetic structure (standard p value in Mantel test).
# p.BM = The p value under simulation of traits from Brownian phylogenetic structure.


EM.mantel<-function(tree, traits, runs = 999, euclidean= TRUE, sqrtPhylo=FALSE, checkdata = TRUE, ...){
  phylo.dist<-cophenetic(tree)
  if(sqrtPhylo){
    phylo.dist<-sqrt(phylo.dist)
  }
  if(checkdata){
    if(is.null(tree$tip.label)){
      stop("\n Error in tip labels of tree\n")
    }
    if(is.null(rownames(traits))){
      stop("\n Error in row names of traits\n")
    }
    match.names <- match(rownames(traits),rownames(phylo.dist))
    if(sum(is.na(match.names)) > 0){
      stop("\n There are species from traits data that are not on phylogenetic tree\n")
    }
    phylo.dist <- phylo.dist[match.names, match.names]
  }
  if(length(tree$tip.label) > dim(traits)[1]){
    warning("Tree have more species that species in traits data")
  }
  if(dim(phylo.dist)[1] != dim(traits)[1] & checkdata == FALSE){
    stop("\n Different number of species in tree and in traits data, use checkdata = TRUE\n")
  }
  gow.dist<-gowdis(traits, ...)
  if(euclidean){
    gow.sim<-1-gow.dist
    gow.dist<-sqrt(1-gow.sim)
  }
  traits.attr<-attr(gow.dist, "Types", exact = TRUE)
  res.mantel<-mantel(phylo.dist,gow.dist,permutations=runs)
  res.BM<-matrix(NA,runs,1)
  for(k in 1:runs){
    traits_sim<-matrix(NA,length(tree$tip.label),dim(traits)[2])
    rownames(traits_sim)<-tree$tip.label
    for(i in 1:dim(traits)[2]){
      traits_sim[,i]<-rTraitCont(tree,model="BM")
    }
    traits_sim<-decostand(traits_sim,method="standardize",MARGIN=2)
    traits_sim<-as.data.frame(traits_sim)	
    for(i in 1:dim(traits)[2]){
      if(traits.attr[i] == "B" | traits.attr[i] == "A"){
        probs<-sum(traits[,i])/dim(traits)[1]
        threshold<-quantile(traits_sim[,i],probs=1-probs)
        traits_sim[,i]<-ifelse(traits_sim[,i]>=threshold,1,0)
      }
      if(traits.attr[i] == "N" | traits.attr[i] == "O"){
        n.levels<-length(levels(traits[,i]))
        traits.levels<-levels(traits[,i])
        probs<-cumsum(table(traits[,i]))/sum(table(traits[,i]))
        probs<-probs[1:(n.levels-1)]
        threshold<-quantile(traits_sim[,i],probs=probs)
        threshold<-c(min(traits_sim[,i]),threshold,max(traits_sim[,i]))
        temp<-matrix(NA,length(traits_sim[,i]),1)
        for(j in 1:n.levels){
          if(j < n.levels){
            temp[1:length(traits_sim[,i]),1]<-ifelse(traits_sim[,i]>=threshold[j] & traits_sim[,i]<threshold[j+1], traits.levels[j],temp)
          }
          if(j == n.levels){
            temp[1:length(traits_sim[,i]),1]<-ifelse(traits_sim[,i]>=threshold[j] & traits_sim[,i]<=threshold[j+1], traits.levels[j],temp)
          }
        }
        traits_sim[,i]<-as.factor(temp)
        if(traits.attr[i] == "O"){
          traits_sim[,i]<-ordered(temp,levels=levels(traits[,i]))
        }
      }
    }
    if(checkdata == TRUE){
      match.names <- match(rownames(traits),rownames(traits_sim))
      traits_sim<-traits_sim[match.names,,drop=FALSE]
    }
    gow.dist.BM<-gowdis(traits_sim, ...)
    if(euclidean){
      gow.sim.BM<-1-gow.dist.BM
      gow.dist.BM<-sqrt(1-gow.sim.BM)
    }
    res.mantel.BM<-mantel(phylo.dist,gow.dist.BM,permutations=0)
    res.BM[k,1]<-res.mantel.BM$statistic
  }
  p.BM<-(sum(ifelse(res.BM[,1]>=res.mantel$statistic,1,0))+1)/(runs+1)
  p.NULL<-res.mantel$signif
  r.Mantel<-res.mantel$statistic
  RES<-list(perm.NULL=res.mantel$perm,perm.BM=res.BM[,1],r.Mantel=r.Mantel,p.NULL=p.NULL,p.BM=p.BM)
  return(RES)
}

## Examples

require(geiger)
tree<-sim.bdtree(b=0.1,d=0,stop="taxa",n=100,extinct=FALSE)
trait<-matrix(rTraitCont(compute.brlen(tree,power=5),model="BM"),100,1)
rownames(trait)<-tree$tip.label
EM.mantel(tree,trait,runs=99)