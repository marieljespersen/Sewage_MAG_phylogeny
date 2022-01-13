#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

library(tidyverse)
library(ape)
library(vegan)
library(stringr)

work_dir <- args[1] #"/home/projects/cge/data/projects/5001/Binning_vamb/clusters/$CLUSTER/city_gene_trees"
out_file <-args[2]

setwd(work_dir)

metadata <- read.delim("/home/projects/cge/data/projects/5001/Binning_vamb/metadata/gs2_metadata.txt", stringsAsFactors = F)
metadata$new_complete_name <- str_replace(metadata$new_complete_name, "-", "_")
who_regions <- read.csv("/home/projects/cge/data/projects/5001/Binning_vamb/metadata/world_bank.csv", stringsAsFactors = F)
merged_meta <- merge(metadata, who_regions, by.x = "country_alpha3", by.y = "Country.Code", sort = F)


#  function to test tree
fit_adonis.f <- function(f) {
  #cat(paste0(f, "\n"))
  
  # load tree
  tree<-read.tree(f)
  
  # get annotation
  md = merged_meta[match(tree$tip.label, merged_meta$new_complete_name),]
  colnames(md) = colnames(merged_meta)
  
  dup_regions <- unique(md$Region[duplicated(md$Region)])
  samples_to_drop <- md$new_complete_name[which(!md$Region %in% dup_regions)]
  
  # drop samples from where only one sample from the region is found
  tree <- drop.tip(tree, samples_to_drop)
  # create new sub metadata with only the samples left in the tree
  md <- md[!md$new_complete_name %in% samples_to_drop,]
  
  # check for multiple regions in tree before permanova
  if(length(unique(md$Region))>1)
  {
    # create  distance matrix from  tree
    PatristicDistMatrix<-cophenetic.phylo(tree)      
    d = as.dist(PatristicDistMatrix)
      
    # testing distance according to conttinent (and longitude and latitude) using permanova
    #fit = adonis2(d ~ Continent+longtitude_inferred+latitude_inferred, data=as.data.frame(md), permutations=9999)
    fit = adonis2(d ~ Region, data=as.data.frame(md))
    return(fit)
  }
}

#  function to get number of samples from tested trees
get_number_of_samples <- function(c) {
  f <- paste0(c, ".tree")
  
  # load tree
  tree<-read.tree(f)
  
  # get annotation
  md = merged_meta[match(tree$tip.label, merged_meta$new_complete_name),]
  colnames(md) = colnames(merged_meta)
  
  dup_regions <- unique(md$Region[duplicated(md$Region)])
  samples_to_drop <- md$new_complete_name[which(!md$Region %in% dup_regions)]
  
  # drop samples from where only one sample from the region is found
  tree <- drop.tip(tree, samples_to_drop)
  # create new sub metadata with only the samples left in the tree
  md <- md[!md$new_complete_name %in% samples_to_drop,]
  
  # check for multiple regions in tree before permanova
  if(length(unique(md$Region))>1)
  {
    return(length(tree$tip.label))
  }
}

#### Select trees with samples from multiple continents ####
# model can only be made when cluster is found in multiple continents
files <- Sys.glob("gene*.tree")
files_for_test <- c()
print("identifying trees with more than one region")
for(f in files)
{
  # load tree
  tree<-read.tree(f)
  
  #create sub matrix of metadata containing only md for samples in this tree
  sub_meta <- merged_meta[which(merged_meta$new_complete_name %in% tree$tip.label),]
  if(length(unique(sub_meta$Region))>1)
  {
    files_for_test <- c(files_for_test, f)
  }
}
print(paste("Number of trees for testing: ", length(files_for_test)))


#### run tree files trough function ####
print("testing trees...")
fit_list <- lapply(files_for_test, fit_adonis.f)
names(fit_list) <- lapply(files_for_test, function(f){strsplit(f, '[.]')[[1]][1]})
fit_list <- compact(fit_list) # exclude NULL values from the trees that could not be tested (due tolacking multiple regions or one sample per region)

# save number of tips from tested trees
samples_in_trees <- sapply(names(fit_list),get_number_of_samples)

# extract p_vals and R2 values
pvals <- matrix(NA, nrow=length(fit_list), ncol=3)
for (i in 1:length(fit_list)) {
  pvals[i,] <- c(fit_list[[i]]$`Pr(>F)`[1], fit_list[[i]]$R2[1], fit_list[[i]]$SumOfSqs[3])
}
rownames(pvals) = names(fit_list)
colnames(pvals) <- c("p_val", "R2","sumOfSqs")

# add number of samples to matrix
out_mat <- cbind(pvals, samples_in_trees)

# adjust p values
out_mat[,1] <- p.adjust(out_mat[,1], method="fdr")

# save p-vals and R2 values to file
out_df <- as.data.frame(out_mat)
out_df$p_val <- round(out_df$p_val, digits = 5)
out_df$R2 <- round(out_df$R2, digits = 5)
out_df$sumOfSqs <- round(out_df$sumOfSqs, digits = 5)
write.table(out_df, file=out_file, quote=FALSE, col.names=F)

