#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

library(tidyverse)
library(ape)
library(vegan)
library(stringr)

work_dir <- args[1] #"/home/projects/cge/data/projects/5001/Binning_vamb/clusters/C5/city_gene_trees"
out_file <- args[2] #"/home/projects/cge/data/projects/5001/Binning_vamb/clusters/C5/reg_entropy.txt"

setwd(work_dir)

metadata <- read.delim("/home/projects/cge/data/projects/5001/Binning_vamb/metadata/gs2_metadata.txt", stringsAsFactors = F)
metadata$new_complete_name <- str_replace(metadata$new_complete_name, "-", "_")
who_regions <- read.csv("/home/projects/cge/data/projects/5001/Binning_vamb/metadata/world_bank.csv", stringsAsFactors = F)
merged_meta <- merge(metadata, who_regions, by.x = "country_alpha3", by.y = "Country.Code", sort = F)


#  function to calculate regional entropy
get_reg_entropy <- function(f) {
  
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
  
  ### Udregn regional entropy
  tree_reg <- table(md$Region)
  # lav det om til fractionen
  reg_frac <- tree_reg/sum(tree_reg)

  # calculate entropy on regions
  x <- 0
  for(fraction in reg_frac)
  {
    x <- x + log(fraction^fraction)
  }
  entr = -x

  return(entr)
}


#### Select trees with samples from multiple continents ####
# model can only be made when cluster is found in multiple continents
files <- Sys.glob("gene*.tree")

#### run tree files trough function ####
print("calculating...")
entr_list <- lapply(files, get_reg_entropy)
names(entr_list) <- lapply(files, function(f){strsplit(f, '[.]')[[1]][1]})


out_dat <- t(as.data.frame(entr_list))
write.table(out_dat, file=out_file, quote=FALSE, col.names=F)
