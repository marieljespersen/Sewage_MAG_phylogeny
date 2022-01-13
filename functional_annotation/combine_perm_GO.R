#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

library(dplyr)
library(ggpubr)
library(tidyr)

dir <- args[1]
cluster <- args[2]

#### Testing on Cellular Compartment ####
#loading file with permanova results from all genetrees
perm_file <- paste0(dir, "/", cluster,"/city_gene_trees/permanova.citygenetrees.res.txt")
CC_perm <- read.delim(perm_file, sep=" ", header=F)
names(CC_perm) <- c("geneid", "pval", "R2")

#loading files with GO annotation for geneids
CC_file <- paste0(dir, "/", cluster,"/", cluster,".CC.geneid.tsv")
CC <- read.delim(CC_file, header=T, stringsAsFactors = F,na.strings = "")

# adding GO annotation to permanova-df
CC_perm$CC <- NA
for(i in 1:length(CC_perm$geneid))
{
  gene <- CC_perm$geneid[i]
  if(gene %in% CC$geneid)
  {
    CC_perm$CC[i] <- CC$annotation[which(CC$geneid==gene)][1]
  }
}
CC_perm <- separate_rows(CC_perm,"CC",sep = "[;]")

#saving permanova R2 results in file together with CC annotation
CC_perm_file <- paste0(dir, "/R2_pvals/", cluster, ".CC_perm.tsv")
write.table(CC_perm,file = CC_perm_file, quote = F, row.names = F, sep = '\t')




#### Testing on Biological Process ####
#loading file with permanova results from all genetrees
perm_file <- paste0(dir, "/", cluster,"/city_gene_trees/permanova.citygenetrees.res.txt")
BP_perm <- read.delim(perm_file, sep=" ", header=F)
names(BP_perm) <- c("geneid", "pval", "R2")

#loading files with GO annotation for geneids
BP_file <- paste0(dir, "/", cluster,"/", cluster,".BP.geneid.tsv")
BP <- read.delim(BP_file, header=T, stringsAsFactors = F,na.strings = "")
 
# adding GO annotation to permanova-df
BP_perm$BP <- NA
for(i in 1:length(BP_perm$geneid))
{
  gene <- BP_perm$geneid[i]
  if(gene %in% BP$geneid)
  {
    BP_perm$BP[i] <- BP$annotation[which(BP$geneid==gene)][1]
  }
}
BP_perm <- separate_rows(BP_perm,"BP",sep = "[;]")

#saving permanova R2 results in file together with CC annotation
BP_perm_file <- paste0(dir, "/R2_pvals/", cluster, ".BP_perm.tsv")
write.table(BP_perm,file = BP_perm_file, quote = F, row.names = F, sep = '\t')
