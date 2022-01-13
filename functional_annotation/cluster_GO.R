#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

library(GO.db)
library(tidyr)
library(dplyr)

cluster <- args[1]
dir <- args[2]
in_file <- paste0(dir, "/", cluster,".GO_annotations.tsv")

#### Read in GO annotations for cluster (eg. C1) ####
#setwd(work_dir)
GO_cluster <- read.delim(in_file, header=F, stringsAsFactors = F,na.strings = "")
names(GO_cluster) <- c("Protein accession","Sequence MD5 digest","Sequence length", "Analysis", "Signature accession", "Signature description", "Start location", "Stop location", "Score", "Status", "Date", "InterPro annotations - accession", "InterPro annotations - description","GO annotations", "Pathways annotations")
GO_sep <- separate_rows(GO_cluster,14,sep = "[|]")
GO_sep <- GO_sep[which(!is.na(GO_sep$`GO annotations`)),]
GO_sep$ontology <- Ontology(GOTERM[GO_sep$`GO annotations`])
GO_sub <- GO_sep[names(GO_sep) %in% c("Protein accession", "GO annotations", "ontology")]

# functions to create list to identify "top layer" term for each GO term
# Seperate functions for each of the three different GO term aspects
get.CC <- function(goid) 
{
  ancestors <- c()
  if(length(GOCCANCESTOR[[goid]])==2)
  {
    ancestors <- c(ancestors,Term(GOTERM[[goid]]))
  }
  else
  {
    # tjek for all ancestors which one only has two ancestors
    for(a in GOCCANCESTOR[[goid]])
    {
      if(length(GOCCANCESTOR[[a]])==2)
      {
        ancestors <- c(ancestors,Term(GOTERM[[a]]))
      }
    }
  }
  return(ancestors)
}

get.BP <- function(goid) 
{
  ancestors <- c()
  if(length(GOBPANCESTOR[[goid]])==2)
  {
    ancestors <- c(ancestors,Term(GOTERM[[goid]]))
  }
  else
  {
    # tjek for all ancestors which one only has two ancestors
    for(a in GOBPANCESTOR[[goid]])
    {
      if(length(GOBPANCESTOR[[a]])==2)
      {
        ancestors <- c(ancestors,Term(GOTERM[[a]]))
      }
    }
  }
  return(ancestors)
}

get.MF <- function(goid) 
{
  ancestors <- c()
  if(length(GOMFANCESTOR[[goid]])==2)
  {
    ancestors <- c(ancestors,Term(GOTERM[[goid]]))
  }
  else
  {
    # tjek for all ancestors which one only has two ancestors
    for(a in GOMFANCESTOR[[goid]])
    {
      if(length(GOMFANCESTOR[[a]])==2)
      {
        ancestors <- c(ancestors,Term(GOTERM[[a]]))
      }
    }
  }
  return(ancestors)
}


### get higher annotation for all go-terms found in CC, BP, or MF
# split df into three, for each of the GO term aspects
GO_ont <- split(GO_sub, GO_sub$ontology)

# identify the GO terms present in the data
CC_goids <- unique(GO_ont$CC$`GO annotations`)
BP_goids <- unique(GO_ont$BP$`GO annotations`)
MF_goids <- unique(GO_ont$MF$`GO annotations`)

# subset df to contain one entry per gene/protein
GO_ont$CC <- GO_ont$CC[ !duplicated(GO_ont$CC$`Protein accession`), ]
GO_ont$BP <- GO_ont$BP[ !duplicated(GO_ont$BP$`Protein accession`), ]
GO_ont$MF <- GO_ont$MF[ !duplicated(GO_ont$MF$`Protein accession`), ]

# create list to look-up higher annotation for GO terms in the dataset
CC_ann <- lapply(CC_goids, get.CC)
names(CC_ann) <- CC_goids
BP_ann <- lapply(BP_goids, get.BP)
names(BP_ann) <- BP_goids
MF_ann <- lapply(MF_goids, get.MF)
names(MF_ann) <- MF_goids

# look-up higher annotation level for each gene/protein
CC.ann <- function(goid) 
{
  return(CC_ann[[goid]])
}
GO_ont$CC$`Cellular Compartment` <- lapply(GO_ont$CC$`GO annotations`, CC.ann)
GO_ont$CC$`Cellular Compartment` <- sapply(GO_ont$CC$`Cellular Compartment`, paste, collapse=";")

BP.ann <- function(goid) 
{
  return(BP_ann[[goid]])
}
GO_ont$BP$`Biological Process` <- lapply(GO_ont$BP$`GO annotations`, BP.ann)
GO_ont$BP$`Biological Process` <- sapply(GO_ont$BP$`Biological Process`, paste, collapse=";")

MF.ann <- function(goid) 
{
  return(MF_ann[[goid]])
}
GO_ont$MF$`Molecular Function` <- lapply(GO_ont$MF$`GO annotations`, MF.ann)
GO_ont$MF$`Molecular Function` <- sapply(GO_ont$MF$`Molecular Function`, paste, collapse=";")

# save the 3 GO tables in different files
out_CC <- paste0(dir, "/", cluster, ".CC.tsv")
out_BP <- paste0(dir, "/", cluster, ".BP.tsv")
out_MF <- paste0(dir, "/", cluster, ".MF.tsv")

write.table(GO_ont$CC,file = out_CC, quote = F, row.names = F, sep = '\t',col.names = F)
write.table(GO_ont$BP,file = out_BP, quote = F, row.names = F, sep = '\t',col.names = F)
write.table(GO_ont$MF,file = out_MF, quote = F, row.names = F, sep = '\t',col.names = F)
