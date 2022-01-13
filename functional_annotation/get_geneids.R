#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

library(stringr)

cluster <- args[1]
dir <- args[2]
sonic_file <- args[3]
genetrees_file <- paste0(dir, "/genetrees.txt")
CC_file <- paste0(dir, "/", cluster, ".CC.tsv")
BP_file <- paste0(dir, "/", cluster, ".BP.tsv")

# read orthologous file
sonic <- read.delim(sonic_file, header=T, stringsAsFactors = F,na.strings = "")
genes <- sonic[seq(5,length(sonic),2)] #genes are every second  column
# insert geneids as first column in genes dataframe
rownames(genes) <- paste0("gene_",sonic$group_id)
# read genetrees and subset genes dataframe to contain only genes that are made as trees
genetrees <- as.vector(as.matrix(read.delim(genetrees_file, header=F)))
genes <- genes[rownames(genes) %in% genetrees,]

# read GO annotation files
CC <- read.delim(CC_file, header=F, stringsAsFactors = F, na.strings = "")
names(CC) <- c("Protein","GO","ontology","annotation")
BP <- read.delim(BP_file, header=F, stringsAsFactors = F, na.strings = "")
names(BP) <- c("Protein","GO","ontology","annotation")

#  function to correct names to samples in stead of bins
correct.names <- function(name_list) {
  new_names <- c()
  for(n in name_list)
  {
    if(str_detect(n, "_vambbin_"))
    {
      sample <- strsplit(n, "_vambbin_")[[1]][1]
    }
    else if(str_detect(n, "_metabin_"))
    {
      sample <- strsplit(n, "_metabin_")[[1]][1]
    }
    new_names <- c(new_names, sample)
  }
  return(new_names)
}
names(genes) <- correct.names(names(genes))

# add geneid to cellular compartment df
CC$geneid <- NA
for(i in 1:length(CC$Protein))
{
  p <- CC$Protein[i]
  sample <- str_replace(strsplit(p, "_NODE_")[[1]][1], "-", "_")
  prot_list <- genes[[sample]]
  if(p %in% prot_list)
  {
    CC$geneid[i] <- rownames(genes[which(prot_list==p),])
  }
}

# add geneid to Biological Process df
BP$geneid <- NA
for(i in 1:length(BP$Protein))
{
  p <- BP$Protein[i]
  sample <- str_replace(strsplit(p, "_NODE_")[[1]][1], "-", "_")
  prot_list <- genes[[sample]]
  if(p %in% prot_list)
  {
    BP$geneid[i] <- rownames(genes[which(prot_list==p),])
  }
}

CC_outfile  <- paste0(dir, "/", cluster, ".CC.geneid.tsv")
BP_outfile <- paste0(dir, "/", cluster, ".BP.geneid.tsv")

write.table(CC,file = CC_outfile, quote = F, row.names = F, sep = '\t')
write.table(BP,file = BP_outfile, quote = F, row.names = F, sep = '\t')
