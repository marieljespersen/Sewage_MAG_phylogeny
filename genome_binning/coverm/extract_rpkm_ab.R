#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

# get inputs
cov_file <- args[1] 
TPM_out_file <- args[2]
read_count_out_file <- args[3] 

cov_table <- read.table(cov_file, sep="\t", header=T, stringsAsFactors=F)
splittet_filename <- strsplit(cov_file,"/" )
samplename <- splittet_filename[[1]][length(splittet_filename[[1]])]

TPM <- c(samplename)
read_count <- c(samplename)
for(b in 2:nrow(cov_table))
{
  if(is.na(cov_table$coverage_ratio[b])) #if no mapped reads
  {
    TPM <- c(TPM, 0)
    read_count <- c(read_count, 0)
  }
  else if(cov_table$coverage_ratio[b]>=0.5)
  {
    TPM <- c(TPM, cov_table$TPM[b])
    read_count <- c(read_count, cov_table$read_count[b])
  }
  else
  {
    TPM <- c(TPM, 0)
    read_count <- c(read_count, 0)
  }
}

write.table(TPM, TPM_out_file, sep="\t", quote=F, row.names=F, col.names=F)
write.table(read_count, read_count_out_file, sep="\t", quote=F, row.names=F, col.names=F)
