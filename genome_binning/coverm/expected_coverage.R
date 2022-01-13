#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

library(stringr)

# get inputs
R1 <- args[1] 
cov_file <- args[2]
out_file <- args[3] 

# read coverm table
cov_table <- read.table(cov_file, sep="\t", header=T, stringsAsFactors=F)
names(cov_table) <- c("genome", "Relative_abundance", "RPKM", "TPM", "read_count", "covered_fraction", "length")


# calculating average sequence length from forward reads
average_seq_length <- as.numeric(str_replace(system(paste("awk '{if(NR%4==2) {count++; bases += length} } END{print bases/count}'", R1), intern=T), ",", "."))

# getting coverage information for each MAG species in this sample
expected_coverage <- c(0)
for(i in 2:nrow(cov_table)) # skip the first row, that is unmapped
{
  genome_length <- cov_table$length[i]
  mapped_reads <- cov_table$read_count[i]
  c <- 1-(1-average_seq_length/genome_length)^mapped_reads
  expected_coverage <- c(expected_coverage, c)
}
cov_table <- cbind(cov_table, expected_coverage)
coverage_ratio <- cov_table$covered_fraction/cov_table$expected_coverage
coverage_difference <- cov_table$covered_fraction-cov_table$expected_coverage
cov_table <- cbind(cov_table, coverage_ratio)

write.table(cov_table, out_file, sep="\t", quote=F, row.names = F)



