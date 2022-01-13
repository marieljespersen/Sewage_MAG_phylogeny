#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

library(ape)

in_file <- args[1]
out_file <-args[2]

tr <- read.tree(in_file)

unrooted_tr <- unroot(tr)
write.tree(unrooted_tr, out_file)
