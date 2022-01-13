library(ggplot2)
library(reshape2)
library(ggExtra)

setwd("/Users/maloj/HCP_Anywhere/GS2/Results/oxy_genes")
GS2_oxy <- read.delim("GS.oxy_genes.txt", stringsAsFactors = F, sep=" ", header=F)
alm_oxy <- read.delim("alm.oxy_genes.txt", stringsAsFactors = F, sep=" ", header=F)


df <- as.data.frame(table(GS2_oxy$V2))
names(df) <- c("genes", "Sewage")
df$human <- c(40,0,0,0,0)

df2 <- melt(df, id.vars="genes")
names(df2) <- c("genes", "source", "frequency")

pdf(file="/Users/maloj/HCP_Anywhere/GS2/Results/oxy_genes/oxy.pdf", width=4, height=3)
ggplot(df2, aes(x=genes, y=frequency, fill=source)) +
  geom_bar(stat='identity', position='dodge') + 
  labs(x="Number of genes observed", y="Frequency") +
  theme_bw() +
  removeGrid() +
  theme(text = element_text(size = 10)) 
dev.off()
