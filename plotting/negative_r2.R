# library
library(ggplot2)
library(plyr)
library(gridExtra)
library(stringr)
library(dplyr)
library(ggpubr)
library(tidyr)
library(ggExtra)


########## køres på C2 ############

clusters <- c("C1","C2","C3","C4","C5","C6","C7","C8","C9","C11","C12","C13","C14","C15","C16","C17","C18","C20","C21","C22","C23","C24","C25","C28","C29","C30","C31","C32","C34","C37","C38","C41")
all_perm_dists <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(all_perm_dists) <- c("geneid", "pval", "R2", "CC", "variance", "Cellular_Component","cluster")

for(i in 1:length(clusters))
{
  c <- clusters[i]
  
  # read in permanova results
  perm_file <- paste0("/home/projects/cge/data/projects/5001/Binning_vamb/clusters/R2_pvals/", c, ".CC_perm.tsv")
  CC_perm <- read.delim(file=perm_file, header=T, stringsAsFactors = F, na.strings = "")
  CC_perm$R2 <- as.numeric(CC_perm$R2)
  
  # 
  div_file <- paste0("/home/projects/cge/data/projects/5001/Binning_vamb/clusters/R2_bias/gene_div/", c, ".gene_div.txt")
  length_file <-  paste0("/home/projects/cge/data/projects/5001/Binning_vamb/clusters/R2_bias/gene_lengths/", c, ".gene_lengths.txt")
  entropy_file <-  paste0("/home/projects/cge/data/projects/5001/Binning_vamb/clusters/R2_bias/reg_entropy/", c, ".reg_entropy.txt")
  sample_file <-  paste0("/home/projects/cge/data/projects/5001/Binning_vamb/clusters/citygenetreee_permanovas/", c, ".perm.cgt.txt")
  
  # read in all tables for one cluster
  gene_div <- read.table(div_file, header=T, sep = " ", fill=TRUE, stringsAsFactors = F)
  gene_div$geneid <- str_replace(gene_div$geneid, ".aln", "")
  gene_lengths <- read.table(length_file, header=T, sep = " ", fill=TRUE, stringsAsFactors = F)
  reg_entropy <- read.table(entropy_file, header=F, sep = " ", fill=TRUE, stringsAsFactors = F)
  names(reg_entropy) <- c("geneid", "entropy")
  no_samples <- read.table(sample_file, header=F, sep = " ", fill=TRUE, stringsAsFactors = F)
  names(no_samples) <- c("geneid", "p_val","R2","samples")
  
  #combine into one df
  combined_gene_prop <- merge(x = CC_perm, y = gene_div,'geneid')
  combined_gene_prop <- merge(x = combined_gene_prop, y = no_samples,'geneid')
  combined_gene_prop <- merge(x = combined_gene_prop, y = gene_lengths,'geneid')
  combined_gene_prop <- merge(x = combined_gene_prop, y = reg_entropy,'geneid')

  
  combined_gene_prop$Cellular_Component <- "Metabolic"
  combined_gene_prop$Cellular_Component[which(combined_gene_prop$CC == "membrane")] <- "Membrane"
  combined_gene_prop$Cellular_Component[which(combined_gene_prop$CC == "organelle")] <- "Organelle"
  combined_gene_prop$Cellular_Component <- as.factor(combined_gene_prop$Cellular_Component)
  combined_gene_prop$cluster <- c
  
  all_perm_dists <- rbind(all_perm_dists, combined_gene_prop)
  
}

write.table(all_perm_dists, file="R2_bias_all_clusters.tsv", quote=FALSE, sep="\t", row.names = F)


########## køres lokalt ############

all_R2_bias <- read.table("/Users/maloj/HCP_Anywhere/GS2/Results/R2_bias/R2_bias_all_clusters.tsv", sep="\t", header=T)

# plot for all clusters
# R2 vs gene variance
p1  <- ggplot(all_R2_bias, aes(x=R2.x, y=variance)) + 
  geom_point() +
  xlab(expression(R^2 ~ "geography")) + ylab("Gene variance") +
  theme_bw() +  removeGrid()
p1  <- p1 + geom_vline(xintercept = 0, color = "red")

# R2 vs number of samples
p2  <- ggplot(all_R2_bias, aes(x=R2.x, y=samples)) + 
  geom_point() +
  xlab(expression(R^2 ~ "geography")) + ylab("Samples")+
  theme_bw() +  removeGrid()
p2  <- p2 + geom_vline(xintercept = 0, color = "red")

# R2 vs gene length
p3  <- ggplot(all_R2_bias, aes(x=R2.x, y=length)) + 
  geom_point() +
  xlab(expression(R^2 ~ "geography")) + ylab("Gene length")+
  theme_bw() +  removeGrid()
p3  <- p3 + geom_vline(xintercept = 0, color = "red")

# R2 vs regional entropy
p4  <- ggplot(all_R2_bias, aes(x=R2.x, y=entropy)) + 
  geom_point() +
  xlab(expression(R^2 ~ "geography")) + ylab("Regional entropy")+
  theme_bw() +  removeGrid()
p4  <- p4 + geom_vline(xintercept = 0, color = "red")

p <- grid.arrange(p1, p2, p3, p4)
ggsave(plot=p,width=10,height=10,filename="/Users/maloj/HCP_Anywhere/GS2/Results/R2_bias/negative_R2.png",device = "png")


# testing seperately for each cluster
#c  <- levels(all_R2_bias$cluster)[27]
all_R2_bias <- all_R2_bias[which(all_R2_bias$R2.x>=0),]
all_cor <- matrix(, nrow=length(levels(all_R2_bias$cluster)), ncol=4)
rownames(all_cor) <- levels(all_R2_bias$cluster)
colnames(all_cor) <- c("gene_variance", "gene_length", "number_of_samples", "regional_entropy")
for(c in levels(all_R2_bias$cluster))
{
  cluster_subset <- all_R2_bias[which(all_R2_bias$cluster==c),]
  gene_var_cor <- cor(cluster_subset$R2.x[which(!is.na(cluster_subset$R2.x))],cluster_subset$variance[which(!is.na(cluster_subset$R2.x))])
  sample_cor <- cor(cluster_subset$R2.x[which(!is.na(cluster_subset$R2.x))],cluster_subset$samples[which(!is.na(cluster_subset$R2.x))])
  entropy_cor <- cor(cluster_subset$R2.x[which(!is.na(cluster_subset$R2.x))],cluster_subset$entropy[which(!is.na(cluster_subset$R2.x))])
  gene_length_cor <- cor(cluster_subset$R2.x[which(!is.na(cluster_subset$R2.x))],cluster_subset$length[which(!is.na(cluster_subset$R2.x))])
  
  all_cor[c,] <- c(gene_var_cor, gene_length_cor, sample_cor, entropy_cor)
}

