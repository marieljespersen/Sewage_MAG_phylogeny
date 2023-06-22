library(stringr)
library(ggplot2)
library(reshape2)

setwd("/Users/maloj/HCP_Anywhere/GS2/Results/strain_mix")

fastANI <- read.delim("fastani.salmonella.tsv", sep="\t", header=F)
names(fastANI) <- c("strain1", "strain2", "ANI", "X", "Y")

fastANI$strain1 <- str_remove_all(fastANI$strain1, "data/references/")
fastANI$strain1 <- str_remove_all(fastANI$strain1, ".fa")

fastANI$strain2 <- str_remove_all(fastANI$strain2, "data/references/")
fastANI$strain2 <- str_remove_all(fastANI$strain2, ".fa")

heat_data <- subset(fastANI, select = -c(X,Y) )

ggplot(heat_data, aes(strain1, strain2)) +
  geom_tile(aes(fill = ANI), colour = "white") +
  scale_fill_gradient(low = "white", high = "red")

strain_df <- dcast(heat_data, strain1 ~ strain2)
rownames(strain_df) <- strain_df$strain1

strain_mat <- as.matrix(subset(strain_df, select = -strain1 ))
#heatmap(strain_mat)

ord <- hclust( dist(strain_mat, method = "euclidean"), method = "ward.D" )$order
heat_data <- melt(strain_mat)
names(heat_data) <- c("strain1", "strain2", "ANI")
str(heat_data)
heat_data$strain1 <- factor(heat_data$strain1, levels=rownames(strain_mat)[ord])
heat_data$strain2 <- factor(heat_data$strain2, levels=rownames(strain_mat)[ord])


p <- ggplot(heat_data, aes(strain1, strain2)) + 
  geom_tile(aes(fill=ANI)) +
  scale_fill_gradient(low="blue", high="red") + 
  theme(text = element_text(size = 25)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(plot=p,width=8,height=6,filename="/Users/maloj/HCP_Anywhere/GS2/Results/strain_mix/ANI_heatmap.pdf", device = "pdf", useDingbats=FALSE)
ggsave(plot=p,width=8,height=6,filename="/Users/maloj/HCP_Anywhere/GS2/Results/strain_mix/ANI_heatmap.jpg", device = "jpg")


ANI_inv <- as.matrix(subset(strain_df, select = -c(strain1, NC_010067, NC_015761, NC_021812, NC_021870)))                     
ANI_df <- as.data.frame(t(ANI_inv))
ANI_inv <- as.matrix(subset(ANI_df, select = -c(NC_010067, NC_015761, NC_021812, NC_021870))) 



