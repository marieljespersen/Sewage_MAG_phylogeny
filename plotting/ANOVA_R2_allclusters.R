library(dplyr)
library(ggpubr)
library(tidyr)
library(openxlsx)
library(gridExtra)
library(ggplot2)
library(ggExtra)

setwd("/Users/maloj/HCP_Anywhere/GS2/Results/go_annotations")

# creating vector of cluster numbers for clusters that do have a gene permanova file
# the clusters not included did not havee multiple cities from multiple regions, and could not be tested
cluster_no <- c(1,2,3,4,5,6,7,8,9,11,12,13,14,15,16,17,18,20,21,22,23,24,25,28,29,30,31,32,34,37,38,41)
clusters <- c()
CC_pvals <- c()
BP_pvals <- c()
memR2 <- c()
naR2 <- c()
orgR2 <- c()
wilcox_pvals <- vector(mode = "list", length = length(cluster_no))
kruskal_pvals <- c()

#### Testing on Cellular Compartment ####
# save plots to file
#pdf(file = "R2.CC.memNAorg.pdf")
# file for saving summary values
#of="CC.memNAorg.sum.xlsx"
#OUT <- createWorkbook()

#cluster <- "C5"
for(i in 1:length(cluster_no))
{
  cluster <- paste0("C",cluster_no[i])
  clusters <- c(clusters, cluster)
  
  #loading file with permanova results from all genetrees
  perm_file <- paste0("R2_pvals/CC/", cluster,".CC_perm.tsv")
  CC_perm <- read.delim(file=perm_file, header=T, stringsAsFactors = F, na.strings = "")
  CC_perm$Cellular_Component <- "Metabolic"
  CC_perm$Cellular_Component[which(CC_perm$CC == "membrane")] <- "membrane"
  CC_perm$Cellular_Component[which(CC_perm$CC == "organelle")] <- "organelle"
  #CC_perm <- CC_perm[which(CC_perm$CC == "membrane" | CC_perm$CC =="organelle"),]

  CC_perm$Cellular_Component <- as.factor(CC_perm$Cellular_Component)
  CC_perm$R2 <- as.numeric(CC_perm$R2)
  CC_perm <- CC_perm[which(CC_perm$R2>=0),]
  levels(CC_perm$Cellular_Component)

  # get summary statistics by phylum
  CC_sum <- group_by(CC_perm,Cellular_Component) %>%
    dplyr::summarise(
      count = n(),
      mean = mean(R2, na.rm = TRUE),
      sd = sd(R2, na.rm = TRUE),
      median = median(R2, na.rm = TRUE),
      IQR = IQR(R2, na.rm = TRUE)
    )
  memR2 <- c(memR2,CC_sum$mean[1])
  naR2 <- c(naR2,CC_sum$mean[2])
  orgR2 <- c(orgR2,CC_sum$mean[3])
  #addWorksheet(OUT, cluster)
  #writeData(OUT, sheet = cluster, x = CC_sum)

  #### Box plots ####
  # Plot R2 by gene group
  #pdf("/Users/maloj/HCP_Anywhere/GS2/Results/C5_R2_boxplot.pdf", width=4, height=4)
  xlabs <- paste(levels(CC_perm$Cellular_Component),"\n(N=",table(CC_perm$Cellular_Component),")",sep="")
  CC_plot <- ggboxplot(CC_perm, x = "Cellular_Component", y = "R2", color = "Cellular_Component", na.rm=T) +
    scale_x_discrete(labels=xlabs) +
    labs(title = cluster) + ylab("R2") + xlab("Cellular component") +
    theme_bw() +
    theme(legend.position = "none") +
    removeGrid() +
    theme(text = element_text(size = 15))
  print(CC_plot)
  #dev.off()

  #### Testing #####
  # test for significant differences between the different phylum
  CC_test <- kruskal.test(R2 ~ Cellular_Component, data = CC_perm)
  kruskal_pvals <- c(kruskal_pvals, format.pval(CC_test$p.value))

  # if  Kruskal-Wallis test shows siginificance, use pairwise wilcox test
  wilcox_pvals[[i]] <- pairwise.wilcox.test(CC_perm$R2, CC_perm$Cellular_Component)
 
}
#saveWorkbook(OUT,of)
#dev.off()

kruskal_pvals <- as.numeric(kruskal_pvals)
kruskal_pvals <- p.adjust(kruskal_pvals, method="fdr")
kruskal_pvals <- round(kruskal_pvals, digits = 3)
length(kruskal_pvals[which(kruskal_pvals<0.05)])
length(kruskal_pvals)
kruskal_df <- as.data.frame(kruskal_pvals)
rownames(kruskal_df) <- clusters
write.table(kruskal_df,"kruskal_pvals.tsv", sep="\t", row.names = TRUE, quote=FALSE)
clusters[which(kruskal_pvals<0.05)]

names(wilcox_pvals) <- clusters
# check wilcox result for clusters that are significant in kruskal wallis test
wp <- as.data.frame(matrix(nrow=9, ncol=3))
names(wp) <- c("mem_NA","mem_org", "org_NA")
rownames(wp) <- clusters[which(kruskal_pvals<0.05)]
wp["C3","mem_NA"] <- wilcox_pvals$C3$p.value[1]
wp["C3","mem_org"] <- wilcox_pvals$C3$p.value[2]
wp["C3","org_NA"] <- wilcox_pvals$C3$p.value[4]
wp["C5","mem_NA"] <- wilcox_pvals$C5$p.value[1]
wp["C5","mem_org"] <- wilcox_pvals$C5$p.value[2]
wp["C5","org_NA"] <- wilcox_pvals$C5$p.value[4]
wp["C8","mem_NA"] <- wilcox_pvals$C8$p.value[1]
wp["C8","mem_org"] <- wilcox_pvals$C8$p.value[2]
wp["C8","org_NA"] <- wilcox_pvals$C8$p.value[4]
wp["C13","mem_NA"] <- wilcox_pvals$C13$p.value[1]
wp["C13","mem_org"] <- wilcox_pvals$C13$p.value[2]
wp["C13","org_NA"] <- wilcox_pvals$C13$p.value[4]
wp["C14","mem_NA"] <- wilcox_pvals$C14$p.value[1]
wp["C14","mem_org"] <- wilcox_pvals$C14$p.value[2]
wp["C14","org_NA"] <- wilcox_pvals$C14$p.value[4]
wp["C24","mem_NA"] <- wilcox_pvals$C24$p.value[1]
wp["C24","mem_org"] <- wilcox_pvals$C24$p.value[2]
wp["C24","org_NA"] <- wilcox_pvals$C24$p.value[4]
wp["C34","mem_NA"] <- wilcox_pvals$C34$p.value[1]
wp["C34","mem_org"] <- wilcox_pvals$C34$p.value[2]
wp["C34","org_NA"] <- wilcox_pvals$C34$p.value[4]
wp["C37","mem_NA"] <- wilcox_pvals$C37$p.value[1]
wp["C37","mem_org"] <- wilcox_pvals$C37$p.value[2]
wp["C37","org_NA"] <- wilcox_pvals$C37$p.value[4]
wp["C38","mem_NA"] <- wilcox_pvals$C38$p.value[1]
wp["C38","mem_org"] <- wilcox_pvals$C38$p.value[2]
wp["C38","org_NA"] <- wilcox_pvals$C38$p.value[4]
wp$cluster <- rownames(wp)

molten_wp <- melt(wp)
molten_wp$value <- p.adjust(molten_wp$value, method="BH")
molten_wp$value <- round(molten_wp$value, digits = 3)
molten_wp[which(molten_wp$variable=="mem_NA"),]

# extract p values from wilcox test
mem_org_pvals <- c()
mem_NA_pvals <- c()
NA_org_pvals <- c()
for(i in which(kruskal_pvals<0.05))
{
  print(i)
  mem_org <- wilcox_pvals$C3$p.value[1,1]
  mem_NA <- wilcox_pvals$C3$p.value[2,1]
  NA_org <- wilcox_pvals$C3$p.value[2,2]
  mem_org_pvals <- c(mem_org_pvals, mem_org)
  mem_NA_pvals <- c(mem_NA_pvals, mem_NA)
  NA_org_pvals <- c(NA_org_pvals, mem_org)
}

mem_org <- wilcox_pvals$C3$p.value[1,1]
mem_NA <- wilcox_pvals$C3$p.value[2,1]
NA_org <- wilcox_pvals$C3$p.value[2,2]


#### Testing on Biological Process ####
pdf(file = "R2.BP.all.pdf")
# file for saving summary values
of="BP_all_sum.xlsx"
OUT <- createWorkbook()

for(i in 1:length(cluster_no))
{
  cluster <- paste0("C",cluster_no[i])
  
  #loading file with permanova results from all genetrees
  perm_file <- paste0("R2_pvals/", cluster,".BP_perm.tsv")
  BP_perm <- read.delim(file=perm_file, header=T, stringsAsFactors = F)
  
  BP_perm$BP <- as.factor(BP_perm$BP)
  BP_perm$R2 <- as.numeric(BP_perm$R2)
  levels(BP_perm$BP)
  
  # get summary statistics by phylum
  BP_sum <- group_by(BP_perm, BP) %>%
    summarise(
      count = n(),
      mean = mean(R2, na.rm = TRUE),
      sd = sd(R2, na.rm = TRUE),
      median = median(R2, na.rm = TRUE),
      IQR = IQR(R2, na.rm = TRUE)
    )
  addWorksheet(OUT, cluster)
  writeData(OUT, sheet = cluster, x = BP_sum)
  
  #### Box plots ####
  # Plot R2 by biological process
  BP_plot <- ggboxplot(BP_perm, x = "BP", y = "R2", 
                       color = "BP", na.rm=T) + 
    theme(axis.text.x=element_text(color = "black", size=11, angle=30, vjust=.8, hjust=0.8)) +
    labs(title = cluster) + ylab("R2") + xlab("Biological Process")
  print(BP_plot)
  
  
  #### Testing #####
  # test for significant differences between the different phylum
  BP_test <- kruskal.test(R2 ~ BP, data = BP_perm)
  BP_pvals <- c(BP_pvals, format.pval(BP_test$p.value))
  
  # if  Kruskal-Wallis test shows siginificance, use pairwise wilcox test
  #pairwise.wilcox.test(CC_perm$R2, CC_perm$CC, p.adjust.method = "BH")
  
}
saveWorkbook(OUT,of)
dev.off()

# adjust p values
BP_pvals <- format.pval(p.adjust(BP_pvals, method="fdr"), digits=3)
CC_pvals <- format.pval(p.adjust(CC_pvals, method="fdr"), digits=3)

# write table with p-values to output file
write.table(cbind(clusters,CC_pvals,BP_pvals), file = "anova_R2_pvals.tsv", quote = F, row.names = F, sep = '\t')
write.table(cbind(clusters,CC_pvals), file = "CC_sub.anova_R2_pvals.tsv", quote = F, row.names = F, sep = '\t')


##### plot histograms of all R2 values
# save plots to file
pdf(file = "R2.hist.pdf")
for(i in 1:length(cluster_no))
{
  cluster <- paste0("C",cluster_no[i])
  clusters <- c(clusters, cluster)
  
  #loading file with permanova results from all genetrees
  perm_file <- paste0("R2_pvals/", cluster,".CC_perm.tsv")
  perm <- read.delim(file=perm_file, header=T, stringsAsFactors = F, na.strings = "")
  perm$R2 <- as.numeric(perm$R2)
  
  # Plot R2 histogram
  qplot(perm$R2, geom="histogram") 
  plot <- qplot(perm$R2, geom="histogram") +
    labs(title = cluster) + ylab("Frequency") + xlab("R2")
  print(plot)
}
dev.off()

