############### reading libaries #############
library(stringr)
library(dplyr)
library(ggpubr)
library(tidyr)
library(gridExtra)
library(ggplot2)
library(plyr)
library(ggExtra)
library(reshape2)
library(forcats)


############### reading data #############
clusters <- c("C3", "C5", "C8", "C13", "C14", "C24", "C34", "C37", "C38")

summaries <- list()
dnds_all <- list()
for(i in 1:length(clusters))
{
  dnds_file <- paste0("~/HCP_Anywhere/GS2/Results/dnds/codeml/", clusters[i], ".dnds.txt")
  dnds_cluster <- read.delim(dnds_file, stringsAsFactors = F, na.strings = "", sep=",", header=F)
  names(dnds_cluster) <- c("geneid", "dnds")
  
  GO_file <- paste0("/Users/maloj/HCP_Anywhere/GS2/Results/go_annotations/R2_pvals/CC/", clusters[i], ".CC_perm.tsv")
  GO_cluster <- read.delim(file=GO_file, header=T, stringsAsFactors = F, na.strings = "")
  # remove genes with negative R2?
  #GO_cluster <- GO_cluster[which(GO_cluster$R2>=0),]
  dnds_all[[i]] <- merge(GO_cluster, dnds_cluster, by="geneid")
  
  # remove extreme values and NAs
  #dnds_all[[i]]$dnds[dnds_all[[i]]$dnds==0.00010] <- NA
  #dnds_all[[i]]$dnds[dnds_all[[i]]$dnds==999.00000] <- NA
  dnds_all[[i]] <- dnds_all[[i]][!is.na(dnds_all[[i]]$dnds),]
  
  # add variable to distinguish membrane, organelle, and metanolic
  dnds_all[[i]]$Cellular_Component <- "Metabolic"
  dnds_all[[i]]$Cellular_Component[which(dnds_all[[i]]$CC == "membrane")] <- "Membrane"
  dnds_all[[i]]$Cellular_Component[which(dnds_all[[i]]$CC == "organelle")] <- "Organelle"
  dnds_all[[i]]$Cellular_Component <- as.factor(dnds_all[[i]]$Cellular_Component)
  
  # add variable with cluster name for later merging
  dnds_all[[i]]$cluster <- clusters[i]
  
  summaries[[i]] <- dnds_all[[i]] %>%
    group_by(Cellular_Component, add=TRUE) %>%
    dplyr::summarise(
      count = n(),
      mean = mean(dnds, na.rm = TRUE),
      sd = sd(dnds, na.rm = TRUE),
      median = median(dnds, na.rm = TRUE),
      IQR = IQR(dnds, na.rm = TRUE)
    )
  
}

all_merged <- rbind(dnds_all[[1]], dnds_all[[2]], dnds_all[[3]], dnds_all[[4]], dnds_all[[5]], dnds_all[[6]], dnds_all[[7]], dnds_all[[8]], dnds_all[[9]])

for(i in 1:length(clusters))
{
  print(clusters[i])
  print(nrow(all_merged[which(all_merged$cluster==clusters[i]),]))
}

################## summary #######################
# get summary statistics by cluster
dnds_diff <- group_by(all_merged,cluster) %>%
group_by(Cellular_Component, add=TRUE) %>%
  dplyr::summarise(
    count = n(),
    mean = mean(dnds, na.rm = TRUE),
    sd = sd(dnds, na.rm = TRUE),
    median = median(dnds, na.rm = TRUE),
    IQR = IQR(dnds, na.rm = TRUE)
  )
cluster_diff <- split(dnds_diff, dnds_diff$cluster)


write.table(file="~/HCP_Anywhere/GS2/Results/dnds/codeml/dnds_nofilter_diff.csv", dnds_diff, sep=",", quote = F, row.names=F) 

#### look #####

# clusters with mem mean larger than organelle
for(i in 1:length(cluster_diff))
{
  #print(cluster_diff[[i]])
  if(cluster_diff[[i]]$mean[1]>cluster_diff[[i]]$mean[3])
  {
    print(cluster_diff[[i]]$cluster)
  }
}

# clusters with mem mean larger than NA
for(i in 1:length(cluster_diff))
{
  #print(cluster_diff[[i]])
  if(cluster_diff[[i]]$mean[1]>cluster_diff[[i]]$mean[2])
  {
    print(cluster_diff[[i]]$cluster)
  }
}

# clusters with mem median larger than organelle
for(i in 1:length(cluster_diff))
{
  #print(cluster_diff[[i]])
  if(cluster_diff[[i]]$median[1]>cluster_diff[[i]]$median[3])
  {
    print(cluster_diff[[i]]$cluster)
  }
}

# clusters with mem median larger NA
for(i in 1:length(cluster_diff))
{
  #print(cluster_diff[[i]])
  if(cluster_diff[[i]]$median[1]>cluster_diff[[i]]$median[2])
  {
    print(cluster_diff[[i]]$cluster)
  }
}


###### plotting #########
all_merged$cluster<- factor(all_merged$cluster, levels=c("C3",  "C5",  "C8", "C13", "C14", "C24", "C34", "C37", "C38")) # 

# boxplot log
p <- ggplot(data = all_merged, aes(x=Cellular_Component, y=dnds, na.rm=T, color=Cellular_Component)) + 
  scale_y_continuous(trans='log10') +
  xlab("") +
  ylab("log(dN/dS)") +
  geom_boxplot() + 
  facet_wrap( ~ cluster, scales="free", nrow = 2, labeller = labeller(cluster = c("C3"="C3*","C5"="C5*","C8"="C8","C13"="C13","C14"="C14***","C24"="C24","C34"="C34***","C37"="C37**","C38"="C38***")))  +
  theme_bw() +
  removeGrid() + 
  theme(legend.position="top") + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  geom_jitter(width = 0.1)
p
ggsave(plot=p,width=18,height=10,units="cm",filename="/Users/maloj/HCP_Anywhere/GS2/Results/dnds/codeml/dnds_boxplots.pdf", useDingbats=FALSE)


#dodge <- position_dodge(width = 1)
#geom_boxplot(width=.1, outlier.colour=NA, position = dodge) + geom_violin(position = dodge) 
# boxplot regular
ggplot(data = all_merged, aes(x=cluster, y=dnds, na.rm=T)) + 
  geom_boxplot(aes(color=Cellular_Component)) +
  xlab("") +
  ylab("log(dN/dS)") +
  theme_bw() +
  removeGrid() + facet_wrap( ~ cluster, scales="free", ncol = 3) 

# density plots for each clsuter
p <- ggplot(all_merged, aes(x=dnds, color=Cellular_Component)) +
  geom_density() +
  theme_bw() +
  removeGrid()
p <- p + facet_wrap( ~ cluster, scales="free", ncol = 3)
p

#### Testing #####
# splitting again to test
all_split <- split(all_merged, all_merged$cluster)

######### test for significant differences between the different phylum ################

kruskal_pvals <- c()
wilcox_pvals <- c()
for(c in clusters)
{
  CC_test <- kruskal.test(dnds ~ Cellular_Component, data = all_split[[c]])
  kruskal_pvals <- c(kruskal_pvals, CC_test$p.value)
  
  # if  Kruskal-Wallis test shows siginificance, use pairwise wilcox test
  wilcox_pvals[[c]] <- pairwise.wilcox.test(all_split[[c]]$dnds, all_split[[c]]$Cellular_Component)
}
kruskal_pvals <- as.numeric(kruskal_pvals)
kruskal_pvals <- p.adjust(kruskal_pvals, method="fdr")
kruskal_pvals <- round(kruskal_pvals, digits = 3)
#kruskal_pvals[which(kruskal_pvals<0.05)]
kruskal_df <- as.data.frame(kruskal_pvals)
rownames(kruskal_df) <- clusters

clusters[which(kruskal_df$kruskal_pvals<0.05)]
#write.table(file="~/HCP_Anywhere/GS2/Results/dnds/kruskal.csv", kruskal_df, sep=",", quote = F) 



# check wilcox result for clusters that are significant in kruskal wallis test
wp <- as.data.frame(matrix(nrow=length(clusters[which(kruskal_df$kruskal_pvals<0.05)]), ncol=3))
names(wp) <- c("mem_NA","mem_org", "org_NA")
rownames(wp) <- clusters[which(kruskal_pvals<0.05)]

wp["C3","mem_NA"] <- wilcox_pvals$C3$p.value[1]
wp["C3","mem_org"] <- wilcox_pvals$C3$p.value[2]
wp["C3","org_NA"] <- wilcox_pvals$C3$p.value[4]

wp["C5","mem_NA"] <- wilcox_pvals$C5$p.value[1]
wp["C5","mem_org"] <- wilcox_pvals$C5$p.value[2]
wp["C5","org_NA"] <- wilcox_pvals$C5$p.value[4]

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

wp["C38","mem_NA"] <- wilcox_pvals$C38$p.value[1]
wp["C38","mem_org"] <- wilcox_pvals$C38$p.value[2]
wp["C38","org_NA"] <- wilcox_pvals$C38$p.value[4]

wp["C37","mem_NA"] <- wilcox_pvals$C37$p.value[1]
wp["C37","mem_org"] <- wilcox_pvals$C37$p.value[2]
wp["C37","org_NA"] <- wilcox_pvals$C37$p.value[4]

wp$cluster <- rownames(wp)


molten_wp <- melt(wp)
molten_wp$value <- p.adjust(molten_wp$value, method="BH")
molten_wp$value <- round(molten_wp$value, digits = 3)
molten_wp[which(molten_wp$variable=="mem_org" & molten_wp$value<0.05 & molten_wp$value>0.01),]
molten_wp[which(molten_wp$variable=="mem_org" & molten_wp$value<0.01 & molten_wp$value>0.001),]
molten_wp[which(molten_wp$variable=="mem_org" & molten_wp$value<0.001),]


