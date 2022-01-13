#read libraries
library(ggplot2)
library(plyr)
library(gridExtra)
library(stringr)
library(dplyr)
library(ggpubr)
library(tidyr)
library(reshape2)

clusters <- c("C1","C2","C3","C4","C5","C6","C7","C8","C9","C11","C12","C13","C14","C15","C16","C17","C18","C20","C21","C22","C23","C24","C25","C28","C29","C30","C31","C32","C34","C37","C38","C41")
all_perm_dists <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(all_perm_dists) <- c("geneid", "pval", "R2", "CC", "variance", "Cellular_Component","cluster")
wilcox_pvals <- vector(mode = "list", length = length(clusters))
kruskal_pvals <- c()

for(i in 1:length(clusters))
{
  c <- clusters[i]
  
  perm_file <- paste0("/Users/maloj/HCP_Anywhere/GS2/Results/go_annotations/R2_pvals/CC/", c, ".CC_perm.tsv")
  CC_perm <- read.delim(file=perm_file, header=T, stringsAsFactors = F, na.strings = "")
  CC_perm$R2 <- as.numeric(CC_perm$R2)
  
  ###### gene variation ######
  div_file <- paste0("/Users/maloj/HCP_Anywhere/GS2/Results/R2_bias/gene_div/", c, ".gene_div.txt")
  gene_div <- read.table(div_file, header=T, sep = " ", fill=TRUE, stringsAsFactors = F)
  gene_div$geneid <- str_replace(gene_div$geneid, ".aln", "")

  perm_dist <- merge(CC_perm, gene_div, 'geneid')
  
  perm_dist$Cellular_Component <- "Metabolic"
  perm_dist$Cellular_Component[which(CC_perm$CC == "membrane")] <- "Membrane"
  perm_dist$Cellular_Component[which(CC_perm$CC == "organelle")] <- "Organelle"
  perm_dist$Cellular_Component <- as.factor(perm_dist$Cellular_Component)
  perm_dist$cluster <- c
  
  all_perm_dists <- rbind(all_perm_dists, perm_dist)
  
  #### Testing #####
  # test for significant differences between the different phylum
  CC_test <- kruskal.test(variance ~ Cellular_Component, data = perm_dist)
  kruskal_pvals <- c(kruskal_pvals, format.pval(CC_test$p.value))
  
  # if  Kruskal-Wallis test shows siginificance, use pairwise wilcox test
  wilcox_pvals[[i]] <- pairwise.wilcox.test(perm_dist$variance, perm_dist$Cellular_Component)
  
  
}


########## variance in organelle genes ############
group_by(all_perm_dists,Cellular_Component) %>%
  dplyr::summarise(
    count = n(),
    mean = mean(variance, na.rm = TRUE),
    sd = sd(variance, na.rm = TRUE),
    median = median(variance, na.rm = TRUE),
    IQR = IQR(variance, na.rm = TRUE)
  )

all_perm_dists$log_var <- log(all_perm_dists$variance)
all_perm_dists$cluster <- as.factor(all_perm_dists$cluster)
levels(all_perm_dists$cluster) <- clusters


##### testing variance ####
kruskal_pvals <- as.numeric(kruskal_pvals)
kruskal_pvals <- p.adjust(kruskal_pvals, method="fdr")
kruskal_pvals <- round(kruskal_pvals, digits = 3)
length(kruskal_pvals[which(kruskal_pvals<0.05)])
length(kruskal_pvals)
kruskal_df <- as.data.frame(kruskal_pvals)
rownames(kruskal_df) <- clusters
clusters[which(kruskal_pvals<0.05)]

### wilcox p values in df ####
names(wilcox_pvals) <- clusters
# check wilcox result for clusters that are significant in kruskal wallis test
wp <- as.data.frame(matrix(nrow=length(clusters[which(kruskal_pvals<0.05)]), ncol=3))
names(wp) <- c("mem_NA","mem_org", "org_NA")
rownames(wp) <- clusters[which(kruskal_pvals<0.05)]
wp["C1","mem_NA"] <- wilcox_pvals$C1$p.value[1]
wp["C1","mem_org"] <- wilcox_pvals$C1$p.value[2]
wp["C1","org_NA"] <- wilcox_pvals$C1$p.value[4]

wp["C2","mem_NA"] <- wilcox_pvals$C2$p.value[1]
wp["C2","mem_org"] <- wilcox_pvals$C2$p.value[2]
wp["C2","org_NA"] <- wilcox_pvals$C2$p.value[4]

wp["C3","mem_NA"] <- wilcox_pvals$C3$p.value[1]
wp["C3","mem_org"] <- wilcox_pvals$C3$p.value[2]
wp["C3","org_NA"] <- wilcox_pvals$C3$p.value[4]

wp["C4","mem_NA"] <- wilcox_pvals$C4$p.value[1]
wp["C4","mem_org"] <- wilcox_pvals$C4$p.value[2]
wp["C4","org_NA"] <- wilcox_pvals$C4$p.value[4]

wp["C5","mem_NA"] <- wilcox_pvals$C5$p.value[1]
wp["C5","mem_org"] <- wilcox_pvals$C5$p.value[2]
wp["C5","org_NA"] <- wilcox_pvals$C5$p.value[4]

wp["C6","mem_NA"] <- wilcox_pvals$C6$p.value[1]
wp["C6","mem_org"] <- wilcox_pvals$C6$p.value[2]
wp["C6","org_NA"] <- wilcox_pvals$C6$p.value[4]

wp["C7","mem_NA"] <- wilcox_pvals$C7$p.value[1]
wp["C7","mem_org"] <- wilcox_pvals$C7$p.value[2]
wp["C7","org_NA"] <- wilcox_pvals$C7$p.value[4]

wp["C8","mem_NA"] <- wilcox_pvals$C8$p.value[1]
wp["C8","mem_org"] <- wilcox_pvals$C8$p.value[2]
wp["C8","org_NA"] <- wilcox_pvals$C8$p.value[4]

wp["C9","mem_NA"] <- wilcox_pvals$C9$p.value[1]
wp["C9","mem_org"] <- wilcox_pvals$C9$p.value[2]
wp["C9","org_NA"] <- wilcox_pvals$C9$p.value[4]

wp["C11","mem_NA"] <- wilcox_pvals$C11$p.value[1]
wp["C11","mem_org"] <- wilcox_pvals$C11$p.value[2]
wp["C11","org_NA"] <- wilcox_pvals$C11$p.value[4]

wp["C12","mem_NA"] <- wilcox_pvals$C12$p.value[1]
wp["C12","mem_org"] <- wilcox_pvals$C12$p.value[2]
wp["C12","org_NA"] <- wilcox_pvals$C12$p.value[4]

wp["C13","mem_NA"] <- wilcox_pvals$C13$p.value[1]
wp["C13","mem_org"] <- wilcox_pvals$C13$p.value[2]
wp["C13","org_NA"] <- wilcox_pvals$C13$p.value[4]

wp["C14","mem_NA"] <- wilcox_pvals$C14$p.value[1]
wp["C14","mem_org"] <- wilcox_pvals$C14$p.value[2]
wp["C14","org_NA"] <- wilcox_pvals$C14$p.value[4]

wp["C15","mem_NA"] <- wilcox_pvals$C15$p.value[1]
wp["C15","mem_org"] <- wilcox_pvals$C15$p.value[2]
wp["C15","org_NA"] <- wilcox_pvals$C15$p.value[4]

wp["C16","mem_NA"] <- wilcox_pvals$C16$p.value[1]
wp["C16","mem_org"] <- wilcox_pvals$C16$p.value[2]
wp["C16","org_NA"] <- wilcox_pvals$C16$p.value[4]

wp["C17","mem_NA"] <- wilcox_pvals$C17$p.value[1]
wp["C17","mem_org"] <- wilcox_pvals$C17$p.value[2]
wp["C17","org_NA"] <- wilcox_pvals$C17$p.value[4]

wp["C18","mem_NA"] <- wilcox_pvals$C18$p.value[1]
wp["C18","mem_org"] <- wilcox_pvals$C18$p.value[2]
wp["C18","org_NA"] <- wilcox_pvals$C18$p.value[4]

wp["C20","mem_NA"] <- wilcox_pvals$C20$p.value[1]
wp["C20","mem_org"] <- wilcox_pvals$C20$p.value[2]
wp["C20","org_NA"] <- wilcox_pvals$C20$p.value[4]

wp["C22","mem_NA"] <- wilcox_pvals$C22$p.value[1]
wp["C22","mem_org"] <- wilcox_pvals$C22$p.value[2]
wp["C22","org_NA"] <- wilcox_pvals$C22$p.value[4]

wp["C24","mem_NA"] <- wilcox_pvals$C24$p.value[1]
wp["C24","mem_org"] <- wilcox_pvals$C24$p.value[2]
wp["C24","org_NA"] <- wilcox_pvals$C24$p.value[4]

wp["C25","mem_NA"] <- wilcox_pvals$C25$p.value[1]
wp["C25","mem_org"] <- wilcox_pvals$C25$p.value[2]
wp["C25","org_NA"] <- wilcox_pvals$C25$p.value[4]

wp["C29","mem_NA"] <- wilcox_pvals$C29$p.value[1]
wp["C29","mem_org"] <- wilcox_pvals$C29$p.value[2]
wp["C29","org_NA"] <- wilcox_pvals$C29$p.value[4]

wp["C30","mem_NA"] <- wilcox_pvals$C30$p.value[1]
wp["C30","mem_org"] <- wilcox_pvals$C30$p.value[2]
wp["C30","org_NA"] <- wilcox_pvals$C30$p.value[4]

wp["C31","mem_NA"] <- wilcox_pvals$C31$p.value[1]
wp["C31","mem_org"] <- wilcox_pvals$C31$p.value[2]
wp["C31","org_NA"] <- wilcox_pvals$C31$p.value[4]

wp["C32","mem_NA"] <- wilcox_pvals$C32$p.value[1]
wp["C32","mem_org"] <- wilcox_pvals$C32$p.value[2]
wp["C32","org_NA"] <- wilcox_pvals$C32$p.value[4]

wp["C37","mem_NA"] <- wilcox_pvals$C37$p.value[1]
wp["C37","mem_org"] <- wilcox_pvals$C37$p.value[2]
wp["C37","org_NA"] <- wilcox_pvals$C37$p.value[4]

wp["C38","mem_NA"] <- wilcox_pvals$C38$p.value[1]
wp["C38","mem_org"] <- wilcox_pvals$C38$p.value[2]
wp["C38","org_NA"] <- wilcox_pvals$C38$p.value[4]

wp["C41","mem_NA"] <- wilcox_pvals$C41$p.value[1]
wp["C41","mem_org"] <- wilcox_pvals$C41$p.value[2]
wp["C41","org_NA"] <- wilcox_pvals$C41$p.value[4]
wp$cluster <- rownames(wp)

molten_wp <- melt(wp)
molten_wp$value <- p.adjust(molten_wp$value, method="BH")
molten_wp$value <- round(molten_wp$value, digits = 3)
nrow(molten_wp[which(molten_wp$variable=="org_NA" & molten_wp$value<0.05),])
nrow(molten_wp[which(molten_wp$variable=="mem_org" & molten_wp$value<0.05),])
nrow(molten_wp[which(molten_wp$variable=="mem_NA" & molten_wp$value<0.05),])

#### plotting ######
sign_clusters <- clusters[which(kruskal_pvals<0.05)]
all_perm_dists$significance <- ""
for(i in 1:length(sign_clusters))
{
  all_perm_dists$significance[which(all_perm_dists$cluster==sign_clusters[i])] <- "*"
}

# Plot variance in one big plot
p <- ggboxplot(all_perm_dists, x = "Cellular_Component", y = "log_var", na.rm=T) + 
  geom_boxplot(aes(fill=Cellular_Component)) +
  theme_bw() +
  removeGrid() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  facet_wrap( ~ cluster, scales="free", ncol = 5, labeller = labeller(cluster = c("C1"="C1***","C2"="C2***","C3"="C3**","C4"="C4**","C5"="C5***","C6"="C6**","C7"="C7***","C8"="C8**","C9"="C9***",
                                                                                  "C11"="C11*","C12"="C12***","C13"="C13***","C14"="C14***","C15"="C15***","C16"="C16***","C17"="C17***","C18"="C18***",
                                                                                  "C20"="C20***","C21"="C21","C22"="C22***","C23"="C23","C24"="C24***","C25"="C25**","C28"="C28","C29"="C29*",
                                                                                  "C30"="C30*","C31"="C31***","C32"="C32***","C34"="C34","C37"="C37***","C38"="C38**","C41"="C41*"))) + 
  xlab("gene groups") + ylab("log mean pi") + 
  theme(legend.position="top")
ggsave(plot=p,width=18,height=22,units="cm",dpi=300,filename="/Users/maloj/HCP_Anywhere/GS2/Results/R2_bias/variance_gene_groups.pdf", useDingbats=FALSE)
p

########## one plot R2 box ################## 

# exclude negative R2 values
nrow(all_perm_dists)
nrow(all_perm_dists[which(all_perm_dists$R2<0),])
all_perm_dists <- all_perm_dists[which(!all_perm_dists$R2<0),]
all_perm_dists$cluster <- factor(all_perm_dists$cluster, levels = c("C1","C2","C3","C4","C5","C6","C7","C8","C9","C11","C12","C13","C14","C15","C16","C17","C18","C20","C21","C22","C23","C24","C25","C28","C29","C30","C31","C32","C34","C37","C38","C41"))

# Plot R2 in one big plot
p <- ggboxplot(all_perm_dists, x = "Cellular_Component", y = "R2", color = "Cellular_Component", na.rm=T) +
  theme_bw() +
  removeGrid()  + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  facet_wrap( ~ cluster, scales="free", ncol = 6, labeller = labeller(cluster = c("C1"="C1","C2"="C2","C3"="C3*","C4"="C4","C5"="C5*","C6"="C6","C7"="C7","C8"="C8*","C9"="C9",
                                                                                         "C11"="C11","C12"="C12","C13"="C13**","C14"="C14***","C15"="C15","C16"="C16","C17"="C17","C18"="C18",
                                                                                         "C20"="C20","C21"="C21","C22"="C22","C23"="C23","C24"="C24***","C25"="C25","C28"="C28","C29"="C29",
                                                                                         "C30"="C30","C31"="C31","C32"="C32","C34"="C34***","C37"="C37**","C38"="C38***","C41"="C41"))) + 
  xlab("gene groups") +   ylab(expression(R^2 ~ "geography"))+ 
  theme(legend.position="top")
p
ggsave(plot=p,width=18,height=18,units="cm",dpi=300,filename="/Users/maloj/HCP_Anywhere/GS2/Results/R2_results/R2_boxplots.pdf", useDingbats=FALSE)



ggplot(all_perm_dists, aes(x=R2, y=variance)) + 
  geom_point()



