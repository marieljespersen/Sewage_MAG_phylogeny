library(reshape2)
library(ggplot2)
library(ggExtra)


########## look at p-values form permanova test of ASTRAL species trees ##########
# using pvals_df from  tree_permanova
pvals_df <-  read.table("/Users/maloj/HCP_Anywhere/GS2/Results/city_permanova_res.txt")
pvals_df$significance <- NA
pvals_df$significance[which(pvals_df$p_val<0.05)] <- "P-val < 0.05"
pvals_df$significance[which(pvals_df$p_val>=0.05)] <- "P-val >= 0.05"
pvals_df$cluster <- rownames(pvals_df)
# reove C26 due to bias 
pvals_df <- pvals_df[!(row.names(pvals_df) %in% c("C26")), ]

pdf("/Users/maloj/HCP_Anywhere/GS2/Results/R2_vs_genomes.pdf", width=6, height=4)
ggplot(pvals_df,aes(x=samples_in_trees, y=R2_model, color=significance)) + 
  theme_minimal() +
  xlab("Genomes in tree") + 
  ylab("Regional R2") +
  geom_text(aes(label=cluster),size = 5) +
  theme(legend.position="top") +
  theme_bw() +
  removeGrid() +
  theme(axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15)) +
  theme(legend.text=element_text(size=15),
        legend.title=element_text(size=15))
dev.off()

########## look at differences in R2 values between organelle and membrane genes ##########
# using memR2 and orgR2 from  ANOVA_R2_all_clusters
R2diff <- as.data.frame(cbind(memR2,naR2,orgR2))
rownames(R2diff) <- clusters
names(R2diff) <- c("membrane", "Other", "organelle")
R2diff <- cbind(R2diff, clusters)
R2diff$difference <- NA
R2diff$difference[which(R2diff$membrane-R2diff$organelle<0)] <- "Organelle larger"
R2diff$difference[which(R2diff$membrane-R2diff$organelle>0)] <- "Membrane larger"
R2diff$clusters <- factor(R2diff$clusters, levels = R2diff$clusters[order(R2diff$membrane)])
moltenR2 <- melt(R2diff)

# exclude C26 because of biased R2 values
moltenR2 <- moltenR2[which(moltenR2$clusters!="C26"),] 
# plot
#pdf(file="/Users/maloj/HCP_Anywhere/GS2/Results/Mean_R2_diff.pdf", width=6, height=4)
p <- ggplot(moltenR2) + 
  geom_point(
    aes(x = clusters, y = value, fill = variable), 
    size=4, pch = 21, colour = alpha("white", 0))+
  theme_bw() +
  geom_linerange(
    aes(x = clusters, ymin = membrane, ymax = organelle, colour = difference), 
    data = spread(moltenR2, variable, value)) +
  scale_color_brewer(palette="Dark2") +
  xlab("") + ylab("Mean R2") +
  theme(text=element_text(size=11)) +
  theme(axis.text.x = element_text(size = 11,angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size = 11)) +     
  geom_text(aes(x = "C3", y = 0.26, label = "*")) +     
  geom_text(aes(x = "C5", y = 0.31, label = "*")) +     
  geom_text(aes(x = "C8", y = 0.33, label = "*")) +     
  geom_text(aes(x = "C13", y = 0.25, label = "**")) +     
  geom_text(aes(x = "C14", y = 0.83, label = "***")) +     
  geom_text(aes(x = "C24", y = 0.53, label = "***")) +     
  geom_text(aes(x = "C34", y = 0.41, label = "***")) +     
  geom_text(aes(x = "C37", y = 0.46, label = "**")) +     
  geom_text(aes(x = "C38", y = 0.6, label = "***")) +
  removeGrid()
ggsave(plot=p,width=6.5,height=3,dpi=200,filename="/Users/maloj/HCP_Anywhere/GS2/Results/Mean_R2_diff.pdf", useDingbats=FALSE)

#print(p)
#dev.off()
######### calculate % less regional variation in organelles ###########

percent_diff <- c()
for(i in 1:nrow(R2diff))
{
  meanR2memNA <- (R2diff[i,]$membrane+R2diff[i,]$Other)/2
  diff_mem_org <- meanR2memNA-R2diff[i,]$organelle
  cluster_p_diff <- diff_mem_org/meanR2memNA
  percent_diff <- c(percent_diff, cluster_p_diff)
  print(i)
}



