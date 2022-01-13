library(gplots)
library(RColorBrewer)
library(vegan)
library(stringr)
library(dplyr)
library(ggpubr)
library(tidyr)
#library(tidyverse)
library(reshape2)
library(viridis)
library(ggExtra)


# read data
setwd("~/HCP_Anywhere/GS2/Results/abundance")
# almeida rpkm matrix seperate
almeida <- read.table("almeida.rpkm.mat.tsv", header=T, row.names = 1)
names(almeida) <- sapply(strsplit(names(almeida), split="[.]"), head, 1)
# GS2 rpkm matrix seperate
GS2 <- read.table("GS2.rpkm.mat.tsv", header=T, row.names = 1)
names(GS2) <- sapply(strsplit(names(GS2), split="[.]"), head, 1)
# combined matrix
combined <- read.table("combined.rpkm.mat.tsv", header=T, row.names = 1)
names(combined) <- sapply(strsplit(names(combined), split="[.]"), head, 1)
# metadata
merged_meta <- read.delim("/Users/maloj/HCP_Anywhere/GS2/GS2_metadata/merged_metadata.txt", stringsAsFactors = F)
merged_meta$new_complete_name <- str_replace(merged_meta$new_complete_name, "-", "_")


#colSums(combined)
#colSums(almeida)

#### Aplha diversity ####
# test example
#data(BCI)
#H <- diversity(BCI)

#### Seperate mappings ####
# calculate shannon diversity index and save in df
div_almeida <- as.data.frame(diversity(t(almeida)))
names(div_almeida) <- "shannon_index"
div_almeida$sample_origin <- "Human gut"

div_GS2 <- as.data.frame(diversity(t(GS2)))
names(div_GS2) <- "shannon_index"
div_GS2$sample_origin <- NA
# insert regions from metadata in "sample_origin" column
for(i in 1:nrow(div_GS2))
{
  samplename <- rownames(div_GS2[i,])
  if(str_detect(samplename, "DTU"))
  {
    #print(samplename)
    div_GS2$sample_origin[i] <- merged_meta$Region[which(merged_meta$new_complete_name==samplename)]
  }
  else
  {
    samplename <- substring(samplename, 2)
    #print(samplename)
    div_GS2$sample_origin[i] <- merged_meta$Region[which(merged_meta$unique_sample_id==samplename)][1]
  }
}

div_seperate <- rbind(div_GS2, div_almeida)

ggboxplot(div_seperate, x = "sample_origin", y = "shannon_index", color = "sample_origin") + 
  theme(axis.text.x=element_text(color = "black", size=11, angle=30, vjust=.8, hjust=0.8)) + 
  labs(title = "Shannon diversity index according to sample origin") + ylab("Shannon diversity index") + xlab("Sample origin")




#### Combined mappings ####
div_combined <- as.data.frame(diversity(t(combined)))
names(div_combined) <- "shannon_index"
div_combined$sample_origin <- NA
# insert regions from metadata in "sample_origin" column
for(i in 1:nrow(div_combined))
{
  samplename <- rownames(div_combined[i,])
  if(str_detect(samplename, "DTU"))
  {
    #print(samplename)
    div_combined$sample_origin[i] <- merged_meta$Region[which(merged_meta$new_complete_name==samplename)]
  }
  else if(str_detect(samplename, "ERR"))
  {
    div_combined$sample_origin[i] <- "Human gut"
  }
  else if(str_detect(samplename, "SRR"))
  {
    div_combined$sample_origin[i] <- "Human gut"
  }
  else
  {
    samplename <- substring(samplename, 2)
    #print(samplename)
    div_combined$sample_origin[i] <- merged_meta$Region[which(merged_meta$unique_sample_id==samplename)][1]
  }
}

# plot
pdf("/Users/maloj/HCP_Anywhere/GS2/Results/abundance/shannon.boxplot.pdf", width=4, height=4)
ggboxplot(div_combined, x = "sample_origin", y = "shannon_index", color = "sample_origin")  + 
  ylab("") + xlab("") +
  theme_bw() +
  removeGrid() + 
  theme(axis.text.x=element_text(color = "black", size=11, angle=30, vjust=.8, hjust=0.8)) +
  theme(legend.position = "None")
dev.off()

# get summaries
group_by(div_combined, sample_origin) %>%
  summarise(
    count = n(),
    mean = mean(shannon_index, na.rm = TRUE),
    sd = sd(shannon_index, na.rm = TRUE),
    median = median(shannon_index, na.rm = TRUE),
    IQR = IQR(shannon_index, na.rm = TRUE)
  )




#### Beta diversity ####
# Calculate a distance matrix
combined.beta.dist <- vegdist(t(combined), method = "bray") # you have to transpose the dataset to get the genera as rows
# Cluster the distance matrix.
combined.hc <- hclust(as.dist(combined.beta.dist))
# Process and melt the distance matrix
combined.beta.dist.long <- combined.beta.dist %>% as.matrix %>% melt %>%
  mutate(Var1 = as.character(Var1), Var2 = as.character(Var2))
# Plot the heatmap
combined.beta.dist.long %>% #as.matrix %>% .[vare.hc$order, vare.hc$order] %>% melt %>%
  ggplot(aes(x = Var1, y = Var2, fill = value)) + geom_tile() + scale_fill_viridis(direction = 1) +
  theme(axis.text.x = element_text(angle = 270, hjust = 0, vjust = 0.5))



# simple plot
heatmap(as.matrix(combined.beta.dist), keep.dendro=F)



# test example
data(varespec)

# Calculate a distance matrix
vare.dist <- vegdist(varespec)

# Cluster the distance matrix.
vare.hc <- hclust(as.dist(vare.dist))

# Process and melt the distance matrix
vare.dist.long <- vare.dist %>% as.matrix %>% melt %>%
  mutate(Var1 = as.character(Var1), Var2 = as.character(Var2))

vare.dist.long1 <- vare.dist.long %>% mutate(Var1 = factor(Var1, levels = unique(Var1)[vare.hc$order]))
vare.dist.long2 <- vare.dist.long %>% mutate(Var2 = factor(Var2, levels = unique(Var2)[vare.hc$order]))


# Plot the heatmap
vare.dist.long2 %>% #as.matrix %>% .[vare.hc$order, vare.hc$order] %>% melt %>%
  ggplot(aes(x = Var1, y = Var2, fill = value)) + geom_tile() + scale_fill_viridis(direction = 1) +
  theme(axis.text.x = element_text(angle = 270, hjust = 0, vjust = 0.5))
                                   





col.clus <- hclust(combined.beta.dist, "aver")

# Do average linkage hierarchical clustering. 
#Other options are 'complete' or 'single'. You'll need to choose the one that best fits the needs of your situation and your data.
row.clus <- hclust(combined.beta.dist, "aver")

# make the heatmap with Rowv = as.dendrogram(row.clus)
heatmap(as.matrix(combined.beta.dist), Rowv = as.dendrogram(row.clus), Colv = as.dendrogram(col.clus), margins = c(10, 3))


drows = vegdist(combined, method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = FALSE)
dcols =vegdist(t(combined), method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = FALSE)
pheatmap(combined, clustering_distance_rows = combined.beta.dist, clustering_distance_cols = combined.beta.dist)


#### PcoA ####
# First step is to calculate a distance matrix. 
# Here we use Bray-Curtis distance metric
dist <- vegdist(varespec,  method = "bray")

# PCoA is not included in vegan. 
# We will use the ape package instead
library(ape)
PCOA <- pcoa(dist)

# plot the eigenvalues and interpret
barplot(PCOA$values$Relative_eig[1:10])
# Can you also calculate the cumulative explained variance of the first 3 axes?

# Some distance measures may result in negative eigenvalues. In that case, add a correction:
PCOA <- pcoa(dist, correction = "cailliez")

# Plot your results
biplot.pcoa(PCOA)

# You see what`s missing? 
# Indeed, there are no species plotted on this biplot. 
# That's because we used a dissimilarity matrix (sites x sites) 
# as input for the PCOA function. 
# Hence, no species scores could be calculated. 
#However, we could work around this problem like this:
biplot.pcoa(PCOA, varespec)


