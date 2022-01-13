library("FactoMineR")
library("factoextra")
library("compositions")
library(stringr)
library(dplyr)
library(ggExtra)

setwd("/Users/maloj/HCP_Anywhere/GS2/Results/abundance")
read_counts <- read.table("combined.read_count.mat.tsv", row.names=1, stringsAsFactors = F, header=T)
names(read_counts) <- str_remove_all(names(read_counts), ".coverage.tsv")


######### get tax table ##########
setwd("/Users/maloj/HCP_Anywhere/GS2/Results/gtdbtk_tree")
bac_bins <- read.delim("tax_all_bins.tsv", stringsAsFactors = F)
alm_bins <- read.delim("clean.almeida.tax.txt", stringsAsFactors = F)
merged_meta <- read.delim("/Users/maloj/HCP_Anywhere/GS2/GS2_metadata/merged_metadata.txt", stringsAsFactors = F)

# OBS! Der er forskel på antallet af phylum her of i gtdb træet fordi archea ikke er med

gs_bins <- bac_bins$bin
bac_bins <- rbind(bac_bins, alm_bins)

bac_bins$Region <- "Human Gut"
bac_bins$Sample <- NA
bac_bins$ID <- NA
for(b in gs_bins)
{
  if(str_detect(b, ".vambbin."))
  {
    sample <- strsplit(b, ".vambbin.")[[1]][1]
    bac_bins$Sample[which(bac_bins$bin==b)] <- sample
    bac_bins$Region[which(bac_bins$bin==b)] <- merged_meta$Region[which(merged_meta$vamb_sample==sample)]
    bac_bins$ID[which(bac_bins$bin==b)] <- merged_meta$unique_sample_id[which(merged_meta$vamb_sample==sample)]
  }
  else if(str_detect(b, ".metabin."))
  {
    sample <- strsplit(b, ".metabin.")[[1]][1]
    bac_bins$Sample[which(bac_bins$bin==b)] <- sample
    bac_bins$Region[which(bac_bins$bin==b)] <- merged_meta$Region[which(merged_meta$new_complete_name==sample)]
    bac_bins$ID[which(bac_bins$bin==b)] <- merged_meta$unique_sample_id[which(merged_meta$new_complete_name==sample)]
  }
}

bac_bins$bin <- str_replace_all(bac_bins$bin,"[.]", "_")

read_counts$phylum <- NA
read_counts$genus <- NA
for(i in 1:nrow(read_counts))
{
  bin <- row.names(read_counts[i,])
  if(bin %in% bac_bins$bin)
  {
    read_counts$phylum[i] <- bac_bins$phylum[which(bac_bins$bin==bin)]
    read_counts$genus[i] <- bac_bins$genus[which(bac_bins$bin==bin)]
  }
}

#length(read_counts$phylum[which(is.na(read_counts$phylum))])
#length(read_counts$phylum)
#table(read_counts$phylum)

############### PCA plot on species level ##########

genus_df <- read_counts %>%
  group_by(genus) %>%
  summarise_if(
    is.numeric,
    sum,
    na.rm = TRUE
  )

genus_df <- genus_df[1:562,]# removing NA count
row.names(genus_df) <- genus_df$genus
genus_with_meta <- as.data.frame(t(genus_df[,-1]))
names(genus_with_meta) <- genus_df$genus

genus_with_meta$region <- NA
for(i in 1:nrow(genus_with_meta))
{
  sample <- rownames(genus_with_meta)[i]
  if(str_detect(sample,"X2"))
  {
    sample <- sub(".", "", sample)
    genus_with_meta$region[i] <- merged_meta$Region[which(merged_meta$unique_sample_id==sample)][1]
  }
  else if(sample %in% merged_meta$vamb_sample)
  {
    genus_with_meta$region[i] <- merged_meta$Region[which(merged_meta$vamb_sample==sample)]
  }
  else
  {
    genus_with_meta$region[i] <- "Human gut"
  }
}

no_zeros <- exp(digamma( genus_with_meta[1:562]+1))
#no_zeros$region <- genus_with_meta$region

#close compostitions in order to compare between samples
closeure <- function(x) { 100/sum(x)*x }
closed_counts <- apply(t(no_zeros),2,closeure)


geomean <- apply(t(closed_counts), 2, geometricmean) # Bliver nul for alle??
var <- variation(acomp(t(closed_counts)))
tot_var <- (sum(var)/2)/562 # number of phylum
scaled_counts <- no_zeros^1/sqrt(tot_var)


scaled_centered_counts <- scale(scaled_counts, center=geomean)
clr_counts <- clr(scaled_centered_counts)

count_PC <- FactoMineR::PCA(clr_counts,scale.unit = FALSE)
fviz_eig(count_PC, addlabels = TRUE, ylim = c(0, 100))
ind.p <- fviz_pca_ind(count_PC, geom.ind = "point", pointshape=20, col.ind = genus_with_meta$region)
ggpubr::ggpar(ind.p,
              title = "Principal Component Analysis",
              subtitle = "Phylum abundance",
              legend.title = "Sample origin", legend.position = "top",
              xlab = "PC1 (21.7%)", ylab = "PC2 (12.1%)",
              addEllipses = TRUE,
              ggtheme = theme_gray(),
)

fviz_pca_biplot(count_PC, 
                col.ind = genus_with_meta$region, 
                repel = TRUE,
                label = "none",
                col.var = "black",
                legend.title = "Countries",
                addEllipses = TRUE,
                xlab = "PC1 (21.7%)", ylab = "PC2 (12.1%)",
                ggtheme = theme_minimal(),
                pointshape=20
) 

fviz_pca_ind(ind.p, label="none", habillage=genus_with_meta$region,
             addEllipses=TRUE)


############### PCA plot on genus level ##################

genus_df <- read_counts %>%
  group_by(genus) %>%
  summarise_if(
    is.numeric,
    sum,
    na.rm = TRUE
  )

genus_df <- genus_df[1:562,]# removing NA count
row.names(genus_df) <- genus_df$genus
genus_with_meta <- as.data.frame(t(genus_df[,-1]))
names(genus_with_meta) <- genus_df$genus

genus_with_meta$region <- NA
for(i in 1:nrow(genus_with_meta))
{
  sample <- rownames(genus_with_meta)[i]
  if(str_detect(sample,"X2"))
  {
    sample <- sub(".", "", sample)
    genus_with_meta$region[i] <- merged_meta$Region[which(merged_meta$unique_sample_id==sample)][1]
  }
  else if(sample %in% merged_meta$vamb_sample)
  {
    genus_with_meta$region[i] <- merged_meta$Region[which(merged_meta$vamb_sample==sample)]
  }
  else
  {
    genus_with_meta$region[i] <- "Human gut"
  }
}

no_zeros <- exp(digamma( genus_with_meta[1:562]+1))
#no_zeros$region <- genus_with_meta$region

#close compostitions in order to compare between samples
closeure <- function(x) { 100/sum(x)*x }
closed_counts <- apply(t(no_zeros),2,closeure)


geomean <- apply(t(closed_counts), 2, geometricmean) # Bliver nul for alle??
var <- variation(acomp(t(closed_counts)))
tot_var <- (sum(var)/2)/562 # number of phylum
scaled_counts <- no_zeros^1/sqrt(tot_var)


scaled_centered_counts <- scale(scaled_counts, center=geomean)
clr_counts <- clr(scaled_centered_counts)

count_PC <- FactoMineR::PCA(clr_counts,scale.unit = FALSE)
fviz_eig(count_PC, addlabels = TRUE, ylim = c(0, 100))
ind.p <- fviz_pca_ind(count_PC, geom.ind = "point", pointshape=20, col.ind = genus_with_meta$region)
ggpubr::ggpar(ind.p,
              title = "Principal Component Analysis",
              subtitle = "Phylum abundance",
              legend.title = "Sample origin", legend.position = "top",
              xlab = "PC1 (21.7%)", ylab = "PC2 (12.1%)",
              addEllipses = TRUE,
              ggtheme = theme_gray(),
)

fviz_pca_biplot(count_PC, 
                col.ind = genus_with_meta$region, 
                repel = TRUE,
                label = "none",
                col.var = "black",
                legend.title = "Countries",
                addEllipses = TRUE,
                xlab = "PC1 (21.7%)", ylab = "PC2 (12.1%)",
                ggtheme = theme_minimal(),
                pointshape=20
) 

fviz_pca_ind(ind.p, label="none", habillage=genus_with_meta$region,
             addEllipses=TRUE)

################# PCA on phylum level ###################
phylum_df <- read_counts %>%
  group_by(phylum) %>%
  summarise_if(
    is.numeric,
    sum,
    na.rm = TRUE
  )

phylum_df <- phylum_df[1:36,]# removing NA count
row.names(phylum_df) <- phylum_df$phylum
phy_with_meta <- as.data.frame(t(phylum_df[,-1]))
names(phy_with_meta) <- phylum_df$phylum

phy_with_meta$region <- NA
for(i in 1:nrow(phy_with_meta))
{
  sample <- rownames(phy_with_meta)[i]
  if(str_detect(sample,"X2"))
  {
    sample <- sub(".", "", sample)
    phy_with_meta$region[i] <- merged_meta$Region[which(merged_meta$unique_sample_id==sample)][1]
  }
  else if(sample %in% merged_meta$vamb_sample)
  {
    phy_with_meta$region[i] <- merged_meta$Region[which(merged_meta$vamb_sample==sample)]
  }
  else
  {
    phy_with_meta$region[i] <- "Human gut"
  }
}

#close compostitions in order to compare between samples
closeure <- function(x) { 100/sum(x)*x }

no_zeros <- exp(digamma(phy_with_meta[,1:36]+1))
closed_counts <- apply(t(no_zeros),2,closeure)


geomean <- apply(t(closed_counts), 2, geometricmean) # Bliver nul for alle??
var <- variation(acomp(t(closed_counts)))
tot_var <- (sum(var)/2)/36 # number of phylum
scaled_counts <- no_zeros^1/sqrt(tot_var)


scaled_centered_counts <- scale(scaled_counts, center=geomean)
clr_counts <- clr(scaled_centered_counts)

# forsøg med factoextra
# http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/
#count_PC <- FactoMineR::PCA(clr_counts,scale.unit = FALSE)
fviz_eig(count_PC, addlabels = TRUE, ylim = c(0, 100))
ind.p <- fviz_pca_ind(count_PC, geom.ind = "point", pointshape=20, col.ind = phy_with_meta$region,habillage=as.factor(phy_with_meta$region))
ggpubr::ggpar(ind.p,
              title = "Principal Component Analysis",
              subtitle = "Phylum abundance",
              legend.title = "Sample origin", legend.position = "top",
              xlab = "PC1 (18.9%)", ylab = "PC2 (12.2%)",
              habillage=as.factor(phy_with_meta$region),
              addEllipses = TRUE,
              ggtheme = theme_gray(),
)

################ Dette plot er med i artiklen!! ###################
pdf(file="/Users/maloj/HCP_Anywhere/GS2/Results/abundance/PCA/PCA.pdf", width=6.5, height=4)
fviz_pca_ind(count_PC, 
             label="none",
             col.ind = as.factor(phy_with_meta$region),
             addEllipses = TRUE,
             title = "",
             xlab = "PC1 (18.9%)", ylab = "PC2 (12.2%)",
             ggtheme = theme_bw()
             ) +
  removeGrid() +
  theme(text = element_text(size = 15))     
dev.off()

###############
fviz_pca_biplot(count_PC, 
                col.ind = phy_with_meta$region, 
                label = "var",
                col.var = "black", 
                repel = TRUE,
                legend.title = "Countries",
                addEllipses = TRUE,
                xlab = "PC1 (18.9%)", ylab = "PC2 (12.2%)",
                ggtheme = theme_minimal(),
                pointshape=20
) 



########### PCA 50 most prevalent phylums #############

# create df of 50 most prevalent phylums
prev_all <- c()
for(i in 1:ncol(phy_with_meta))
{
  prevalence <- length(phy_with_meta[,i][which(phy_with_meta[,i]!=0)])
  prev_all <- c(prev_all, prevalence)
}

phy_50 <- phy_with_meta[,which(prev_all>50)]
test <- exp(digamma(phy_50[1:22]+1))
test$region <- phy_50$region

#close compostitions in order to compare between samples
closeure <- function(x) { 100/sum(x)*x }
closed_counts <- apply(t(test[1:22]),2,closeure)


geomean <- apply(t(closed_counts), 2, geometricmean) # Bliver nul for alle??
var <- variation(acomp(t(closed_counts)))
tot_var <- (sum(var)/2)/22 # number of phylum
scaled_counts <- test[1:22]^1/sqrt(tot_var)


scaled_centered_counts <- scale(scaled_counts, center=geomean)
clr_counts <- clr(scaled_centered_counts)

library(ChemometricsWithR)
count_PC <- ChemometricsWithR::PCA(clr_counts)

# Plotting result from PCA
# Show score and loadings plot coloured according to country
#par(mfrow=c(1,2))
#par(mfrow=c(1,1))
ChemometricsWithR::scoreplot(count_PC, main = "Scores",col = as.integer(as.factor(phy_50$region)))
legend("topleft", inset=.05, c("cod", "dnk", "ind", "usa"), lwd=1, col=1:4, lty=1, cex = 0.5)

loadingplot(count_PC, show.names = TRUE, main='Loadings',xlim=c(-0.7,0.7), ylim=c(-0.7,0.7))
plot(count_PC$var, cex=2, main="Scree type plot: variances", xlab='Number of Components', ylab='Eigenvalue')
lines(count_PC$var)

#colour according to year
scoreplot(AMR_PC, main = "Scores",col = as.integer(as.factor(AMR$YEAR)))
legend("topleft", inset=.05, c("2016", "2017", "2018"), lwd=1, col=2016:2018, lty=1, cex = 0.5)
loadingplot(AMR_PC, show.names = TRUE, main='Loadings')



# forsøg med factoextra
# http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/
count_PC <- FactoMineR::PCA(clr_counts,scale.unit = FALSE)
fviz_eig(count_PC, addlabels = TRUE, ylim = c(0, 100))
ind.p <- fviz_pca_ind(count_PC, geom.ind = "point", pointshape=20, col.ind = phy_50$region)
ggpubr::ggpar(ind.p,
              title = "Principal Component Analysis",
              subtitle = "Phylum abundance",
              legend.title = "Sample origin", legend.position = "top",
              xlab = "PC1 (21.7%)", ylab = "PC2 (12.1%)",
              addEllipses = TRUE,
              ggtheme = theme_gray(),
)

fviz_pca_biplot(count_PC, 
                col.ind = phy_50$region, 
                label = "var",
                col.var = "black", 
                repel = TRUE,
                legend.title = "Countries",
                addEllipses = TRUE,
                xlab = "PC1 (21.7%)", ylab = "PC2 (12.1%)",
                ggtheme = theme_minimal(),
                pointshape=20
) 





