library(stringr)
library(RColorBrewer)
library(ggplot2)
library(ggExtra)
library(gridExtra)

setwd("/Users/maloj/HCP_Anywhere/GS2/Results/gtdbtk_tree")
bac_bins <- read.delim("tax_all_bins.tsv", stringsAsFactors = F)
alm_bins <- read.delim("clean.almeida.tax.txt", stringsAsFactors = F)
merged_meta <- read.delim("/Users/maloj/HCP_Anywhere/GS2/GS2_metadata/merged_metadata.txt", stringsAsFactors = F)
phy_list <- read.delim("/Users/maloj/HCP_Anywhere/GS2/Results/gtdbtk_tree/phyla.txt",  sep=" ", header=F)
phy_col <- read.delim("/Users/maloj/HCP_Anywhere/GS2/Results/gtdbtk_tree/phyla_colors.txt",  sep=" ", header=F)
phy_df <- as.data.frame(cbind(t(phy_col),t(phy_list)))
names(phy_df) <- c("color", "Phylum")


# OBS! Der er forskel på antallet af phylum her of i gtdb træet fordi archea ikke er med

gs_bins <- bac_bins$bin
bac_bins <- rbind(bac_bins, alm_bins)

bac_bins$Region <- "Human Gut"
bac_bins$Sample <- NA
for(b in gs_bins)
{
  if(str_detect(b, ".vambbin."))
  {
    sample <- strsplit(b, ".vambbin.")[[1]][1]
    bac_bins$Sample[which(bac_bins$bin==b)] <- sample
    bac_bins$Region[which(bac_bins$bin==b)] <- merged_meta$Region[which(merged_meta$vamb_sample==sample)]
  }
  else if(str_detect(b, ".metabin."))
  {
    sample <- strsplit(b, ".metabin.")[[1]][1]
    bac_bins$Sample[which(bac_bins$bin==b)] <- sample
    bac_bins$Region[which(bac_bins$bin==b)] <- merged_meta$Region[which(merged_meta$new_complete_name==sample)]
  }
}

# create color vector
col_vector <- as.character(phy_df$color[order(phy_df$Phylum)])
bacbac <- bac_bins[which(bac_bins$superkingdom=="Bacteria"),]
bacbac$Region <- factor(bacbac$Region, levels = c("East Asia & Pacific", "Europe & Central Asia", "Latin America & Caribbean", "Middle East & North Africa", "North America","South Asia","Sub-Saharan Africa", "Human Gut"))


######### Plot alle ##############
bacbac$count <- 1
p <- ggplot(bacbac, aes(fill = phylum, y=count, x=Region)) + 
  geom_bar(position="fill", stat="identity") + 
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=12, angle=90, hjust=1, vjust=0.5),
        axis.text.y = element_blank()) +
  theme(axis.ticks.y = element_blank()) +
  scale_fill_manual(values = col_vector) +
  theme(legend.title=element_text(size=12), 
        legend.text=element_text(size=12)) +
  removeGrid()
print(p)
ggsave(plot=p,width=8,height=6,filename="/Users/maloj/HCP_Anywhere/GS2/Results/phylum_barplots/all_phylum_per_region.png", device = "png")
ggsave(plot=p,width=7,height=7,filename="/Users/maloj/HCP_Anywhere/GS2/Results/phylum_barplots/all_phylum_per_region_legend.png", device = "png")




######## plot kun phylum over threshold ###########

# udvælg phylum som er tilstede i en vis mængde
morethan20 <- 0
phylum_plot <- c()
for(p in unique(bac_bins$phylum[which(bac_bins$superkingdom=="Bacteria")]))
{
  count <- length(bac_bins$phylum[which(bac_bins$phylum==p)])
  if(count>20)
  {
    print(p)
    morethan20 <- morethan20 +1
    phylum_plot <- c(phylum_plot, p)
  }
}
print(morethan20)
print(phylum_plot)

bacbac$Phylum <- "others"
bacbac$Phylum[which(bacbac$phylum=="Proteobacteria")] <- "Proteobacteria"
bacbac$Phylum[which(bacbac$phylum=="Bacteroidetes")] <- "Bacteroidetes"
bacbac$Phylum[which(bacbac$phylum=="Actinobacteriota")] <- "Actinobacteriota"
bacbac$Phylum[which(bacbac$phylum=="Firmicutes")] <- "Firmicutes"
bacbac$Phylum[which(bacbac$phylum=="Synergistota")] <- "Synergistota"
bacbac$Phylum[which(bacbac$phylum=="Fusobacteriota")] <- "Fusobacteriota"
bacbac$Phylum[which(bacbac$phylum=="Chloroflexota")] <- "Chloroflexota"
bacbac$Phylum[which(bacbac$phylum=="Desulfobacterota")] <- "Desulfobacterota"
bacbac$Phylum[which(bacbac$phylum=="Verrucomicrobiota")] <- "Verrucomicrobiota"
bacbac$Phylum[which(bacbac$phylum=="Planctomycetota")] <- "Planctomycetota"
bacbac$Phylum[which(bacbac$phylum=="Acidobacteriota")] <- "Acidobacteriota"


# create color vector with same colors only for the selected phylum
col_df <- phy_df[which(phy_df$Phylum %in% phylum_plot),]
col_mat <- as.matrix(col_df)
col_mat <- rbind(col_mat, c("#AAAFAA", "others"))
col_df <- as.data.frame(col_mat)
col_vector <- as.character(col_df$color[order(col_df$Phylum)])

pdf("/Users/maloj/HCP_Anywhere/GS2/Results/phylum_barplots/phylum.barplot.pdf", width=4, height=4)
ggplot(bacbac, aes(fill = Phylum, y=count, x=Region)) + 
  geom_bar(position="fill", stat="identity") + 
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 50, hjust = 1,size = 7),
        axis.text.y = element_blank()) +
  theme(axis.ticks.y = element_blank()) +
  scale_fill_manual(values = col_vector) +
  removeGrid() +
  theme(legend.title=element_text(size=7), legend.text=element_text(size=6))
dev.off()
  
############ Test difference between human and GS ###########
bac_bins$phy_val <- as.numeric(factor(bac_bins$Phylum))
alm_phy <- bac_bins$phy_val[which(bac_bins$Region=="Human Gut")]
GS2_phy <- bac_bins$phy_val[which(bac_bins$Region!="Human Gut")]
ks.test(GS2_phy, alm_phy)


