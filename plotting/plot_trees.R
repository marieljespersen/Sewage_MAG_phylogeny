library(treeio)
library(ggplot2)
library(ggtree)
library(stringr)
library(network)
library(RColorBrewer) 
library(gridExtra)

### read in data
setwd("/Users/maloj/HCP_Anywhere/GS2/Results/trees")
#treeFiles <- lapply(Sys.glob("C*.astral.branch_length.clean.tree"), read.tree) # read as C1, C10, C11 ...
#treeList <- Sys.glob("C*.astral.branch_length.clean.tree")
metadata <- read.delim("/Users/maloj/HCP_Anywhere/GS2/GS2_metadata/gs2_metadata.txt", stringsAsFactors = F)
metadata$new_complete_name <- str_replace(metadata$new_complete_name, "-", "_")
who_regions <- read.csv("/Users/maloj/HCP_Anywhere/GS2/GS2_metadata/world_bank.csv", stringsAsFactors = F)
merged_meta <- merge(metadata, who_regions, by.x = "country_alpha3", by.y = "Country.Code", sort = F)
cluster_no <- c(1,2,3,4,5,6,7,8,9,11,12,13,14,15,16,17,18,20,21,22,23,24,25,26,28,29,30,31,32,34,37,38,41)
clusters <- c()


######## create color palette to all trees #########
levels(as.factor(merged_meta$Region))
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n = length(levels(as.factor(merged_meta$Region)))
cols = gg_color_hue(n)
color_df <- as.data.frame(cols)
color_df$reg <- levels(as.factor(merged_meta$Region))


########### functions #############
#  function to correct node names in tree
correct.treenames <- function(t) {
  
  bins <- t$tip.label
  for(j in 1:length(bins))
  {
    t$tip.label[j] <- str_replace_all(bins[j], "-", "N")
    if(str_detect(t$tip.label[j], "_vambbin_"))
    {
      sample <- strsplit(t$tip.label[j], "_vambbin_")[[1]][1]
    }
    else if(str_detect(t$tip.label[j], "_metabin_"))
    {
      sample <- strsplit(t$tip.label[j], "_metabin_")[[1]][1]
    }
    #print(sample)
    t$tip.label[j] <- sample
  }
  return(t)
}

# function to drop tips from tree if multiple tips from same loaction
drop.multiple.locations <- function(t, md, c) 
{
  dRep_path <- paste("/Users/maloj/HCP_Anywhere/GS2/Results/dRep_GS2/clusters/dRep_scores/",c,".Sdb.csv", sep="")
  dRep_score <- read.csv(dRep_path)
  
  # Removing ".fna" from genome name
  dRep_score$genome <- str_replace(dRep_score$genome, ".fna",  "")
  
  # Creating column with sample names
  for(j in 1:nrow(dRep_score))
  {
    if(str_detect(dRep_score$genome[j], "_vambbin_"))
    {
      sample <- strsplit(dRep_score$genome[j], "_vambbin_")[[1]][1]
    }
    else if(str_detect(dRep_score$genome[j], "_metabin_"))
    {
      sample <- strsplit(dRep_score$genome[j], "_metabin_")[[1]][1]
    }
    #print(sample)
    dRep_score$sample[j] <- sample
  }
  
  
  dupl_cities <- unique(md$city[duplicated(md$city)])
  for(c in dupl_cities)
  {
    city_samples <- md$new_complete_name[which(md$city==c)]
    sub_drep <- dRep_score[dRep_score$sample %in% city_samples,]
    samples_to_drop <- sub_drep$sample[which(!sub_drep$score==max(sub_drep$score))]
    t <- drop.tip(t, samples_to_drop)
  }
  return(t)
}

###### plotting #########
#pdf(file = "circular_trees.pdf")

# cluster <- "C5"
#plot_list <- list()
for(i in 1:length(cluster_no))
{
  # read in tree
  cluster <- paste0("C",cluster_no[i])
  clusters <- c(clusters, cluster)
  #cluster <- strsplit(treeList[i], '[.]')[[1]][1]
  tree <- read.tree(paste0(cluster, ".astral.branch_length.clean.tree"))
  
  # correct names in tree
  tree <- correct.treenames(tree)
  
  # get subst of metadata
  sub_meta <- merged_meta[which(merged_meta$new_complete_name %in% tree$tip.label),]
  
  # reduce tree to include only one sample per city
  dropped_tree <- drop.multiple.locations(tree, sub_meta, cluster)
  
  if(dropped_tree$Nnode >1)
  {
    region_tree <-rename_taxa(dropped_tree, sub_meta, new_complete_name, Region)
    Region <- as.factor(c(region_tree$tip.label, rep(NA,dropped_tree$Nnode)))
    tree_col <- color_df[which(color_df$reg %in% levels(Region)),]
    tree_color <- as.character(tree_col$cols)
    
    
    # regular tree
    #p1 <- ggtree(region_tree)  + 
      #geom_tippoint(aes(color=Region), size=3) + 
      #ggtitle(cluster)
    #p <- ggtree(continent_tree, layout="equal_angle")  +  geom_tippoint(color=as.color(continent_tree$tip.label, 0.5)) + ggtitle(treeList[cluster])
    
    # circular tree
    #pdf(file = paste0("plots/", cluster, ".circular_tree.pdf"))
    #p2 <- ggtree(region_tree, layout="circular")  + 
      #geom_tippoint(aes(color=Region), size=3) + 
      #ggtitle(cluster) + 
      #scale_color_manual(values = as.character(tree_col$cols)) + 
      #theme(text = element_text(size = 15)) +
      #theme(legend.position = "none")
    #print(p2)
    #dev.off()
    
    #### unrooted tree 1
    #pdf(file = paste0("plots/equalangle/", cluster,".equalangle_tree.pdf"))
    #p3 <- ggtree(region_tree,layout="equal_angle")  + 
      #geom_tippoint(aes(color=Region), size=5) + 
      #ggtitle(cluster) + 
      #scale_color_manual(values = as.character(tree_col$cols)) +
      #theme(text = element_text(size = 20)) +
      #theme(legend.position = "none")
    #print(p3)
    #dev.off()
    
    #### unrooted tree 2
    pdf(file = paste0("plots/unroot/", cluster, ".unroot_tree.pdf"))
    p4 <- ggtree(region_tree, layout="daylight")  + 
      geom_tippoint(aes(color=Region), size=5) + 
      ggtitle(cluster) + 
      scale_color_manual(values = as.character(tree_col$cols)) + 
      theme(text = element_text(size = 20)) +
      theme(legend.position = "none")
    print(p4)
    dev.off()
    
    #grid.arrange(p1, p2, p3, p4, nrow=2)
    
    #plot_list[[i]] = p2
    
    #print(p1)
    #print(p2)
    #print(p3)
    #print(p4)
    print(cluster)
  }
  else
  {
    print(paste(cluster,"has only one node left"))
  }
}


#dev.off()


############## testing ##############
cluster_no <- c(1,2,3,4,5,6,7,8,9,11,12,13,14,15,16,17,18,20,21,22,23,24,25,28,29,30,31,32,34,37,38,41)
clusters <- c()

grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8,p9, p10, p11, p12,p13, p14, p15, p16,p17, p18, p19, p20,p21, p22, p23, p24,p25, p26, p27, p28,p29, p30, p31, p32, p33, nrow=6)


plot_list

plist <- list(p2, p4)
plist <- list(plist)
n <- length(plist)
nCol <- floor(sqrt(n))
do.call("grid.arrange", c(ggtree(plot_list), ncol=2))
multiplot(p1, p2, p3, p4, p5, p6, p7, p8,p9, p10, p11, p12,p13, p14, p15, p16,p17, p18, p19, p20,p21, p22, p23, p24,p25, p26, p27, p28,p29, p30, p31, p32, p33, ncol=2)


plist <- list(p2, p2)





