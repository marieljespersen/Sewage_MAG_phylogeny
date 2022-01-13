library(dplyr)
library(ggpubr)
library(tidyr)
library(ggExtra)

setwd("/Users/maloj/HCP_Anywhere/GS2/Results/")
tax_R2 <- read.delim("cluster_tax.tsv", sep="\t", header=T)

levels(tax_R2$genus)

# get summary statistics by phylum
group_by(tax_R2, genus) %>%
  summarise(
    count = n(),
    mean = mean(R2, na.rm = TRUE),
    sd = sd(R2, na.rm = TRUE),
    median = median(R2, na.rm = TRUE),
    IQR = IQR(R2, na.rm = TRUE)
  )

#### Box plots ####
# Plot R2 by phylum and color by phylum
xlabs <- paste(levels(droplevels(tax_R2$phylum[which(!is.na(tax_R2$R2))])),"\n(N=",table(droplevels(tax_R2$phylum[which(!is.na(tax_R2$R2))])),")",sep="")
p  <- ggboxplot(tax_R2[which(!is.na(tax_R2$R2) & !is.na(tax_R2$phylum)),], x = "phylum", y = "R2", color = "phylum") + 
  scale_x_discrete(labels=xlabs) +
  theme_bw() +
  theme(axis.text.x=element_text(color = "black", size=11, angle=30, vjust=.8, hjust=0.8)) + 
  ylab(expression(R^2 ~ "geography")) + xlab("") +
  removeGrid() +
  theme(legend.position="none")
ggsave(plot=p,width=5,height=5,dpi=200,filename="/Users/maloj/HCP_Anywhere/GS2/Results/R2_results/R2_phylum.pdf", useDingbats=FALSE)


#### Testing #####
# test for significant differences between the different phylum
kruskal.test(R2 ~ phylum, data = tax_R2)


