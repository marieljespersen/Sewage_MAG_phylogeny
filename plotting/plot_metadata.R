library(ggplot2)
library(lubridate, warn.conflicts = FALSE)

setwd("~/HCP_Anywhere/GS2/GS2_metadata")
merged_meta <- read.delim("/Users/maloj/HCP_Anywhere/GS2/GS2_metadata/merged_metadata.txt", stringsAsFactors = F)

#metadata_with_nc <- read.delim("gs_2_metadata.csv", header = TRUE, sep = ";")
#metadata <- metadata_with_nc[which(metadata_with_nc$country!="Negative control" & metadata_with_nc$country!=""),]

sample_time <- ymd(merged_meta$collection_date)
df <- as.data.frame(sample_time)
df$country <- merged_meta$country
df$Region <- merged_meta$Region

#qplot(sample_time, country, data=df, colour = Region)

# Ordering the dataframe according to continent for nice visualisation
df$Region<- factor(df$Region,levels=unique(df$Region))
df = df[order(df$Region), ]
df$country = factor(df$country, levels = rev(unique(df$country)))

# Plotting
p <- ggplot(df, aes(x=sample_time, y=country, color=Region)) + 
  geom_point() +
  theme_bw() +
  labs(y = "Country", x = "Sample time") +
  theme(legend.title=element_text(size=14), 
        legend.text=element_text(size=14),
        axis.title = element_text(size=14),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12)) + 
  labs(y = "Country", x = "Sample time") +
  theme(legend.position = "top")
ggsave(plot=p,width=5,height=15,filename="/Users/maloj/HCP_Anywhere/GS2/Results/metadata_plot.pdf", useDingbats=FALSE)
ggsave(plot=p,width=15,height=15,filename="/Users/maloj/HCP_Anywhere/GS2/Results/metadata_plot_legend.pdf", useDingbats=FALSE)


