library("ggmap")
library(maptools)
library("maps")

merged_meta <- read.delim("/Users/maloj/HCP_Anywhere/GS2/GS2_metadata/merged_metadata.txt", stringsAsFactors = F)
merged_meta$latitude_inferred <- str_replace(merged_meta$latitude_inferred, ",", ".")
merged_meta$longtitude_inferred <- str_replace(merged_meta$longtitude_inferred, ",", ".")

visit.y <- as.numeric(merged_meta$latitude_inferred)
visit.x <- as.numeric(merged_meta$longtitude_inferred)

#Using GGPLOT, plot the Base World Map
mp <- NULL
mapWorld <- borders("world", colour="gray50", fill="gray50") # create a layer of borders
mp <- ggplot() +   mapWorld

#Now Layer the cities on top
p <- mp + geom_point(data=merged_meta, 
             aes(x=as.numeric(longtitude_inferred), 
                 y=as.numeric(latitude_inferred),
                 color=Region), 
             size=2) +
  theme_bw()+
  theme(axis.text.x = element_text(color = "grey20", size = 8, angle = 0, hjust = 1, vjust = 0, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 8, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain"))
p + theme(legend.position="top")
