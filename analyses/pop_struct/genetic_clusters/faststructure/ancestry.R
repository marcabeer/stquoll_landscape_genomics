#plotting
library(ggpubr)
library(ggplot2)
library(cowplot)
library(egg)
library(grid)
library(scatterpie)

#spatial analysis
library(ggspatial)
library(sf)
library(rgeos)
library(rnaturalearth)
library(rnaturalearthdata)
#High resolution naturalearth data requires installation from GitHub
library(devtools)
#devtools::install_github("ropensci/rnaturalearthhires")
library(rnaturalearthhires)
library(mapdata)
library(tidyr)
library(dplyr)
library(reshape2)

#genetic analysis
library(vcfR)
library(adegenet)

##############################
#Load in  metadata
##############################
latlongs <- read.table("latlongs", sep="\t", header=TRUE, stringsAsFactors=FALSE)
md <- read.table("metadata", sep="\t", stringsAsFactors=FALSE, header=TRUE)


##############################
#Load in coancestry matrices for select values of K (i.e., Q matrices)
##############################

#load in Q matrix
k3_q <- read.table("test.r1.3.meanQ", header=FALSE, sep=" ")
k4_q <- read.table("test.r1.4.meanQ", header=FALSE, sep=" ")
k5_q <- read.table("test.r1.5.meanQ", header=FALSE, sep=" ")
k6_q <- read.table("test.r1.6.meanQ", header=FALSE, sep=" ")
k7_q <- read.table("test.r1.7.meanQ", header=FALSE, sep=" ")
k8_q <- read.table("test.r1.8.meanQ", header=FALSE, sep=" ")
k9_q <- read.table("test.r1.9.meanQ", header=FALSE, sep=" ")

#create Q matrix metadata dataframes
k3_df <- data.frame(id=md$X, longitude=latlongs[,1], latitude=latlongs[,2], radius=0.075, fs=k3_q)
k4_df <- data.frame(id=md$X, longitude=latlongs[,1], latitude=latlongs[,2], radius=0.075, fs=k4_q)
k5_df <- data.frame(id=md$X, longitude=latlongs[,1], latitude=latlongs[,2], radius=0.075, fs=k5_q)
k6_df <- data.frame(id=md$X, longitude=latlongs[,1], latitude=latlongs[,2], radius=0.075, fs=k6_q)
k7_df <- data.frame(id=md$X, longitude=latlongs[,1], latitude=latlongs[,2], radius=0.075, fs=k7_q)
k8_df <- data.frame(id=md$X, longitude=latlongs[,1], latitude=latlongs[,2], radius=0.075, fs=k8_q)
k9_df <- data.frame(id=md$X, longitude=latlongs[,1], latitude=latlongs[,2], radius=0.075, fs=k9_q)

###
#Assign samples to genetic clusters based on majority contribution
##creates a table analogous to the one output in the DAPC script

#order individuals by longitude
k3_df_ordered <- k3_df[order(k3_df$longitude),]
k3_df_ordered$integer <- seq(1:nrow(k3_df_ordered))

#label three genetic clusters as western Tasmania (WT), central Tas (CT), and eastern Tas (ET)
for (i in 1:nrow(k3_df_ordered)){
  cluster_index<-which(k3_df_ordered[i,5:7]==max(k3_df_ordered[i,5:7]))
  if (cluster_index==3){
    k3_df_ordered$cluster[i]<-"WT"
  }else{
    if (cluster_index==2){
      k3_df_ordered$cluster[i]<-"CT"
    }else{
      k3_df_ordered$cluster[i]<-"ET"
    }
  }
}

#save output in table
#write.table(k3_df_ordered, "k3_df_ordered", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

#if not re-running the analyses, you can load our results below:
k3_df_orderd <- read.table("k3_df_ordered", header=TRUE, sep="\t")


##############################
#Plot individual FastStructure ancestry coefficients for K=3
##############################

#convert to long format for ggplot
k3_df_ordered_long <- reshape2::melt(data=k3_df_ordered, id.vars=c(1:4,8), measure.vars=5:7)

#Create ancestry stackplot
ggplot(k3_df_ordered_long, aes(fill=variable, y=value, x=as.factor(integer)))+
  geom_bar( inherit.aes=TRUE, stat="identity")+
  scale_fill_manual(values=c("#359B73", "#E69F00", "#2271B2"))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.text=element_blank(),
        legend.title=element_blank(),
        legend.position="none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())


###
#plot ancestries on map
world_hires <- ggplot2::map_data("worldHires")
tas_2 <- ggplot2::map_data("worldHires", xlim=c(144.5,149), ylim=c(-44,-40.5))

#make Tasmania basemap
p <- ggplot(tas_2, aes(x=long, y=lat, group=group))+
  geom_polygon(fill="grey92", color="black")+
  coord_quickmap()+
  xlab("Longitude")+
  ylab("Latitude")+
  theme_light()

#create pie charts for each  value of K
k3_pie.list <- k3_df %>% 
  tidyr::gather(type, value, -longitude, -latitude, -radius, -id) %>%
  tidyr::nest(type, value) %>%
  
  # make a pie chart from each row, & convert to grob
  mutate(pie.grob = purrr::map(data,
                               function(d) ggplotGrob(ggplot(d, 
                                                             aes(x = 1, y = value, fill = type)) +
                                                        geom_col(color = "black",
                                                                 show.legend = FALSE) +
                                                        coord_polar(theta = "y") +
                                                        scale_fill_manual(values=c("#359B73", "#E69F00", "#2271B2"))+
                                                        theme_void()))) %>%
  
  # convert each grob to an annotation_custom layer. I've also adjusted the radius
  # value to a reasonable size (based on my screen resolutions).
  rowwise() %>%
  mutate(radius = 0.075) %>%
  mutate(subgrob = list(annotation_custom(grob = pie.grob,
                                          xmin = longitude - radius, xmax = longitude + radius,
                                          ymin = latitude - radius, ymax = latitude + radius)))

#add individual ancestry pies to basemap
k3_map <- p + 
  geom_tile(data = k3_df %>% tidyr::gather(type, value, -longitude, -latitude, -radius, -id),
            aes(x = longitude, y = latitude, fill = type), 
            color = "black", width = 0.01, height = 0.01, 
            inherit.aes = FALSE) +
  k3_pie.list$subgrob+
  theme(legend.position="none")+
  labs(fill="Cluster")

k3_map


##############################
#Plots for other values of K
##############################

k4_df_ordered <- k4_df[order(k4_df$longitude),]
k4_df_ordered$integer <- seq(1:nrow(k4_df_ordered))
k4_df_ordered_long <- reshape2::melt(data=k4_df_ordered, id.vars=c(1:4,9), measure.vars=5:8)

k5_df_ordered <- k5_df[order(k5_df$longitude),]
k5_df_ordered$integer <- seq(1:nrow(k5_df_ordered))
k5_df_ordered_long <- reshape2::melt(data=k5_df_ordered, id.vars=c(1:4,10), measure.vars=5:9)

k6_df_ordered <- k6_df[order(k6_df$longitude),]
k6_df_ordered$integer <- seq(1:nrow(k6_df_ordered))
k6_df_ordered_long <- reshape2::melt(data=k6_df_ordered, id.vars=c(1:4,11), measure.vars=5:10)

k7_df_ordered <- k7_df[order(k7_df$longitude),]
k7_df_ordered$integer <- seq(1:nrow(k7_df_ordered))
k7_df_ordered_long <- reshape2::melt(data=k7_df_ordered, id.vars=c(1:4,12), measure.vars=5:11)

k8_df_ordered <- k8_df[order(k8_df$longitude),]
k8_df_ordered$integer <- seq(1:nrow(k8_df_ordered))
k8_df_ordered_long <- reshape2::melt(data=k8_df_ordered, id.vars=c(1:4,13), measure.vars=5:12)

k9_df_ordered <- k9_df[order(k9_df$longitude),]
k9_df_ordered$integer <- seq(1:nrow(k9_df_ordered))
k9_df_ordered_long <- reshape2::melt(data=k9_df_ordered, id.vars=c(1:4,14), measure.vars=5:13)


#create individual ancestry pies for mapping
k4_pie.list <- k4_df %>% 
  tidyr::gather(type, value, -longitude, -latitude, -radius, -id) %>%
  tidyr::nest(type, value) %>%
  
  # make a pie chart from each row, & convert to grob
  mutate(pie.grob = purrr::map(data,
                               function(d) ggplotGrob(ggplot(d, 
                                                             aes(x = 1, y = value, fill = type)) +
                                                        geom_col(color = "black",
                                                                 show.legend = FALSE) +
                                                        coord_polar(theta = "y") +
                                                        scale_fill_manual(values=c("#2271B2", "#E69F00", "#359B73", "#D55E00"))+
                                                        theme_void()))) %>%
  
  # convert each grob to an annotation_custom layer. I've also adjusted the radius
  # value to a reasonable size (based on my screen resolutions).
  rowwise() %>%
  mutate(radius = 0.075) %>%
  mutate(subgrob = list(annotation_custom(grob = pie.grob,
                                          xmin = longitude - radius, xmax = longitude + radius,
                                          ymin = latitude - radius, ymax = latitude + radius)))

k5_pie.list <- k5_df %>% 
  tidyr::gather(type, value, -longitude, -latitude, -radius, -id) %>%
  tidyr::nest(type, value) %>%
  
  # make a pie chart from each row, & convert to grob
  mutate(pie.grob = purrr::map(data,
                               function(d) ggplotGrob(ggplot(d, 
                                                             aes(x = 1, y = value, fill = type)) +
                                                        geom_col(color = "black",
                                                                 show.legend = FALSE) +
                                                        coord_polar(theta = "y") +
                                                        scale_fill_manual(values=c("#CC79A7", "#E69F00", "#D55E00", "#2271B2", "#359B73"))+
                                                        theme_void()))) %>%
  
  # convert each grob to an annotation_custom layer. I've also adjusted the radius
  # value to a reasonable size (based on my screen resolutions).
  rowwise() %>%
  mutate(radius = 0.075) %>%
  mutate(subgrob = list(annotation_custom(grob = pie.grob,
                                          xmin = longitude - radius, xmax = longitude + radius,
                                          ymin = latitude - radius, ymax = latitude + radius)))


#make maps

k4_map <- p + 
  geom_tile(data = k4_df %>% tidyr::gather(type, value, -longitude, -latitude, -radius, -id),
            aes(x = longitude, y = latitude, fill = type), 
            color = "black", width = 0.01, height = 0.01, 
            inherit.aes = FALSE) +
  k4_pie.list$subgrob+
  theme(legend.position="none")+
  labs(fill="Cluster")
k4_map

k5_map<- p + 
  geom_tile(data = k5_df %>% tidyr::gather(type, value, -longitude, -latitude, -radius, -id),
            aes(x = longitude, y = latitude, fill = type), 
            color = "black", width = 0.01, height = 0.01, 
            inherit.aes = FALSE) +
  k5_pie.list$subgrob+
  theme(legend.position="none")+
  labs(fill="Cluster")
k5_map
