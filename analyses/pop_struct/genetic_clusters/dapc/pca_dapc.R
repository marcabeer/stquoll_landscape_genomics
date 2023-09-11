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
#Load in genetic data and metadata
##############################

#load in genind
stq_genind <- readRDS("stq_genind")

#longitude & latitude coordinates
latlongs <- read.table("latlongs", sep="\t", header=TRUE, stringsAsFactors=FALSE)

#other metadata
md <- read.table("metadata", sep="\t", stringsAsFactors=FALSE, header=TRUE)


##############################
#Carry out PCA
##############################

#carry out PCA on genetic data
stq_scaled <- scaleGen(stq_genind, NA.method="mean")
stq_pca <- dudi.pca(stq_scaled,  center=TRUE, scale=FALSE, scannf=FALSE, nf=6)
stq_pca_df <- data.frame(PC1=stq_pca$li[,1], PC2=stq_pca$li[,2], PC3=stq_pca$li[,3], long=latlongs$X, lat=latlongs$Y, year=md$year)

#percent variance explained by PCs
pc1_var <- stq_pca$eig[1]/sum(stq_pca$eig)
pc2_var <- stq_pca$eig[2]/sum(stq_pca$eig)

##############################
#Plot PCA results
##############################

#Plot PC1 vs PC2 and color by longitude and latitude
pca_long <- ggplot()+
  geom_point(data=stq_pca_df, mapping=aes(x=PC1, y=PC2, fill=long), shape=21, size=2)+
  xlim(c(min(stq_pca_df$PC1), max(stq_pca_df$PC1)))+
  ylim(c(min(stq_pca_df$PC2), max(stq_pca_df$PC2)))+
  scale_fill_gradient(low = "blue", high = "yellow")+
  labs(fill="Longitude", x="PC1 (6.96%)", y="PC2 (3.14%)")+
  theme_light()+
  theme(legend.position=c(0.15,0.2))

pca_lat<-ggplot()+
  geom_point(data=stq_pca_df, mapping=aes(x=PC1, y=PC2, fill=lat), shape=21, size=2)+
  xlim(c(min(stq_pca_df$PC1), max(stq_pca_df$PC1)))+
  ylim(c(min(stq_pca_df$PC2), max(stq_pca_df$PC2)))+
  scale_fill_gradient(low = "blue", high = "yellow")+
  labs(fill="Latitude", x="PC1 (6.96%)", y="PC2 (3.14%)")+
  theme_light()+
  theme(legend.position=c(0.11,0.15))+
  theme(legend.position=c(0.15,0.2))

pca_long
pca_lat

#Plot 

###
#plot with basemap
library(ggplot2)
library(ggspatial)
library(sf)
library(rgeos)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)

#get basemap for Tasmania
#https://r-spatial.org/r/2018/10/25/ggplot2-sf.html
aus <- ne_countries(scale = "large", country="Australia", returnclass = "sf")
tas<-st_crop(aus, xmin=144.5, xmax=149, ymin=-44, ymax=-40.5)

pc1_map<-ggplot(data=tas)+
  geom_sf(fill="grey92") +
  geom_point(data=stq_pca_df, aes(x=long, y=lat, fill=PC1), shape=21, size=3)+
  scale_fill_gradient(low = "blue", high = "yellow")+
  labs(fill="PC1")+
  xlab("Longitude")+
  ylab("Latitude")+
  theme_light()+
  theme(legend.position=c(0.11,0.15))

pc2_map<-ggplot(data=tas)+
  geom_sf(fill="grey92") +
  geom_point(data=stq_pca_df, aes(x=long, y=lat, fill=PC2), shape=21, size=3)+
  scale_fill_gradient(low = "blue", high = "yellow")+
  labs(fill="PC2")+
  xlab("Longitude")+
  ylab("Latitude")+
  theme_light()+
  theme(legend.position=c(0.11,0.15))


##############################
#Discriminant Analysis of Principal Components (DAPC)
##############################

###
#determine number of clusters with the most support (lowest BIC)
stq_clusters <- find.clusters(stq_genind, n.pca=345)

#BIC plot
dapc_bic_df<-data.frame(K=seq(1:20), BIC=stq_clusters$Kstat[1:20])

ggplot(dapc_bic_df, aes(x=K, y=BIC)) + 
  geom_point(size=1.25)+
  geom_line(size=0.6)+
  xlim(0,20)+
  scale_x_continuous(breaks=seq(0, 20, 2))+
  xlab("Number of clusters")+
  ylab("BIC")+
  theme_bw()+
  theme(legend.position="none")


###
#carry out an initial dapc to determine number of discriminant axes to retain
##should show biggest F-statistic for one DA, and small statistic for two DA; we'll retain two DA for now
stq_dapc_test <- dapc(stq_genind, stq_clusters$grp, n.pca=345, n.da=2)

#use a-score optimization to determine number of PCs to retain
##higher a-score is better
optim.a.score(stq_dapc_test, n.pca=50)

#re-run DAPC using the optimal number of PC axes
stq_dapc_ascore <- dapc(stq_genind, stq_clusters$grp, n.pca=1, n.da=2)
scatter(stq_dapc_ascore)

#create dapc dataframe to store results
dapc_df<-data.frame(id=rownames(stq_dapc_ascore$ind.coord), longitude=latlongs[,1], latitude=latlongs[,2], radius=0.075, dapc=stq_dapc_ascore$posterior)
#write.table(dapc_df, "dapc_df", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

#if not re-running the analyses, you can load our results below:
dapc_df <- read.table("dapc_df", header=TRUE, sep="\t")


###
#Assign samples to genetic clusters based on majority contribution

#order samples by longitude (which is the geographic axis on which genetic variation is largely distributed)
dapc_df_ordered<-dapc_df[order(dapc_df$longitude),]
dapc_df_ordered$integer<-seq(1:nrow(dapc_df_ordered))

#label three genetic clusters as western Tasmania (WT), central Tas (CT), and eastern Tas (ET)
for (i in 1:nrow(dapc_df_ordered)){
  cluster_index<-which(dapc_df_ordered[i,5:7]==max(dapc_df_ordered[i,5:7]))
  if (cluster_index==1){
    dapc_df_ordered$cluster[i]<-"WT"
  }else{
    if (cluster_index==2){
      dapc_df_ordered$cluster[i]<-"CT"
    }else{
      dapc_df_ordered$cluster[i]<-"ET"
    }
  }
}

#save results to tables
#write.table(dapc_df_ordered, "dapc_df_ordered", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

#if not re-running the analyses, you can load our results below:
dapc_df_ordered <- read.table("dapc_df_ordered", header=TRUE, sep="\t")


##############################
#Plot individual DAPC ancestry coefficients
##############################

#convert dataframe to long format
dapc_df_ordered_long<-reshape2::melt(data=dapc_df_ordered, id.vars=c(1:4,8), measure.vars=5:7)

###
#ancestry stackplot
ggplot(dapc_df_ordered_long, aes(fill=variable, y=value, x=as.factor(integer)))+
  geom_bar( inherit.aes=TRUE, stat="identity")+
  scale_fill_manual(values=c("#2271B2", "#E69F00", "#359B73"))+
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

world_hires<-ggplot2::map_data("worldHires")
tas_2<-ggplot2::map_data("worldHires", xlim=c(144.5,149), ylim=c(-44,-40.5))

#make Tasmania basemap
p<-ggplot(tas_2, aes(x=long, y=lat, group=group))+
  geom_polygon(fill="grey92", color="black")+
  coord_quickmap()+
  xlab("Longitude")+
  ylab("Latitude")+
  theme_light()

#create pie charts for each  value of K
dapc_pie.list <- dapc_df %>% 
  tidyr::gather(type, value, -longitude, -latitude, -radius, -id) %>%
  tidyr::nest(type, value) %>%
  
  # make a pie chart from each row, & convert to grob
  mutate(pie.grob = purrr::map(data,
                               function(d) ggplotGrob(ggplot(d, 
                                                             aes(x = 1, y = value, fill = type)) +
                                                        geom_col(color = "black",
                                                                 show.legend = FALSE) +
                                                        coord_polar(theta = "y") +
                                                        scale_fill_manual(values=c("#2271B2", "#E69F00", "#359B73"))+
                                                        theme_void()))) %>%
  
  # convert each grob to an annotation_custom layer. I've also adjusted the radius
  # value to a reasonable size (based on my screen resolutions).
  rowwise() %>%
  mutate(radius = 0.075) %>%
  mutate(subgrob = list(annotation_custom(grob = pie.grob,
                                          xmin = longitude - radius, xmax = longitude + radius,
                                          ymin = latitude - radius, ymax = latitude + radius)))


dapc_map <- p + 
  geom_tile(data = dapc_df %>% tidyr::gather(type, value, -longitude, -latitude, -radius, -id),
            aes(x = longitude, y = latitude, fill = type), 
            color = "black", width = 0.01, height = 0.01, 
            inherit.aes = FALSE) +
  dapc_pie.list$subgrob+
  theme(legend.position="none")+
  labs(fill="Cluster")

dapc_map