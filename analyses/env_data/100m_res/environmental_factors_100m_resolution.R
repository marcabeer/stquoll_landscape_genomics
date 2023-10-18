setwd("C:/Users/icepi/Documents/GitHub/stquoll_landscape_genomics/analyses/env_data/100m_res")
library(raster)
library(geobuffer)
library(geodist)
library(sp)
library(corrplot)


##############################
#Import devil density and DFTD data
##############################

#import DFTD arrival raster
new_prj<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

dftd <- raster("DFTD_spread_raster_fig3B.tif")
dftd_reproj <- projectRaster(dftd, crs=new_prj, method="ngb")

#import devil density rasterstack
devil_rasterstack<-stack("predictionStack_devils_1985to2035.tif")
devil_reproj<-projectRaster(devil_rasterstack, crs=new_prj, method="ngb")


##############################
#Import binary raster for each TASVEG landcover type
##############################
#binary rasters made using QGIS
#https://gis.stackexchange.com/questions/121532/how-to-reclass-a-raster-with-reclassify-grid-values-in-qgis

#modified
tasveg_100m_mod<-raster(x="tasveg_100m_mod.tif")
#dry eucalypt forest
tasveg_100m_def<-raster(x="tasveg_100m_def.tif")
#highland treeless
tasveg_100m_htl<-raster(x="tasveg_100m_htl.tif")
#moorland
tasveg_100m_mol<-raster(x="tasveg_100m_mol.tif")
#native grassland
tasveg_100m_ngl<-raster(x="tasveg_100m_ngl.tif")
#non-eucalypt forest
tasveg_100m_nef<-raster(x="tasveg_100m_nef.tif")
#rainforest
tasveg_100m_raf<-raster(x="tasveg_100m_raf.tif")
#saltmarsh wetland
tasveg_100m_swl<-raster(x="tasveg_100m_swl.tif")
#scrub
tasveg_100m_scr<-raster(x="tasveg_100m_scr.tif")
#wet eucalypt forest
tasveg_100m_wef<-raster(x="tasveg_100m_wef.tif")
#other
tasveg_100m_oth<-raster(x="tasveg_100m_oth.tif")


##############################
#Import WorldClim data
##############################
#already has WGS84 CRS, so no need to reproject

#bio1 = MAT
bio1 <- raster("wc2.1_30s_bio_1.tif")
#bio2 = mean diurnal range
bio2 <- raster("wc2.1_30s_bio_2.tif")
#bio3 = isothermality
bio3 <- raster("wc2.1_30s_bio_3.tif")
#bio4 = temp seasonality
bio4 <- raster("wc2.1_30s_bio_4.tif")
#bio12 = annual precip
bio12 <- raster("wc2.1_30s_bio_12.tif")
#bio15 = precip seasonality
bio15 <- raster("wc2.1_30s_bio_15.tif")
#elevation
elev <- raster("wc2.1_30s_elev.tif")

###
#crop WorldClim data to make manipulation easier
s_extent<-raster::extent(143.4833,148.7417,-43.8,-39.18333 )

bio1_scrop<-raster::crop(bio1, s_extent)
bio2_scrop<-raster::crop(bio2, s_extent)
bio3_scrop<-raster::crop(bio3, s_extent)
bio4_scrop<-raster::crop(bio4, s_extent)
bio12_scrop<-raster::crop(bio12, s_extent)
bio15_scrop<-raster::crop(bio15, s_extent)
elev_scrop<-raster::crop(elev, s_extent)
bioclim_cropped<-stack(bio1_scrop, bio2_scrop, bio3_scrop, bio4_scrop, bio12_scrop, bio15_scrop, elev_scrop)

###
#if not re-running previous WorldClim-related code, load in the cropped rasters
bioclim_cropped <- stack("bioclim_cropped.tif")


##############################
#Resample data to enable stacking and subsequent data extraction
##############################

dftd_resamp<-resample(x=dftd_reproj, y=tasveg_100m_mod, method="bilinear")
devil_resamp<-resample(x=devil_reproj, y=tasveg_100m_mod, method="bilinear")
bioclim_resamp<-resample(x=bioclim_cropped, y=tasveg_100m_mod, method="bilinear")

#stack all resampled rasters together
all_vars<-stack(bioclim_resamp,
                   tasveg_100m_mod,
                   tasveg_100m_def,
                   tasveg_100m_htl,
                   tasveg_100m_mol,
                   tasveg_100m_ngl,
                   tasveg_100m_nef,
                   tasveg_100m_raf,
                   tasveg_100m_swl,
                   tasveg_100m_scr,
                   tasveg_100m_wef,
                   tasveg_100m_oth,
                   dftd_resamp,
                   devil_resamp)

writeRaster(x=all_vars, filename="all_vars_100m", format="GTiff")


##############################
#Data extraction at sampling localities
##############################

#load in rasterstack
all_vars<-stack("all_vars_100m.tif")

#load in metadata
latlongs<-read.csv("latlongs", sep="\t", header=TRUE)

#calculate buffer radius based on published home range size estimates
#Andersen et al. (2020): 626 hectare = 6.26 sqkm
rad_and<-sqrt(6.26/pi)

#buffer latlongs
latlong_buffer<-geobuffer::geobuffer_pts(xy=latlongs, dist_m=(rad_and*1000), step_dg=10)

#extract environmental data within buffers
#env_mean<-raster::extract(x=all_vars, y=latlong_buffer, fun="mean", na.rm=TRUE, df=TRUE, layer=1, nl=70, sp=FALSE)
write.table(x=env_mean, file="env_data_extracted", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

#read in environmental data
env_mean<-read.table("env_data_extracted", header=TRUE, sep="\t")

#some localities may have NA values due to rasterization of data
##obtain latlongs of points with NA values
na_index<-which(is.na(env_mean[,21]))
na_latlongs<-latlongs[na_index,]

#convert one of the raster layers to a df
r_df<-raster::as.data.frame(x=all_vars[[1]], xy=TRUE, na.rm=TRUE)

#calculate distance between each of the NA sample coordinates and coordinates that do have raster values
na_dist<-geodist::geodist(x=na_latlongs, y=r_df[,1:2], measure="geodesic")

#replace NA sample coordinates with coordinates of nearest raster cell with environmental values
##very minorly changes the location of 3 samples
na_dist_min_index<-c()
for (i in 1:nrow(na_latlongs)){
  na_dist_vec<-as.vector(na_dist[i,])
  na_dist_min_index[i]<-which(na_dist_vec==min(na_dist_vec))
}

na_new_pts<-r_df[na_dist_min_index,1:2]

#buffer the new points and extract data at them
na_new_buffer<-geobuffer::geobuffer_pts(xy=na_new_pts, dist_m=(rad_and*1000), step_dg=10)
na_env_mean<-raster::extract(x=all_vars, y=na_new_buffer, fun="mean", na.rm=TRUE, df=TRUE, layer=1, nl=70, sp=FALSE)

#input new non-na values into original env_mean df
for (i in 1:length(na_index)){
  env_mean[na_index[i], 2:ncol(na_env_mean)]<-na_env_mean[i,-1]
}

write.table(x=env_mean, file="env_data_nafilled_extracted", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)


##############################
#Further process environmental data
##############################

###
#load in environmental and spatiotemporal data
latlongs<-read.table("latlongs", sep="\t", header=TRUE)
metadata<-read.table("metadata", sep="\t", header=TRUE)
env_mean<-read.table("env_data_nafilled_extracted", header=TRUE, sep="\t")

###
#calculate years/generations diseased
gen_dftd<-metadata$year-env_mean$DFTD_spread_raster_fig3B

#convert negative values to zeroes
gen_dftd[which(gen_dftd<0)]<-0

#replace dftd arrival year column with gens diseased variable
env_mean_v2<-env_mean
env_mean_v2$DFTD_spread_raster_fig3B<-gen_dftd


###
#isolate devil densities for years 1985 - 2035
devil_densities<-env_mean_v2[,21:ncol(env_mean_v2)]
colnames(devil_densities)<-1985:2035

#determine what year needs to be extracted for 5, 10, and 15 generations pre-collection year
gen5_year<-metadata$year-5
gen10_year<-metadata$year-10
gen15_year<-metadata$year-15

#isolate devil densities corresponding to the 5, 10, 15-generation lags
gen5_devil<-c(rep(NA, nrow(metadata)))
gen10_devil<-c(rep(NA, nrow(metadata)))
gen15_devil<-c(rep(NA, nrow(metadata)))

for (i in 1:nrow(metadata)){
  gen5_devil[i]<-devil_densities[i, which(colnames(devil_densities)==as.character(gen5_year[i]))]
  gen10_devil[i]<-devil_densities[i, which(colnames(devil_densities)==as.character(gen10_year[i]))]
  gen15_devil[i]<-devil_densities[i, which(colnames(devil_densities)==as.character(gen15_year[i]))]
}


###
#create new df of environmental predictors
env<-data.frame(env_mean_v2[,2:19], gen_dftd=env_mean_v2[,20], gen5_devil, gen10_devil, gen15_devil)

#check that variation exists for each environmental variable
env_var<-apply(env, 2, var)

#remove variable(s) that have zero variance
env<-env[,-which(env_var==0)]

#Bio3 and Bio4 are multiplied by 100 in the original WorldClim dataset
env$wc2.1_30s_bio_3<-env$wc2.1_30s_bio_3/100
env$wc2.1_30s_bio_4<-env$wc2.1_30s_bio_4/100

#output data
#write.table(env, "env_all", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

###
#check inter-environmental correlations
#practical guide to landscape genomics cites paper removing factors based on |r|>=0.8

#calculate correlation coefficient
env_cor_abs<-abs(cor(env))
env_cor_abs[upper.tri(x=abs(cor(env)), diag=TRUE)]<-NA

#the lagged devil density variables are highly collinear (|r| > 0.90)
##keep the longest-lagged devil density because shorter lags are less likely to be reflected in genetic data
env<-subset(env, select=-c(gen5_devil, gen10_devil))

#bio1 and bio2 are highly correlated with other variables
env<-subset(env, select=-c(wc2.1_30s_bio_1, wc2.1_30s_bio_2))

#save final set of moderately correlated environmental variables
write.table(x=env, "env_r80", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t ")


##############################
#Create correlation plot
##############################

#create new df
env_new<-data.frame(env_mean_v2[,2:19], gen_dftd=env_mean_v2[,20], gen5_devil, gen10_devil, gen15_devil)
env_new<-env_new[,-which(env_var==0)]

#load in decimal year
dates <- read.table("stq_dates_yeardec", sep="\t", header=TRUE)
env_new$year <- dates$yeardec
colnames(env_new)<-c("bio1", "bio2", "bio3", "bio4", "bio12", "bio15", "elev", "tasveg_mod", "tasveg_def", "tasveg_mol", "tasveg_ngl", "tasveg_nef", "tasveg_raf", "tasveg_swl", "tasveg_scr", "tasveg_wef", "tasveg_oth", "gen_dftd", "gen5_devil", "gen10_devil", "gen15_devil", "year")

#create correlation plot
par(mfrow=c(1,1))
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(cor(env_new), method="color", col=col(200), type="upper", addgrid.col="gray70", addCoef.col="gray10", tl.cex=1, tl.col="black", number.cex=0.7, diag=TRUE)
