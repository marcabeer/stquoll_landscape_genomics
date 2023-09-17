library(raster)
library(geobuffer)
library(geodist)


##############################
#Import devil density data
##############################

#import DFTD arrival raster
new_prj<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

dftd <- raster("DFTD_spread_raster_fig3B.tif")
dftd_reproj <- projectRaster(dftd, crs=new_prj, method="ngb")

#import devil density rasterstack
devil_rasterstack<-stack("predictionStack_devils_1985to2035.tif")
plot(devil_rasterstack)
devil_reproj<-projectRaster(devil_rasterstack, crs=new_prj, method="ngb")
plot(devil_reproj)


##############################
#Import TASVEG raster (reclassified; landcover data)
##############################

tasveg<-raster("tasveg_1km_final.tif")


##############################
#Import WorldClim data
##############################

#bio1 = MAT
bio1<-raster("wc2.1_30s_bio_1.tif")
#bio2 = mean diurnal range
bio2<-raster("wc2.1_30s_bio_2.tif")
#bio3 = isothermality
bio3<-raster("wc2.1_30s_bio_3.tif")
#bio4 = temp seasonality
bio4<-raster("wc2.1_30s_bio_4.tif")
#bio12 = annual precip
bio12<-raster("wc2.1_30s_bio_12.tif")
#bio15 = precip seasonality
bio15<-raster("wc2.1_30s_bio_15.tif")
#elevation
elev<-raster("wc2.1_30s_elev.tif")


##############################
#Crop WorldClim data
##############################

bio1_scrop<-raster::crop(bio1, devil_reproj[[1]])
bio2_scrop<-raster::crop(bio2, devil_reproj[[1]])
bio3_scrop<-raster::crop(bio3, devil_reproj[[1]])
bio4_scrop<-raster::crop(bio4, devil_reproj[[1]])
bio12_scrop<-raster::crop(bio12, devil_reproj[[1]])
bio15_scrop<-raster::crop(bio15, devil_reproj[[1]])
elev_scrop<-raster::crop(elev, devil_reproj[[1]])

##############################
#Further process devil density and DFTD rasters
##############################

#resampling necessary to get same extents
dftd_resamp<-resample(x=dftd_reproj, y=bio1_scrop, method="ngb")
devil_resamp<-resample(x=devil_reproj, y=bio1_scrop, method="ngb")

#replace missing data in devil density rasters with zeroes
devil_resamp[is.na(devil_resamp[])] <- 0 
spplot(devil_resamp[[1]])
devil_resamp_remove_nas<-devil_resamp

#remove non-overlapping NAs between bioclimatic rasters and density rasters
devil_resamp_remove_nas[is.na(bio1_scrop[])]<-NA
spplot(devil_resamp_remove_nas[[5]])
writeRaster(x=devil_resamp_remove_nas, filename="devil_resampled_processed", format="GTiff")

#write dftd arrival raster
##gets processed in QGIS to interpolate DFTD arrival raster
writeRaster(x=dftd_resamp, filename="dftd_resampled", format="GTiff")

#load in dftd arrival raster with NAs having been interpolated in QGIS
dftd_resamp_nafilled<-raster("dftd_resampled_fillnodata_10pix.tif")

#remove non-overlapping NAs between bioclimatic rasters and dftd arrival rasters
dftd_resamp_nafilled[is.na(bio1_scrop[])]<-NA
spplot(dftd_resamp_nafilled)
writeRaster(x=dftd_resamp_nafilled, filename="dftd_resampled_processed", format="GTiff")

#read in resampled dftd and devil rasters
dftd_resamp_processed<-raster("dftd_resampled_processed.tif")
devil_resamp_processed<-stack("devil_resampled_processed.tif")


##############################
#Further process devil density and DFTD rasters
##############################

roads<-raster("roads_1km_v2.tif")
rivers<-raster("rivers_1km.tif")


##############################
#stack 1km rasters and export
##############################

all_vars_1km<-stack(bio1_scrop,
                    bio2_scrop,
                    bio3_scrop,
                    bio4_scrop,
                    bio12_scrop,
                    bio15_scrop,
                    elev_scrop,
                    tasveg,
                    roads,
                    rivers,
                    dftd_resamp_processed,
                    devil_resamp_processed)

#writeRaster(x=all_vars_1km, filename="all_vars_1km_v2", format="GTiff")


##############################
#Trim rasters to extent needed for ResistanceGA analysis
##############################

#load in 1km-res raster
all_vars_1km_res_raw<-stack("all_vars_1km_v2.tif")

#divide bio3 and bio3 by 100
##the raw WorldClim values are multipled by 100
##this reverts them to their true values
all_vars_1km_res_raw$all_vars_1km_v2.3<-all_vars_1km_res_raw$all_vars_1km_v2.3/100
all_vars_1km_res_raw$all_vars_1km_v2.4<-all_vars_1km_res_raw$all_vars_1km_v2.4/100

#average devil density rasters for given time intervals
##earliest devil sample included is in 2006
##average 2001-2007 for 5 gen
##average 1996-2002 for 10 gen
##average 1991-1997 for 15 gen
##average 1986-1992 for 20 gen
all_vars_devil<-all_vars_1km_res_raw[[12:62]]
names(all_vars_devil)<-paste("devil", 1985:2035, sep="_")
names(all_vars_devil)

devil_5gen<-mean(all_vars_devil[[17:23]])
devil_10gen<-mean(all_vars_devil[[12:18]])
devil_15gen<-mean(all_vars_devil[[7:13]])
devil_20gen<-mean(all_vars_devil[[2:8]])

#re-stack variables with the averaged devil densities; omit years diseased
all_vars_trimmed<-stack(all_vars_1km_res_raw[[1:10]], devil_5gen, devil_10gen, devil_15gen, devil_20gen)
names(all_vars_trimmed)<-c("bio1", "bio2", "bio3", "bio4", "bio12", "bio15", "elev", "tasveg", "roads", "rivers", "devil_5gen", "devil_10gen", "devil_15gen", "devil_20gen")

###
#trim the rasters

#first turn every raster value for lat > max(sample lat) to NA
samp_lat_max <- max(samp_latlongs$Latitude)
lat_max_row <- rowFromY(devil_5gen, samp_lat_max)
lat_max_cell <- max(cellFromRow(devil_5gen, lat_max_row))

all_vars_trimmed_maxlat <- all_vars_trimmed
all_vars_trimmed_maxlat[1:lat_max_cell] <- NA

#next turn every raster value for long > max(sample long) + buffer to NA
samp_long_min<-min(samp_latlongs$Longitude)
samp_long_max<-max(samp_latlongs$Longitude)
samp_long_edge<-samp_long_max+((samp_long_max-samp_long_min)/5)
long_edge_col<-colFromX(devil_5gen, samp_long_edge)

all_vars_trimmed_maxlong<-all_vars_trimmed_maxlat
all_vars_trimmed_maxlong[,long_edge_col:631]<-NA

#next turn every raster value for lat < min(sample lat) + buffer to NA
samp_lat_min<-min(samp_latlongs$Latitude)
samp_lat_edge<-samp_lat_min-((samp_lat_max-samp_lat_min)/5)
lat_edge_row<-rowFromY(devil_5gen, samp_lat_edge)

all_vars_trimmed_minlat<-all_vars_trimmed_maxlong
all_vars_trimmed_minlat[lat_edge_row:554,]<-NA

#crop extent of raster to remove unnecessary NAs
all_vars_trimmed_cropped<-trim(all_vars_trimmed_minlat)

#output files
writeRaster(x=all_vars_trimmed_cropped, filename="rasters_for_rga", format="GTiff")


##############################
#Make sure all samples are located on cells that have data (might not be true because 1km resolution is quite coarse)
##############################

#load in metadata for sample subset
samp_metadata<-read.table("metadata_latlong_year_subsamp", header=TRUE, sep="\t")
samp_latlongs<-samp_metadata[,c(4:5)]

#extract environmental data at sample latlongs
env_mean<-raster::extract(x=all_vars_trimmed_cropped, y=samp_latlongs, fun="mean", na.rm=TRUE, df=TRUE, layer=1, nl=14, sp=FALSE)

#identify samples that are not on raster cells due to raster resolution
na_index<-which(is.na(env_mean[,2]))
na_latlongs<-samp_latlongs[na_index,]

#convert one of the raster layers to a df
r_df<-raster::as.data.frame(x=all_vars_trimmed_cropped, xy=TRUE, na.rm=TRUE)

#calculate distance between each of the NA sample coordinates and coordinates that do have raster values
na_dist<-geodist::geodist(x=na_latlongs, y=r_df[,1:2], measure="geodesic")

#replace NA sample coordinates with coordinates of nearest raster cell with environmental values
##very minorly changes the location of 8 samples (range: 0.541 - 1.052km change)
na_dist_min_index<-c()
for (i in 1:nrow(na_latlongs)){
  na_dist_vec<-as.vector(na_dist[i,])
  na_dist_min_index[i]<-which(na_dist_vec==min(na_dist_vec))
}

na_new_pts<-r_df[na_dist_min_index,1:2]

samp_latlongs_adjusted<-samp_latlongs
for (i in 1:length(na_index)){
  samp_latlongs_adjusted[na_index[i], 1:2]<-na_new_pts[i,]
}

#write.table(data.frame(ID=samp_metadata[,1], samp_latlongs_adjusted), "latlongs_for_rga", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

##############################
#Fill NA values (generally the ocean) to a unique value because NAs are not allowed by ResistanceGA
##############################

#load environmental data
env<-stack("rasters_for_rga.tif")

#code environmental layers with convenient names for later scripting
##using non-stereotypic filenames might become unwieldly for parallel processing
names(env)<-paste("e", seq(1:14), sep="")
names(env)

#Fill NA values with -1 value
##this value is less than any value found in the non-NA data for all environmental factors
for (i in 1:dim(env)[3]){
  env_temp<-env[[i]]
  env_na<-env_temp
  env_na[which(is.na(env_na[])==TRUE)]<- -1
  
  writeRaster(env_na, paste(names(env)[i], "_na", sep=""), format="ascii")
}

