setwd("D:/stq/env_data2")

#load in packages
library(raster)
library(sp)

##########
###
#download WorldClim data
url_bioclim<-"https://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_30s_bio.zip"
url_elev<-"https://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_30s_elev.zip"
url_precip<-"https://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_30s_prec.zip"
url_solar<-"https://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_30s_srad.zip"

#download.file(url_bioclim, destfile="wc2.1_30s_bio.zip")
#download.file(url_elev, destfile="wc2.1_30s_elev.zip")
#download.file(url_precip, destfile="wc2.1_30s_prec.zip")
#download.file(url_solar, destfile="wc2.1_30s_srad.zip")

#unzip("wc2.1_30s_bio.zip")
#unzip("wc2.1_30s_elev.zip")
#unzip("wc2.1_30s_prec.zip")
#unzip("wc2.1_30s_srad.zip")


##########
###
#import DFTD arrival raster
dftd<-raster("DFTD_spread_raster_fig3B.tif")
plot(dftd)

new_prj<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

dftd_reproj<-projectRaster(dftd, crs=new_prj, method="ngb")
plot(dftd_reproj)

##########
###
#import devil density rasterstack
devil_rasterstack<-stack("predictionStack_devils_1985to2035.tif")
plot(devil_rasterstack)
devil_reproj<-projectRaster(devil_rasterstack, crs=new_prj, method="ngb")
plot(devil_reproj)

##########
###
#import TASVEG raster (reclassified)
tasveg_raster_raw<-raster(x="tasveg_raster_v4.tif")
tasveg_raster_raw_unique<-unique(tasveg_raster_raw)
tasveg_raster_raw[tasveg_raster_raw==(tasveg_raster_raw_unique[2])]<-NA
tasveg_reproj<-projectRaster(tasveg_raster_raw, crs=new_prj, method="ngb")

#import higher resolution (10m) TASVEG raster


#make binary raster for each TASVEG landcover type
#binary rasters made using QGIS
#https://gis.stackexchange.com/questions/121532/how-to-reclass-a-raster-with-reclassify-grid-values-in-qgis
tasveg_100m_raw<-raster(x="tasveg_reclass_100m.tif")
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

#plot categorical raster
plot(tasveg_100m_raw)

#plot individual binary rasters
plot(tasveg_100m_mod)
plot(tasveg_100m_def)
plot(tasveg_100m_htl)
plot(tasveg_100m_mol)
plot(tasveg_100m_ngl)
plot(tasveg_100m_nef)
plot(tasveg_100m_raf)
plot(tasveg_100m_swl)
plot(tasveg_100m_scr)
plot(tasveg_100m_wef)
plot(tasveg_100m_oth)


##########
###
#import WorldClim data
#already has WGS84 CRS, so no need to reproject

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

##########
###
#attempt to overcome differing extents problem
#https://gis.stackexchange.com/questions/217082/handling-multiple-extent-problem-to-create-raster-stack-in-r

#cropping extent
#xlim=c(144.5,149), ylim=c(-44,-40.5)
#custom_extent<-raster::extent(144.5,149,-43.7,-40.5)
s_extent<-raster::extent(142,150,-45,-39)

bio1_scrop<-raster::crop(bio1, s_extent)
bio2_scrop<-raster::crop(bio2, s_extent)
bio3_scrop<-raster::crop(bio3, s_extent)
bio4_scrop<-raster::crop(bio4, s_extent)
bio12_scrop<-raster::crop(bio12, s_extent)
bio15_scrop<-raster::crop(bio15, s_extent)
elev_scrop<-raster::crop(elev, s_extent)


plot(devil_reproj[[1]])
plot(bio1_scrop)

#resampling necessary to get same extents
tasveg_resamp<-resample(x=tasveg_reproj, y=bio1_scrop, method="ngb")
dftd_resamp<-resample(x=dftd_reproj, y=bio1_scrop, method="ngb")
devil_resamp<-resample(x=devil_reproj, y=bio1_scrop, method="ngb")

#resampling retrial (11/26/2021)
bio1_resamp<-resample(x=bio1_scrop, y=devil_reproj, method="bilinear")
bio2_resamp<-resample(x=bio2_scrop, y=devil_reproj, method="bilinear")
bio3_resamp<-resample(x=bio3_scrop, y=devil_reproj, method="bilinear")
bio4_resamp<-resample(x=bio4_scrop, y=devil_reproj, method="bilinear")
bio12_resamp<-resample(x=bio12_scrop, y=devil_reproj, method="bilinear")
bio15_resamp<-resample(x=bio15_scrop, y=devil_reproj, method="bilinear")
elev_resamp<-resample(x=elev_scrop, y=devil_reproj, method="bilinear")
dftd_resamp<-resample(x=dftd_reproj, y=devil_reproj, method="bilinear")

tasveg_10m_mod_resamp<-resample(x=tasveg_10m_mod, y=devil_reproj, method="bilinear")
tasveg_10m_def_resamp<-resample(x=tasveg_10m_def, y=devil_reproj, method="bilinear")
tasveg_10m_htl_resamp<-resample(x=tasveg_10m_htl, y=devil_reproj, method="bilinear")
tasveg_10m_mol_resamp<-resample(x=tasveg_10m_mol, y=devil_reproj, method="bilinear")
tasveg_10m_ngl_resamp<-resample(x=tasveg_10m_ngl, y=devil_reproj, method="bilinear")
tasveg_10m_nef_resamp<-resample(x=tasveg_10m_nef, y=devil_reproj, method="bilinear")
tasveg_10m_raf_resamp<-resample(x=tasveg_10m_raf, y=devil_reproj, method="bilinear")
tasveg_10m_swl_resamp<-resample(x=tasveg_10m_swl, y=devil_reproj, method="bilinear")
tasveg_10m_scr_resamp<-resample(x=tasveg_10m_scr, y=devil_reproj, method="bilinear")
tasveg_10m_wef_resamp<-resample(x=tasveg_10m_wef, y=devil_reproj, method="bilinear")
tasveg_10m_oth_resamp<-resample(x=tasveg_10m_oth, y=devil_reproj, method="bilinear")

plot(tasveg_10m_htl)

all_vars_v2<-stack(bio1_resamp,
      bio2_resamp,
      bio3_resamp,
      bio4_resamp,
      bio12_resamp,
      bio15_resamp,
      elev_resamp,
      tasveg_10m_mod_resamp,
      tasveg_10m_def_resamp,
      tasveg_10m_htl_resamp,
      tasveg_10m_mol_resamp,
      tasveg_10m_ngl_resamp,
      tasveg_10m_nef_resamp,
      tasveg_10m_raf_resamp,
      tasveg_10m_swl_resamp,
      tasveg_10m_scr_resamp,
      tasveg_10m_wef_resamp,
      tasveg_10m_oth_resamp,
      dftd_resamp,
      devil_reproj)

plot(bio1_resamp)
plot(devil_reproj[[1]])
plot(tasveg_10m_mod)
plot(tasveg_10m_mod_resamp)

#resampling retrial v2 (11/26/2021)
dftd_resamp<-resample(x=dftd_reproj, y=tasveg_100m_mod, method="bilinear")
devil_resamp<-resample(x=devil_reproj, y=tasveg_100m_mod, method="bilinear")
bio1_resamp<-resample(x=bio1_scrop, y=tasveg_100m_mod, method="bilinear")
bio2_resamp<-resample(x=bio2_scrop, y=tasveg_100m_mod, method="bilinear")
bio3_resamp<-resample(x=bio3_scrop, y=tasveg_100m_mod, method="bilinear")
bio4_resamp<-resample(x=bio4_scrop, y=tasveg_100m_mod, method="bilinear")
bio12_resamp<-resample(x=bio12_scrop, y=tasveg_100m_mod, method="bilinear")
bio15_resamp<-resample(x=bio15_scrop, y=tasveg_100m_mod, method="bilinear")
elev_resamp<-resample(x=elev_scrop, y=tasveg_100m_mod, method="bilinear")

plot(devil_resamp[[1]])
plot(devil_reproj[[1]])
plot(bio1_resamp)


all_vars<-stack(bio1_resamp,
                   bio2_resamp,
                   bio3_resamp,
                   bio4_resamp,
                   bio12_resamp,
                   bio15_resamp,
                   elev_resamp,
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

##########
###

#load in rasterstack
all_vars<-stack("all_vars_100m.tif")


#create geodesic buffers
#https:/gis.stackexchange.com/questions/250389/euclidean-and-geodesic-buffering-in-r
#devtools::install_github("valentinitnelav/geobuffer")
library(geobuffer)

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




#in case of NA values:
#https://stackoverflow.com/questions/26652629/extracting-a-value-from-a-raster-for-a-specific-point-based-on-the-closest-cell

#obtain latlongs of points with NA values
na_index<-which(is.na(env_mean[,21]))
na_latlongs<-latlongs[na_index,]

#https://stackoverflow.com/questions/26652629/extracting-a-value-from-a-raster-for-a-specific-point-based-on-the-closest-cell
move.points <- function(r, pts, spatial=FALSE) {
  require(raster)
  require(sp)
  
  if (is(pts, 'SpatialPoints')) pts <- coordinates(pts)
  if (is(!r, 'Raster')) r <- raster(r)
  
  loc <- colSums(sapply(pts[, 1], '>', bbox(r)[1, ])) * 3 + 
    colSums(sapply(pts[, 2], '>', bbox(r)[2, ]))
  
  L <- split(as.data.frame(pts), loc)
  
  new.pts <- lapply(names(L), function(x) {
    switch(x, 
           '0' = xyFromCell(r, ncell(r) - ncol(r) + 1)[rep(1, nrow(L[[x]])), ],
           '1' = xyFromCell(r, cellFromXY(r, cbind(xmin(r), L[[x]][, 2]))),
           '2' = xyFromCell(r, 1)[rep(1, nrow(L[[x]])), ],
           '3' = xyFromCell(r, cellFromXY(r, cbind(L[[x]][, 1], ymin(r)))),
           '4' = {
             xy <- as.matrix(L[[x]])
             dimnames(xy) <- list(NULL, c('x', 'y'))
             xy
           },
           '5' = xyFromCell(r, cellFromXY(r, cbind(L[[x]][, 1], ymax(r)))),
           '6' = xyFromCell(r, ncell(r))[rep(1, nrow(L[[x]])), ],
           '7' = xyFromCell(r, cellFromXY(r, cbind(xmax(r), L[[x]][, 2]))),
           '8' = xyFromCell(r, ncol(r))[rep(1, nrow(L[[x]])), ]
    )
  })
  
  new.pts <- unsplit(mapply(function(x, y) {
    row.names(x) <- row.names(y)
    as.data.frame(x)
  }, new.pts, L, SIMPLIFY=FALSE), loc)
  
  colnames(new.pts) <- colnames(pts)
  if(isTRUE(spatial)) new.pts <- SpatialPoints(new.pts)  
  return(new.pts)
}

plot(all_vars[[21]])
r<-all_vars[[21]]

na_pts_moved<-move.points(r=r, pts=na_latlongs, spatial=FALSE)
na_latlong_buffer<-geobuffer::geobuffer_pts(xy=na_pts_moved, dist_m=(rad_and*1000), step_dg=10)
na_env_mean<-raster::extract(x=all_vars, y=na_latlong_buffer, fun="mean", na.rm=TRUE, df=TRUE, layer=1, nl=70, sp=FALSE)


#https://stackoverflow.com/questions/27562076/if-raster-value-na-search-and-extract-the-nearest-non-na-pixel
#a function for sampling
sample_raster_NA <- function(r, xy){
  apply(X = xy, MARGIN = 1, 
        FUN = function(xy) r@data@values[which.min(replace(distanceFromPoints(r, xy), is.na(r), NA))])
  
}

na_new_pts<-sample_raster_NA(r=r, xy=na_latlongs)


#try converting raster to df, omitting NA values
r_df<-raster::as.data.frame(x=r, xy=TRUE, na.rm=TRUE)

library(geodist)

na_dist<-geodist::geodist(x=na_latlongs, y=r_df[,1:2], measure="geodesic")

na_dist_min_index<-c()
for (i in 1:nrow(na_latlongs)){
  na_dist_vec<-as.vector(na_dist[i,])
  na_dist_min_index[i]<-which(na_dist_vec==min(na_dist_vec))
}

na_new_pts<-r_df[na_dist_min_index,1:2]

na_new_buffer<-geobuffer::geobuffer_pts(xy=na_new_pts, dist_m=(rad_and*1000), step_dg=10)

na_env_mean<-raster::extract(x=all_vars, y=na_new_buffer, fun="mean", na.rm=TRUE, df=TRUE, layer=1, nl=70, sp=FALSE)

#input new non-na values into original env_mean df
for (i in 1:length(na_index)){
  env_mean[na_index[i], 2:ncol(na_env_mean)]<-na_env_mean[i,-1]
}

write.table(x=env_mean, file="env_data_nafilled_extracted", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
