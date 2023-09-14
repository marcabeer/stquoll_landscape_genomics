
##########
###
#plotting in R instructions
#https://github.com/dipetkov/eems/tree/master/plotting


##############################
#Plot EEMS results
##############################

library(rEEMSplots)
library(rgdal)
#install.packages("rworldmap")
#install.packages("rworldxtra")
library("rworldmap")
library("rworldxtra")

###
#250-deme lattice

#specify path to EEMS output
eems_results <- file.path("./250deme_run1")

#specify path to export plots
plot_dir <- file.path("./250deme_run1", "250dr1")

#specify input projection
projection_none <- "+proj=longlat +datum=WGS84"

#specify output projection
projection_merc <- "+proj=merc +datum=WGS84"

#create plots
eems.plots(mcmcpath=eems_results,
           xpd=TRUE,
           plotpath=plot_dir,
           projection.in="+proj=longlat +datum=WGS84",
           projection.out=projection_merc,
           add.map=TRUE,
           longlat=TRUE,
           add.grid=TRUE,
           col.grid="gray",
           add.demes=TRUE,
           col.demes="black",
           add.colbar=TRUE,
           m.colscale=1,
           q.colscale=1,
           add.title=FALSE,
           out.png=FALSE)

###
#500-deme lattice

#specify path to EEMS output
eems_results <- file.path("./500deme_run1")

#specify path to export plots
plot_dir <- file.path("./500deme_run1", "500dr1")

#specify input projection
projection_none <- "+proj=longlat +datum=WGS84"

#specify output projection
projection_merc <- "+proj=merc +datum=WGS84"

#create plots
eems.plots(mcmcpath=eems_results,
           xpd=TRUE,
           plotpath=plot_dir,
           projection.in="+proj=longlat +datum=WGS84",
           projection.out=projection_merc,
           add.map=TRUE,
           longlat=TRUE,
           add.grid=TRUE,
           col.grid="gray",
           add.demes=TRUE,
           col.demes="black",
           add.colbar=TRUE,
           m.colscale=1,
           q.colscale=1,
           add.title=FALSE,
           out.png=FALSE)

###
#1000-deme lattice

#specify path to EEMS output
eems_results <- file.path("./1000deme_run1")

#specify path to export plots
plot_dir <- file.path("./1000deme_run1", "1000dr1")

#specify input projection
projection_none <- "+proj=longlat +datum=WGS84"

#specify output projection
projection_merc <- "+proj=merc +datum=WGS84"

#create plots
eems.plots(mcmcpath=eems_results,
           xpd=TRUE,
           plotpath=plot_dir,
           projection.in="+proj=longlat +datum=WGS84",
           projection.out=projection_merc,
           add.map=TRUE,
           longlat=TRUE,
           add.grid=TRUE,
           col.grid="gray",
           add.demes=TRUE,
           col.demes="black",
           add.colbar=TRUE,
           m.colscale=1,
           q.colscale=1,
           add.title=FALSE,
           out.png=FALSE)