
library(ResistanceGA)
library(vcfR)
library(adegenet)
library(raster)


##############################
#Load in data
##############################

#read in sample latlongs
latlongs_raw <- read.table("latlongs_for_rga", header=TRUE, sep="\t")
latlongs <- as.matrix(latlongs_raw[,-1])

#load genetic data
stq_genind <- readRDS("stq_genind_rga")

##############################
#Calculate genetic distances
##############################

#calculate proportion of shared alleles
stq_ps_all <- adegenet::propShared(stq_genind)
stq_ps_all[upper.tri(stq_ps_all, diag=TRUE)] <- NA

#calculate Dps (i.e., 1 - proportion of shared alleles)
stq_dps_all <- 1-stq_ps_all
stq_dps_vec <- as.vector(stq_dps_all)
stq_dps_vec <- stq_dps_vec[!is.na(stq_dps_vec)]


##############################
#Prepare data for gdist usage in resistanceGA
##############################
stq_GD_prep <-  gdist.prep(n.Pops = 189, response = stq_dps_vec, samples = latlongs, longlat=FALSE, method="commuteDistance")


##############################
#Run resistanceGA on individual environmental factors
##############################

ascii_dir <- "C:/Users/icepi/Documents/GitHub/stquoll_landscape_genomics/analyses/pop_struct/ibr/resistanceGA/rerun/env_layers_all/"
ssoptim_dir <- "C:/Users/icepi/Documents/GitHub/stquoll_landscape_genomics/analyses/pop_struct/ibr/resistanceGA/rerun/ssoptim_results/"

ssoptim_trans <- rep(list("A"), 14)
ssoptim_gaprep <- GA.prep(ASCII.dir = ascii_dir, Results.dir = ssoptim_dir, parallel = 7, select.trans=ssoptim_trans, maxiter=1)

SS_optim(gdist.inputs = stq_GD_prep, GA.inputs=ssoptim_gaprep)

print("completed univariate model optimization")


##############################
#Run resistanceGA on full model (all environmental factors)
##############################
ascii_dir <- "C:/Users/icepi/Documents/GitHub/stquoll_landscape_genomics/analyses/pop_struct/ibr/resistanceGA/rerun/env_layers_all/"
fulloptim_dir <- "C:/Users/icepi/Documents/GitHub/stquoll_landscape_genomics/analyses/pop_struct/ibr/resistanceGA/rerun/fullmodel_results/"


fulloptim_trans <- rep(list("A"), 14)
fulloptim_gaprep <- GA.prep(ASCII.dir = ascii_dir, Results.dir = fulloptim_dir, parallel = 10, select.trans=fulloptim_trans)
MS_optim(gdist.inputs = stq_GD_prep, GA.inputs = ga_input_full)

print("completed full model optimization")


##############################
#Run resistanceGA on combos of 2 and 3 environmental factors from the reduced set
##############################
msoptim_ascii_dir <- "C:/Users/icepi/Documents/GitHub/stquoll_landscape_genomics/analyses/pop_struct/ibr/resistanceGA/rerun/env_layers_reduced/"
msoptim_dir <- "C:/Users/icepi/Documents/GitHub/stquoll_landscape_genomics/analyses/pop_struct/ibr/resistanceGA/rerun/allcomb/"

msoptim_trans <- rep(list("A"), 8)
msoptim_gaprep <- GA.prep(ASCII.dir=msoptim_ascii_dir, Results.dir = "all.comb", parallel=10, select.trans=msoptim_trans)

all_comb(gdist.inputs = stq_GD_prep, GA.inputs = msoptim_gaprep, results.dir=msoptim_dir, max.combination=c(2,3), iters=1, sample.prop=1)


##############################
#Bootstrapping
##############################

#requires creation of inputs from runs of resistanceGA carried out in previous lines of this script

#mod_names is a list of models
#mat_list is a list of cost distance matrices, named by model
#n_param is a vector of parameter counts for each model
#stq_dps is the matrix of DPS values (genetic distances)


test<-Resist.boot(
  mod.names = mod_names,
  dist.mat = mat_list,
  n.parameters = n_param,
  sample.prop = 0.70,
  iters = 10000,
  n.cores=10,
  obs = 189,
  rank.method = "AICc",
  genetic.mat = stq_dps_all
)


