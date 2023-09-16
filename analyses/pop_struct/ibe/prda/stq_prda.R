library(vegan)
library(adegenet)
library(spdep)
library(adespatial)
library(splancs)

##############################
#conduct sPCA using adegenet
##############################

#read genind of quoll genotypes
stq <- readRDS("stq_genind")

#Read XY coordinates, they must all be unique, use the jitter function#
XY <- read.table(file="latlongs", header=T, sep="\t")

# Matrix needs diff name than coordinates#
xy.jitter <- matrix(ncol=2, nrow=345)

#The XY[,1] is  adding noise to the first column (i.e., X coordinates) and XY[,2] is the second column/Y coordinates#
xy.jitter[,1] <- jitter(XY[,1], a=0.01)
xy.jitter[,2] <- jitter(XY[,2], a=0.01)

#sPCA function, nfposi is the number of positive axes to retain. Take 2 positive and one negative axis#
mySpca <- spca.genind(stq, xy=xy.jitter, type=5, d1=0, d2=1.99, nfposi=2, nfnega=1)

#if not re-running above code, load in our results
mySpca <- readRDS("spca_results.gz")

#Scores for resultant sPCA, these are the response variables in the RDA. Can copy and paste them from R into a text document, or write a table in R#
spca_scores <- data.frame(id=row.names(mySpca$li), mySpca$li)
#write.table(spca_scores, "spca_scores", row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")


##############################
#load environmental data
##############################

#load in env data
##removed col 10 because zero variance
env_all <- read.table(file="env_all_yeardec", header=T, sep="\t")[,-10]

#load in collection years
dates <- read.table("stq_dates_yeardec", sep="\t", header=TRUE)

#format data
env_new <- data.frame(env_all, year=dates$yeardec)
colnames(env_new) <- c("bio1", "bio2", "bio3", "bio4", "bio12", "bio15", "elev", "tasveg_mod", "tasveg_def", "tasveg_mol", "tasveg_ngl", "tasveg_nef", "tasveg_raf", "tasveg_swl", "tasveg_scr", "tasveg_wef", "tasveg_oth", "gen_dftd", "gen5_devil", "gen10_devil", "gen15_devil", "year")


##############################
#Partial redundancy analysis (pRDA)
##############################

#bind genetic data, latlongs, and environmental data together
stq_env_df<-data.frame(spca_scores[,1:3], XY, env_new)


###
#Variance inflation factor process

#remove variable with max VIF one at a time until all remaining variables have VIF<10
try <- rda(stq_env_df[,2]+stq_env_df[,3]~stq_env_df[,6]+stq_env_df[,7]+stq_env_df[,8]+stq_env_df[,9]+stq_env_df[,10]+stq_env_df[,11]+stq_env_df[,12]+stq_env_df[,13]+stq_env_df[,14]+stq_env_df[,15]+stq_env_df[,16]+stq_env_df[,17]+stq_env_df[,18]+stq_env_df[,19]+stq_env_df[,20]+stq_env_df[,21]+stq_env_df[,23]+stq_env_df[,24]+stq_env_df[,25]+stq_env_df[,26]+Condition(stq_env_df[,4]+stq_env_df[,5]+stq_env_df[,27]))
vif.cca(try)
max(vif.cca(try), na.rm=TRUE)

try2 <- rda(stq_env_df[,2]+stq_env_df[,3]~stq_env_df[,6]+stq_env_df[,8]+stq_env_df[,9]+stq_env_df[,10]+stq_env_df[,11]+stq_env_df[,12]+stq_env_df[,13]+stq_env_df[,14]+stq_env_df[,15]+stq_env_df[,16]+stq_env_df[,17]+stq_env_df[,18]+stq_env_df[,19]+stq_env_df[,20]+stq_env_df[,21]+stq_env_df[,23]+stq_env_df[,24]+stq_env_df[,25]+stq_env_df[,26]+Condition(stq_env_df[,4]+stq_env_df[,5]+stq_env_df[,27]))
vif.cca(try2)
max(vif.cca(try2), na.rm=TRUE)

try3 <- rda(stq_env_df[,2]+stq_env_df[,3]~stq_env_df[,8]+stq_env_df[,9]+stq_env_df[,10]+stq_env_df[,11]+stq_env_df[,12]+stq_env_df[,13]+stq_env_df[,14]+stq_env_df[,15]+stq_env_df[,16]+stq_env_df[,17]+stq_env_df[,18]+stq_env_df[,19]+stq_env_df[,20]+stq_env_df[,21]+stq_env_df[,23]+stq_env_df[,24]+stq_env_df[,25]+stq_env_df[,26]+Condition(stq_env_df[,4]+stq_env_df[,5]+stq_env_df[,27]))
vif.cca(try3)
max(vif.cca(try3), na.rm=TRUE)

try4 <- rda(stq_env_df[,2]+stq_env_df[,3]~stq_env_df[,8]+stq_env_df[,9]+stq_env_df[,10]+stq_env_df[,11]+stq_env_df[,12]+stq_env_df[,13]+stq_env_df[,14]+stq_env_df[,15]+stq_env_df[,16]+stq_env_df[,17]+stq_env_df[,18]+stq_env_df[,19]+stq_env_df[,20]+stq_env_df[,21]+stq_env_df[,23]+stq_env_df[,24]+stq_env_df[,25]+Condition(stq_env_df[,4]+stq_env_df[,5]+stq_env_df[,27]))
vif.cca(try4)
max(vif.cca(try4), na.rm=TRUE)

try5 <- rda(stq_env_df[,2]+stq_env_df[,3]~stq_env_df[,8]+stq_env_df[,9]+stq_env_df[,10]+stq_env_df[,11]+stq_env_df[,12]+stq_env_df[,13]+stq_env_df[,14]+stq_env_df[,15]+stq_env_df[,16]+stq_env_df[,17]+stq_env_df[,18]+stq_env_df[,19]+stq_env_df[,20]+stq_env_df[,21]+stq_env_df[,23]+stq_env_df[,25]+Condition(stq_env_df[,4]+stq_env_df[,5]+stq_env_df[,27]))
vif.cca(try5)
max(vif.cca(try5), na.rm=TRUE)

try6 <- rda(stq_env_df[,2]+stq_env_df[,3]~stq_env_df[,8]+stq_env_df[,9]+stq_env_df[,10]+stq_env_df[,11]+stq_env_df[,12]+stq_env_df[,14]+stq_env_df[,15]+stq_env_df[,16]+stq_env_df[,17]+stq_env_df[,18]+stq_env_df[,19]+stq_env_df[,20]+stq_env_df[,21]+stq_env_df[,23]+stq_env_df[,25]+Condition(stq_env_df[,4]+stq_env_df[,5]+stq_env_df[,27]))
vif.cca(try6)
max(vif.cca(try6), na.rm=TRUE)

try7 <- rda(stq_env_df[,2]+stq_env_df[,3]~stq_env_df[,8]+stq_env_df[,9]+stq_env_df[,10]+stq_env_df[,12]+stq_env_df[,14]+stq_env_df[,15]+stq_env_df[,16]+stq_env_df[,17]+stq_env_df[,18]+stq_env_df[,19]+stq_env_df[,20]+stq_env_df[,21]+stq_env_df[,23]+stq_env_df[,25]+Condition(stq_env_df[,4]+stq_env_df[,5]+stq_env_df[,27]))
vif.cca(try7)
#Max VIF is 11.45 and is Longitude - decide to keep it
max(vif.cca(try7), na.rm=TRUE)


###
#pRDA modeling

#establish null model (only geography and sampling year)
mod0_geog <- rda(stq_env_df[,2]+stq_env_df[,3] ~ 1+Condition(stq_env_df[,4]+stq_env_df[,5]+stq_env_df[,27]))

#establish full model
mod1 <- try7

#carry out model selection
sel3 <- ordiR2step(mod0_geog, scope=formula(mod1), perm.max = 10000)

#if not re-running above code, load in our results
sel3 <- readRDS("prda_sel3_wtime.gz")

#Analysis variance of final selection model
RDAc_sel3_marg <- anova.cca(sel3, by="margin")
rownames(RDAc_sel3_marg) <- c(colnames(stq_env_df)[c(12,10,20,15,23,8,19,14,9)], "Residual")
RDAc_sel3_marg

#get adjusted Rsquared of final model
RsquareAdj(sel3)
