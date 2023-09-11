
library(gdm)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(cowplot)

##############################
#Load in data
##############################

#read in jittered latlongs
latlongs_jit <- read.table("latlongs_jit", header=TRUE, sep="\t")

#read in collection dates
dates <- read.table("stq_dates_yeardec", sep="\t", header=TRUE)

#read in env_r80 data
env <- read.table("env_r80_yeardec", header=TRUE, sep="\t")
env <- data.frame(env, year=dates$yeardec)

#read in Dps matrix (individual pairwise genetic distances)
stq_dps <- read.table("stq_dps", header=FALSE, sep="\t")


##############################
#reformat data for GDM
##############################

#create site variable
site <- 1:nrow(latlongs_jit)

#mimic "envTab" df from GDM tutorial
geoenv <- data.frame(site, latlongs_jit, env)

#attach site column to DPS matrix
stq_dps_site <- data.frame(site, stq_dps)

#use formatesitepair() from gdm package
gdmTab.dis<-formatsitepair(bioData=stq_dps_site, bioFormat=3, XColumn="Long", YColumn="Lat", predData=geoenv, siteColumn="site")


##############################
#run GDM matrix permutation
##############################

#matrix permutation to assess variable importance
gdm.1.matperm <- gdm.varImp(spTable=gdmTab.dis, geo=TRUE, fullModelOnly=FALSE, nPerm=500, parallel=TRUE, cores=5, sampleSites=1, sampleSitePairs=1, outFile="dps_gdm")



##############################
#further evaluate GDM outputs
##############################
load("dps_gdm.Rdata")


#plot deviance explained of models
gdm_perm_res<-data.frame(model=seq(1:19), perc_dev=outObject[[1]][2,])

p <- ggplot()+
  geom_line(data=gdm_perm_res, mapping=aes(x=model, y=perc_dev))+
  geom_point(data=gdm_perm_res, mapping=aes(x=model, y=perc_dev))+
  geom_abline(slope=0, intercept=outObject[[1]][2,1], colour="red")+
  ylim(35,60)+
  scale_x_continuous(breaks=seq(0:20))+
  theme_bw()+
  theme(panel.grid.minor.x=element_blank())+
  labs(x="Model", y="Deviance explained (%)")

p_zoomed <- ggplot()+
  geom_line(data=gdm_perm_res, mapping=aes(x=model, y=perc_dev))+
  geom_point(data=gdm_perm_res, mapping=aes(x=model, y=perc_dev))+
  geom_abline(slope=0, intercept=outObject[[1]][2,1], colour="red")+
  ylim(55,57)+
  scale_x_continuous(breaks=seq(0:20))+
  theme_bw()+
  theme(panel.grid.minor.x=element_blank())+
  labs(x="Model", y="Deviance explained (%)")

cowplot::plot_grid(p, p_zoomed, ncol=1, nrow=2)


#based on deviance explained:
##largest model with only significant effects is Model 10 (i.e., full model - 9 covariates) (which has only 0.36025 % lower deviance explained than full model)
#Model variable importance and significance df
mod10_df <- data.frame(imp=outObject[[2]][,11], p=outObject[[3]][,11], nperm=outObject[[4]][,11])

write.csv(mod10_df, "mod10_df.csv")


##############################
#re-fit best model to obtain splines
##############################

#load raw data
#read in jittered latlongs
latlongs_jit <- read.table("latlongs_jit", header=TRUE, sep="\t")

#read in collection dates
dates <- read.table("stq_dates_yeardec", sep="\t", header=TRUE)

#read in environmental data
env <- read.table("env_r80_yeardec", header=TRUE, sep="\t")
env <- data.frame(env, year=dates$yeardec)

###
#format data for GDM

#create site variable
site <- 1:nrow(latlongs_jit)

#mimic "envTab" df from GDM tutorial
geoenv <- data.frame(site, latlongs_jit, env)

#attach site column to DPS matrix
stq_dps_site <- data.frame(site, stq_dps)

#use formatesitepair() from gdm package
gdmTab.dis <- formatsitepair(bioData=stq_dps_site, bioFormat=3, XColumn="Long", YColumn="Lat", predData=geoenv, siteColumn="site")

###
#subset data to remove extraneous variables
m00 <- gdmTab.dis
m10 <- subset(x=m00, select=-c(s1.gen15_devil, s2.gen15_devil, s1.tasveg_100m_swl, s2.tasveg_100m_swl, s1.tasveg_100m_ngl, s2.tasveg_100m_ngl, s1.wc2.1_30s_elev, s2.wc2.1_30s_elev, s1.tasveg_100m_scr, s2.tasveg_100m_scr, s1.tasveg_100m_mol, s2.tasveg_100m_mol, s1.tasveg_100m_oth, s2.tasveg_100m_oth, s1.tasveg_100m_def, s2.tasveg_100m_def, s1.tasveg_100m_nef, s2.tasveg_100m_nef, s1.year, s2.year))

#re-fit model and plot response curves
m10_refit <- gdm(m10, geo=TRUE)



#####
#plotting splines for re-fit model

#plot splines
length(m10_refit$predictors)
plot(m10_refit, plot.layout=c(4,3))

#extract spline plots
m10_splineDat <- isplineExtract(m10_refit)

#normalize x-axes
m10_spline_x <- scale(m10_splineDat$x, center=TRUE, scale=TRUE)
geog <- data.frame(x=m10_spline_x[,1], y=m10_splineDat$y[,1], var="Geographic distance (GD)")
bio3 <- data.frame(x=m10_spline_x[,2], y=m10_splineDat$y[,2], var="Isothermality (IT)")
bio4 <- data.frame(x=m10_spline_x[,3], y=m10_splineDat$y[,3], var="Temperature seasonality (TS)")
bio12 <- data.frame(x=m10_spline_x[,4], y=m10_splineDat$y[,4], var="Annual precipitation (AP)")
bio15 <- data.frame(x=m10_spline_x[,4], y=m10_splineDat$y[,5], var="Precipitation seasonality (PS)")
tasveg_mod <- data.frame(x=m10_spline_x[,5], y=m10_splineDat$y[,6], var="Modified landcover (MOD)")
tasveg_raf <- data.frame(x=m10_spline_x[,7], y=m10_splineDat$y[,7], var="Rainforest (RAF)")
tasveg_wef <- data.frame(x=m10_spline_x[,8], y=m10_splineDat$y[,8], var="Wet eucalypt forest (WEF)")
gen_dftd <- data.frame(x=m10_spline_x[,9], y=m10_splineDat$y[,9], var="Generations of DFTD presence (DFTD)")


spline_df <- rbind(geog, bio3, bio4, bio12, bio15, tasveg_mod, tasveg_raf, tasveg_wef, gen_dftd)

spline_plot <- ggplot(data=spline_df, aes(x=x, y=y, group=var, color=var)) +
  geom_line(size=1)+
  scale_color_brewer(type="qual", palette="Paired", name="")+
  labs(x="Predictor dissimilarity", y="Genetic dissimilarity")+
  theme(legend.position="top")+
  theme_bw()
spline_plot


#####
#variance partitioning for re-fit model

#variable sets
varSet_env <- list(Env=c("wc2.1_30s_bio_3", "wc2.1_30s_bio_4", "wc2.1_30s_bio_12", "wc2.1_30s_bio_15", "tasveg_100m_mod", "tasveg_100m_raf", "tasveg_100m_wef", "gen_dftd"))
varSet_bio3 <- list(Bio3="wc2.1_30s_bio_3", Env=c("wc2.1_30s_bio_4", "wc2.1_30s_bio_12", "wc2.1_30s_bio_15", "tasveg_100m_mod", "tasveg_100m_raf", "tasveg_100m_wef", "gen_dftd"))
varSet_bio4 <- list(Bio4="wc2.1_30s_bio_4", Env=c("wc2.1_30s_bio_3", "wc2.1_30s_bio_12", "wc2.1_30s_bio_15", "tasveg_100m_mod", "tasveg_100m_raf", "tasveg_100m_wef", "gen_dftd"))
varSet_bio12 <- list(Bio12="wc2.1_30s_bio_12", Env=c("wc2.1_30s_bio_3", "wc2.1_30s_bio_4", "wc2.1_30s_bio_15", "tasveg_100m_mod", "tasveg_100m_raf", "tasveg_100m_wef", "gen_dftd"))
varSet_bio15 <- list(Bio15="wc2.1_30s_bio_15", Env=c("wc2.1_30s_bio_3", "wc2.1_30s_bio_4", "wc2.1_30s_bio_12", "tasveg_100m_mod", "tasveg_100m_raf", "tasveg_100m_wef", "gen_dftd"))
varSet_mod <- list(TASVEG_mod="tasveg_100m_mod", Env=c("wc2.1_30s_bio_3", "wc2.1_30s_bio_4", "wc2.1_30s_bio_12", "wc2.1_30s_bio_15", "tasveg_100m_raf", "tasveg_100m_wef", "gen_dftd"))
varSet_raf <- list(TASVEG_raf="tasveg_100m_raf", Env=c("wc2.1_30s_bio_3", "wc2.1_30s_bio_4", "wc2.1_30s_bio_12", "wc2.1_30s_bio_15", "tasveg_100m_mod", "tasveg_100m_wef", "gen_dftd"))
varSet_wef <- list(TASVEG_wef="tasveg_100m_wef", Env=c("wc2.1_30s_bio_3", "wc2.1_30s_bio_4", "wc2.1_30s_bio_12", "wc2.1_30s_bio_15", "tasveg_100m_mod", "tasveg_100m_raf", "gen_dftd"))
varSet_dftd <- list(DFTD="gen_dftd", Env=c("wc2.1_30s_bio_3", "wc2.1_30s_bio_4", "wc2.1_30s_bio_12", "wc2.1_30s_bio_15", "tasveg_100m_mod", "tasveg_100m_raf", "tasveg_100m_wef"))

#deviance partititioning
dev_env <- gdm.partition.deviance(sitePairTable=m10, varSets=varSet_env, partSpace=TRUE)
dev_bio3 <- gdm.partition.deviance(sitePairTable=m10, varSets=varSet_bio3, partSpace=TRUE)
dev_bio4 <- gdm.partition.deviance(sitePairTable=m10, varSets=varSet_bio4, partSpace=TRUE)
dev_bio12 <- gdm.partition.deviance(sitePairTable=m10, varSets=varSet_bio12, partSpace=TRUE)
dev_bio15 <- gdm.partition.deviance(sitePairTable=m10, varSets=varSet_bio15, partSpace=TRUE)
dev_mod <- gdm.partition.deviance(sitePairTable=m10, varSets=varSet_mod, partSpace=TRUE)
dev_raf <- gdm.partition.deviance(sitePairTable=m10, varSets=varSet_raf, partSpace=TRUE)
dev_wef <- gdm.partition.deviance(sitePairTable=m10, varSets=varSet_wef, partSpace=TRUE)
dev_dftd <- gdm.partition.deviance(sitePairTable=m10, varSets=varSet_dftd, partSpace=TRUE)

#determine total deviance that is explained by the environment
dev_env$DEVIANCE[4]

#determine deviance that is uniquely explained by the individual environmental factors
env_part_unconfounded <- sum(c(dev_bio3$DEVIANCE[9],
                               dev_bio4$DEVIANCE[9],
                               dev_bio12$DEVIANCE[9],
                               dev_bio15$DEVIANCE[9],
                               dev_mod$DEVIANCE[9],
                               dev_raf$DEVIANCE[9],
                               dev_wef$DEVIANCE[9],
                               dev_dftd$DEVIANCE[9]))

#determine deviance that is explained by either geography or the environment
geo_env_confounded <- dev_env$DEVIANCE[3] - (dev_env$DEVIANCE[5] + dev_env$DEVIANCE[4])

#create dataframe for variance partitioning visualization
dev_part_df <- data.frame(FACTOR=c("Bio3", "Bio4", "Bio12", "Bio15", "TASVEG_MOD", "TASVEG_RAF", "TASVEG_WEF", "DFTD", "ENV_CONFOUND", "GEOG", "GEOG_ENV_CONFOUND", "UNEXPLAINED"), DEVIANCE_EXP=c(dev_bio3$DEVIANCE[9],
                                                                                                                                                                                                   dev_bio4$DEVIANCE[9],
                                                                                                                                                                                                   dev_bio12$DEVIANCE[9],
                                                                                                                                                                                                   dev_bio15$DEVIANCE[9],
                                                                                                                                                                                                   dev_mod$DEVIANCE[9],
                                                                                                                                                                                                   dev_raf$DEVIANCE[9],
                                                                                                                                                                                                   dev_wef$DEVIANCE[9],
                                                                                                                                                                                                   dev_dftd$DEVIANCE[9],
                                                                                                                                                                                                   dev_env$DEVIANCE[4]-env_part_unconfounded,
                                                                                                                                                                                                   dev_env$DEVIANCE[5],
                                                                                                                                                                                                   geo_env_confounded,
                                                                                                                                                                                                   dev_env$DEVIANCE[6]),
                          FACTOR_TYPE=c("Climate", "Climate", "Climate", "Climate", "Landcover", "Landcover", "Landcover", "Biotic", "ENV_CONFOUND", "GEOG", "GEOG_ENV_CONFOUND", "UNEXPLAINED"))

dev_part_df <- dev_part_df [c(2,4,3,1,7,6,5,8,9,11,10,12),]
dev_part_df$FACTOR <- factor(dev_part_df$FACTOR, levels=dev_part_df$FACTOR)


#explained variance only
dev_part_df_exp <- dev_part_df[-12,]
dev_part_df_exp$percent <- dev_part_df_exp$DEVIANCE_EXP / sum(dev_part_df_exp$DEVIANCE_EXP)*100
dev_part_df_exp$FACTOR <- factor(dev_part_df_exp$FACTOR, levels=c("Bio4", "Bio15", "Bio12", "Bio3", "TASVEG_WEF", "TASVEG_RAF", "TASVEG_MOD", "DFTD", "ENV_CONFOUND", "GEOG", "GEOG_ENV_CONFOUND"))

#summarize by factor type
dev_part_df_summ <- summarise(dev_part_df_exp, )
dev_part_df_summ <- dev_part_df %>%
  group_by(FACTOR_TYPE) %>%
  summarise(DEV_EXP_SUM = sum(DEVIANCE_EXP))
dev_part_df_summ<-dev_part_df_summ[c(2,6,1,3,5,4,7),]
dev_part_df_summ$FACTOR_TYPE <- factor(dev_part_df_summ$FACTOR_TYPE, levels=dev_part_df_summ$FACTOR_TYPE)

#check that % deviance adds up to 100
sum(dev_part_df$DEVIANCE_EXP)

###
#create pie charts

#set color palettes

pal_a90 <- c("#88CCEE90", "#CC667790", "#DDCC7790", "#11773390", "#33228890", "#AA449990", 
             "#44AA9990", "#99993390", "#88225590", "#66110090", "#6699CC90", "#88888890")

pal_general <- c("#88CCEE", "#332288", "#999933", "#88225598", "#661100", "#6699CC", "#888888")

pal_exp <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
             "#44AA99", "#999933", "#88225598", "#6699CC","#661100")


#pie chart for all individual environmental factors
df2 <- dev_part_df %>% 
  mutate(csum = rev(cumsum(rev(DEVIANCE_EXP))), 
         pos = DEVIANCE_EXP/2 + lead(csum, 1),
         pos = if_else(is.na(pos), DEVIANCE_EXP/2, pos))

summ_df2 <- dev_part_df_summ %>% 
  mutate(csum = rev(cumsum(rev(DEV_EXP_SUM))), 
         pos = DEV_EXP_SUM/2 + lead(csum, 1),
         pos = if_else(is.na(pos), DEV_EXP_SUM/2, pos))

pie_specific <- ggplot(dev_part_df, aes(x="", y=DEVIANCE_EXP, fill=FACTOR))+
  geom_bar(stat="identity", width=1, colour="black", size=0.75)+
  scale_fill_manual(values=pal_a90)+
  coord_polar("y", start=0, direction=-1)+
  geom_label_repel(data = df2,
                   aes(y = pos, x=1.5, label = paste0(DEVIANCE_EXP, "%")),
                   size = 7.5, nudge_x = 1.5, show.legend = FALSE, force=20, colour="black", box.padding=0.5, direction="both") +
  theme_bw()+
  theme(axis.text.x=element_blank(), axis.title.x=element_blank())+
  theme(axis.text.y=element_blank(), axis.title.y=element_blank())+
  theme(axis.ticks = element_blank())+
  theme(panel.border = element_blank())+
  theme(panel.grid=element_blank())
pie_specific

#group individual env factors into climate, landcover, and biotic factors
pie_general <- ggplot(summ_df2, aes(x="", y=DEV_EXP_SUM, fill=FACTOR_TYPE))+
  geom_bar(stat="identity", width=1, colour="black")+
  scale_fill_manual(name="Factor", labels=c("Climate", "Landcover", "DFTD", "Env. confound", "Geog.-Env. confound,", "Geog.", "Unexplained"), values=pal_general)+
  coord_polar("y", start=0, direction=-1)+
  theme_bw()+
  theme(axis.text.x=element_blank(), axis.title.x=element_blank())+
  theme(axis.text.y=element_blank(), axis.title.y=element_blank())+
  theme(axis.ticks = element_blank())+
  theme(panel.border = element_blank())+
  theme(panel.grid=element_blank())
pie_general


#now consider only explained deviance
pie_exp <- ggplot(dev_part_df_exp, aes(x="", y=percent, fill=FACTOR))+
  geom_bar(stat="identity", width=1, colour="black")+
  scale_fill_manual(name="Factor",
                    labels=c("Temp. seaonality", "Precip. seasonality", "Annual precip.", "Isothermality", "% WEF", "% RAF", "% MOD", "DFTD", "Env. confound", "Geog.", "Geog.-Env. confound"),
                    values=pal_exp)+
  coord_polar("y", start=0, direction=-1)+
  theme_bw()+
  theme(axis.text.x=element_blank(), axis.title.x=element_blank())+
  theme(axis.text.y=element_blank(), axis.title.y=element_blank())+
  theme(axis.ticks = element_blank())+
  theme(panel.border = element_blank())+
  theme(panel.grid=element_blank())
pie_exp


