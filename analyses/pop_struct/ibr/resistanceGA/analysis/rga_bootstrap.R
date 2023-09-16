
library(ResistanceGA)
library(vcfR)
library(adegenet)
library(raster)
library(stringr)
library(tidyr)
library(ggplot2)
library(patchwork)

##############################
#Load in resistanceGA model outputs
##############################

#read in commute distance matrices
commute_dist_path <- "C:/Users/icepi/Documents/Github/stquoll_landscape_genomics/analyses/pop_struct/ibr/resistanceGA/analysis/outputs/commute_distance_matrices"
commdist_list <- list.files(path=commute_dist_path, pattern=".csv")

commdist <- list()
for (i in 1:length(commdist_list)){
  commdist[[i]] <- read.csv(paste(commute_dist_path, commdist_list[i], sep="/"), header=FALSE)
}


#read in model parameter counts (k values)
model_kval <- read.table("model_kvals", header=TRUE, sep="\t")
model_kval$mod_name[which(model_kval$mod_name=="Geographic distance")] <- "Distance"
rownames(model_kval) <- model_kval$mod_name

#match commdist and k value order
commdist_models <- str_remove(commdist_list, "_commuteDistance_distMat.csv")
model_kval_sorted <- model_kval[c(commdist_models),]

names(commdist) <- commdist_models


##############################
#Calculate DPS between individuals (pairwise genetic distances)
##############################
stq_genind <- readRDS("C:/Users/icepi/Documents/Github/stquoll_landscape_genomics/analyses/pop_struct/ibr/resistanceGA/analysis/stq_genind_rga")
stq_ps_all <- adegenet::propShared(stq_genind)
stq_dps_all <- 1-stq_ps_all


##############################
#Bootstrapping optimized models
##############################
boots <- Resist.boot(
  mod.names = commdist_models,
  dist.mat = commdist,
  n.parameters = model_kval_sorted$k,
  sample.prop = 0.70,
  iters = 20,
  n.cores=10,
  obs = 189,
  rank.method = "AICc",
  genetic.mat = stq_dps_all
)

#read in table relating coded variable names with environmental factors
env_key <- read.csv("C:/Users/icepi/Documents/Github/stquoll_landscape_genomics/analyses/pop_struct/ibr/resistanceGA/analysis/env_layers_key.csv")

#replace coded variable names with names of actual environmental factors
for (i in 1:nrow(boots)){
  layer_simplify <- str_replace_all(boots$surface[i], pattern="_na", replacement="")
  layer_simplify2 <- unlist(str_split(layer_simplify, "\\."))
  
  for (j in 1:length(layer_simplify2)){
    layer_simplify2[j] <- env_key$env_factor[which(env_key$layer==layer_simplify2[j])]
  }
  
  boots$surface[i] <- paste(layer_simplify2, collapse=" + ")
}

#if not re-running above code, load in our results below:
boots <- read.table("C:/Users/icepi/Documents/Github/stquoll_landscape_genomics/analyses/pop_struct/ibr/resistanceGA/analysis/100mod_10k_bs", header=TRUE, sep="\t")


##############################
#Plot model performance in order of increasing AICc
##############################

model_support <- boots
model_support <- model_support[order(model_support$avg.AICc),]

#anonymize models
mod_stats <- data.frame(model=seq(1:nrow(model_support)), aicc=model_support$avg.AICc, mr2=model_support$avg.R2m)

#isolate env vars included in models
model_vars <- data.frame(model=model_support[,1])

#create inclusion matrix
tempdf <- model_vars

#split models into vectors containing env var names
mod_list<-list()
for (i in 1:nrow(tempdf)){
  vars <- str_replace_all(tempdf$model[i], pattern=" \\+ ", replacement=";")
  vars<-unlist(str_split(vars, ";"))
  mod_list[[i]]<-vars
}

#isolate unique env var names
env_uniq <- unique(unlist(mod_list))
env_uniq_manual <- c("TDR", "TASVEG", "Gen5_devil", "Gen10_devil", "Gen15_devil", "Gen20_devil", "AP", "MAT", "TS", "IT", "PS", "Elev", "Rivers", "Roads", "Geographic distance")

#populate matrix (rows = variables, columns = model index) with binary values indicating whether a variable is included in a model
var_bin<-data.frame(matrix(data=NA, nrow=length(env_uniq), ncol=nrow(model_support)))
rownames(var_bin)<-env_uniq_manual
for (i in 1:length(mod_list)){
  for (j in 1:nrow(var_bin)){
    #if a particular variable (represented in the rownames of var_bin) is present in a model's variables, change cell value to 1
    if (length(intersect(rownames(var_bin)[j], mod_list[[i]]))==1){
      var_bin[j,i]<-1
    }else{
      var_bin[j,i]<-0
    }
  }
}


#reformat binary variable inclusion matrix for plotting
var_bin_plot<-data.frame(model=seq(1:length(mod_list)), t(var_bin))
colnames(var_bin_plot)[which(colnames(var_bin_plot)=="Geographic.distance")] <- "Geog.distance"
var_bin_long <- gather(var_bin_plot, env_var, inclusion, TDR:Geog.distance, factor_key=TRUE)

#revise so that geographic distance is included in every model
##ResistanceGA implicitly accounts for geographic distance because the minimum resistance value for a given surface is 1
var_bin_long$inclusion[which(var_bin_long$env_var=="Geog.distance")] <- 1

#plot AICc
aicc_plot <- ggplot(data=mod_stats, mapping=aes(x=as.factor(model), y=aicc, group=1), colour="black", size=1) +
  geom_point(size=2)+
  geom_line(size=1)+
  theme_bw()+
  labs(y="AICc")+
  theme(panel.grid.minor.x = element_blank())+
  theme(axis.text=element_text(size=12))+
  theme(axis.title=element_text(size=16))+
  theme(axis.text.x=element_blank())+
  theme(axis.title.x=element_blank())
aicc_plot

#plot marginal R-squared
mr2_plot <- ggplot(data=mod_stats, mapping=aes(x=as.factor(model), y=mr2)) +
  geom_col(fill="tan1", colour="black", size=1, alpha=0.666, width=1) +
  scale_x_discrete(expand=c(0,0))+
  scale_y_continuous(limits=c(0,0.75), expand=c(0,0))+
  theme_bw()+
  labs(y="Marginal R-squared")+
  theme(axis.text.x=element_blank())+
  theme(axis.title.x=element_blank())+
  theme(axis.text=element_text(size=12))+
  theme(axis.title=element_text(size=16))+
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
mr2_plot

#plot variable inclusion matrix
mod_tile <- ggplot(var_bin_long, aes(x = as.factor(model), y = env_var, fill = factor(inclusion))) +
  geom_tile(color="gray30", lwd=1, linetype=1)+
  scale_x_discrete(expand=c(0,0), breaks=seq(from=1, to=100, by=2))+
  scale_y_discrete(expand=c(0,0))+
  scale_fill_manual(name="Inclusion", labels=c("No", "Yes"), values=c("white", "slateblue3"))+
  theme_bw()+
  labs(x="Model (ordered by increasing AICc)", y="Environmental factor")+
  theme(axis.text=element_text(size=12))+
  theme(axis.title=element_text(size=16))+
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
mod_tile

#combine plots
aicc_plot / mr2_plot / mod_tile + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 25))


##############################
#Plot model performance in order of decreasing marginal R-squared
##############################


#load raw model performance stats
model_support_mr2 <- boots
model_support_mr2 <- model_support_mr2[order(model_support_mr2$avg.R2m, decreasing=TRUE),]

#anonymize models
mod_stats_mr2 <- data.frame(model=seq(1:nrow(model_support_mr2)), aicc=model_support_mr2$avg.AICc, mr2=model_support_mr2$avg.R2m)

#isolate env vars included in models
model_vars_mr2 <- data.frame(model=model_support_mr2[,1])

#create inclusion matrix
tempdf_mr2 <- model_vars_mr2

#split models into vectors containing env var names
mod_list_mr2 <- list()
for (i in 1:nrow(tempdf_mr2)){
  vars <- str_replace_all(tempdf_mr2$model[i], pattern=" \\+ ", replacement=";")
  vars <- unlist(str_split(vars, ";"))
  mod_list_mr2[[i]]<-vars
}


#isolate unique env var names
env_uniq_mr2 <- unique(unlist(mod_list_mr2))
env_uniq_manual_mr2 <- c("TDR", "TASVEG", "Gen5_devil", "Gen10_devil", "Gen15_devil", "Gen20_devil", "AP", "MAT", "TS", "IT", "PS", "Elev", "Rivers", "Roads", "Geographic distance")

#populate matrix (rows = variables, columns = model index) with binary values indicating whether a variable is included in a model
var_bin_mr2<-data.frame(matrix(data=NA, nrow=length(env_uniq_mr2), ncol=nrow(model_support_mr2)))
rownames(var_bin_mr2)<-env_uniq_manual_mr2
for (i in 1:length(mod_list_mr2)){
  for (j in 1:nrow(var_bin_mr2)){
    #if a particular variable (represented in the rownames of var_bin) is present in a model's variables, change cell value to 1
    if (length(intersect(rownames(var_bin_mr2)[j], mod_list_mr2[[i]]))==1){
      var_bin_mr2[j,i]<-1
    }else{
      var_bin_mr2[j,i]<-0
    }
  }
}

#reformat binary variable inclusion matrix for plotting
var_bin_plot_mr2<-data.frame(model=seq(1:length(mod_list_mr2)), t(var_bin_mr2))
colnames(var_bin_plot_mr2)[which(colnames(var_bin_plot_mr2)=="Geographic.distance")] <- "Geog.distance"
var_bin_long_mr2 <- gather(var_bin_plot_mr2, env_var, inclusion, TDR:Geog.distance, factor_key=TRUE)

#revise so that geographic distance is included in every model
var_bin_long_mr2$inclusion[which(var_bin_long_mr2$env_var=="Geog.distance")] <- 1

#plot AICc
aicc_plot_mr2 <- ggplot(data=mod_stats_mr2, mapping=aes(x=as.factor(model), y=aicc, group=1), colour="black", size=1) +
  geom_point(size=2)+
  geom_line(size=1)+
  theme_bw()+
  labs(y="AICc")+
  theme(panel.grid.minor.x = element_blank())+
  theme(axis.text=element_text(size=12))+
  theme(axis.title=element_text(size=16))+
  theme(axis.text.x=element_blank())+
  theme(axis.title.x=element_blank())
aicc_plot_mr2

#plot marginal R-squared
mr2_plot_mr2 <- ggplot(data=mod_stats_mr2, mapping=aes(x=as.factor(model), y=mr2)) +
  geom_col(fill="tan1", colour="black", size=1, alpha=0.666, width=1) +
  scale_x_discrete(expand=c(0,0))+
  scale_y_continuous(limits=c(0,0.75), expand=c(0,0))+
  theme_bw()+
  labs(y="Marginal R-squared")+
  theme(axis.text.x=element_blank())+
  theme(axis.title.x=element_blank())+
  theme(axis.text=element_text(size=12))+
  theme(axis.title=element_text(size=16))+
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
mr2_plot_mr2

#plot variable inclusion matrix
mod_tile_mr2 <- ggplot(var_bin_long_mr2, aes(x = as.factor(model), y = env_var, fill = factor(inclusion))) +
  geom_tile(color="gray30", lwd=1, linetype=1)+
  scale_x_discrete(expand=c(0,0), breaks=seq(from=1, to=100, by=2))+
  scale_y_discrete(expand=c(0,0))+
  scale_fill_manual(name="Inclusion", labels=c("No", "Yes"), values=c("white", "slateblue3"))+
  theme_bw()+
  labs(x="Model (ordered by decreasing marginal R-squared)", y="Environmental factor")+
  theme(axis.text=element_text(size=12))+
  theme(axis.title=element_text(size=16))+
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
mod_tile_mr2

#combine plots
aicc_plot_mr2 / mr2_plot_mr2 / mod_tile_mr2 + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 25))


