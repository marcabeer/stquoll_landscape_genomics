
library(pegas)
library(ggplot2)
library(raster)
library(rgdal)
library(LEA)
library(rnaturalearth)
library(rnaturalearthdata)
library(RColorBrewer)
library(ggpubr)
library(vegan)
library(qvalue)
library(robust)
library(WMDB)
library(ggVennDiagram)
library(cowplot)
library(corrplot)
library(rgeos)
library(vcfR)
library(tseries)
library(ggrepel)


##############################
#Load in genetic data and metadata
##############################

#load in imputed genotype matrix (individual genotypes, not population allele frequencies)
imputed_geno <- read.matrix("stq.recode.lfmm_imputed.lfmm", header=FALSE, sep=" ")

#get locus names
locus_info <- read.table("locus_info", header=TRUE, sep="\t")
locus_ids <- paste(as.character(locus_info[,1]), as.character(locus_info[,2]), sep="-")

#import individual PC1 and PC2 scores
pc_scores <- read.table("imputed_pca_scores", sep="\t", header=TRUE)

#import latlongs
latlongs <- read.table("latlongs", sep="\t", header=TRUE)

#import environmental data
##column 10 removed because environmental factor is invariant across sites
env_all <- read.table("env_all_yeardec", sep="\t", header=TRUE)[,-10]

#import collection year
year <- read.table("stq_dates_yeardec", sep="\t", header=TRUE)[,4]

#import sample IDs
ids <- read.table("stq_dates_yeardec", sep="\t", header=TRUE)[,1]


##############################
#Reformat data
##############################

#follow RDA tutorial by Dr. Brenna Forester
#file:///C:/Users/marc_/Downloads/RDA_landscape_genomics.html

#make sample IDs the rownames of genotype matrix
row.names(imputed_geno) <- ids

#standardize environmental data
env_all_scaled <- scale(env_all, center=TRUE, scale=TRUE)

#combine different data
variables <- data.frame(ids, PC1=pc_scores[,1], PC2=pc_scores[,2], longitude=latlongs[,1], latitude=latlongs[,2], year, env_all_scaled)
colnames(variables) <- c("ID", "PC1", "PC2", "Longitude", "Latitude", "Collection_Year",
                       "Bio1",
                       "Bio2",
                       "Bio3",
                       "Bio4",
                       "Bio12",
                       "Bio15",
                       "Elev",
                       "Mod",
                       "Def",
                       "Mol",
                       "Ngl",
                       "Nef",
                       "Raf",
                       "Swl",
                       "Scr",
                       "Wef",
                       "Oth",
                       "Gen_DFTD",
                       "Gen05_Devil",
                       "Gen10_Devil",
                       "Gen15_Devil")


##############################
#Run RDA
##############################

#run initial RDA, conditioning on pop structure and collection year
##vif.cca identifies the variable with the largest inflation factor, which is subsequently removed and the process repeated
rda_env <- rda(imputed_geno ~ Bio1 + Bio2 + Bio3 + Bio4 + Bio12 + Bio15 + Elev + Mod + Def + Mol + Ngl + Nef + Raf + Swl + Scr + Wef + Oth + Gen_DFTD + Gen05_Devil + Gen10_Devil + Gen15_Devil + Condition(PC1 + PC2 + Collection_Year),  variables)
which.max(vif.cca(rda_env))
rda_env2 <- rda(imputed_geno ~ Bio1 + Bio3 + Bio4 + Bio12 + Bio15 + Elev + Mod + Def + Mol + Ngl + Nef + Raf + Swl + Scr + Wef + Gen_DFTD + Gen05_Devil + Gen10_Devil + Gen15_Devil + Condition(PC1 + PC2 + Collection_Year),  variables)
which.max(vif.cca(rda_env2))
rda_env3 <- rda(imputed_geno ~ Bio1 + Bio3 + Bio4 + Bio12 + Bio15 + Mod + Def + Mol + Ngl + Nef + Raf + Swl + Scr + Wef + Gen_DFTD + Gen05_Devil + Gen10_Devil + Gen15_Devil + Condition(PC1 + PC2 + Collection_Year),  variables)
which.max(vif.cca(rda_env3))
rda_env4 <- rda(imputed_geno ~ Bio1 + Bio3 + Bio4 + Bio12 + Bio15 + Mod + Def + Mol + Ngl + Nef + Raf + Swl + Scr + Wef + Gen_DFTD + Gen05_Devil + Gen10_Devil + Condition(PC1 + PC2 + Collection_Year),  variables)
which.max(vif.cca(rda_env4))
rda_env5 <- rda(imputed_geno ~ Bio1 + Bio3 + Bio4 + Bio12 + Bio15 + Mod + Def + Mol + Ngl + Nef + Raf + Swl + Scr + Wef + Gen_DFTD + Gen10_Devil + Condition(PC1 + PC2 + Collection_Year),  variables)
which.max(vif.cca(rda_env5))
rda_env6 <- rda(imputed_geno ~ Bio1 + Bio3 + Bio4 + Bio12 + Bio15 + Def + Mol + Ngl + Nef + Raf + Swl + Scr + Wef + Gen_DFTD + Gen10_Devil + Condition(PC1 + PC2 + Collection_Year),  variables)
which.max(vif.cca(rda_env6))

#use the covariate set of rda_env6 as rda_env; lets the rest of code stay the same
rda_env <- rda(imputed_geno ~ Bio1 + Bio3 + Bio4 + Bio12 + Bio15 + Def + Mol + Ngl + Nef + Raf + Swl + Scr + Wef + Gen_DFTD + Gen10_Devil + Condition(PC1 + PC2 + Collection_Year),  variables)

#if not re-running above analyses, can load in rdadapt output for K=4 below:
rda_env <- readRDS("rda_env.gz")


#screeplot suggests we should keep first four or five axes
screeplot(rda_env, main="Eigenvalues of constrained axes")
rda_axis_names <- factor(names(rda_env$CCA$eig)[1:15], levels=c("RDA1", "RDA2", "RDA3", "RDA4", "RDA5", "RDA6", "RDA7", "RDA8", "RDA9", "RDA10", "RDA11", "RDA12", "RDA13", "RDA14", "RDA15"))
rda_eigval <- data.frame(Axis = rda_axis_names, Eigenvalue = rda_env$CCA$eig[1:15])
p_rda_screeplot <- ggplot()+
  geom_bar(data=rda_eigval, aes(x=Axis, y=Eigenvalue, fill=Axis), stat="identity")+
  ylim(0,10)+
  scale_fill_manual(values=c(rep("gray40", 5), rep("gray70", 10)), guide="none")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45, vjust=0.5))
p_rda_screeplot

#Function to conduct a RDA based genome scan
rdadapt <- function(rda,K)
{
  zscores <- rda$CCA$v[,1:as.numeric(K)]
  resscale <- apply(zscores, 2, scale)
  resmaha <- covRob(resscale, distance = TRUE, na.action= na.omit, estim="pairwiseGK")$dist
  lambda <- median(resmaha)/qchisq(0.5,df=K)
  reschi2test <- pchisq(resmaha/lambda,K,lower.tail=FALSE)
  qval <- qvalue(reschi2test)
  q.values_rdadapt<-qval$qvalues
  return(data.frame(p.values=reschi2test, q.values=q.values_rdadapt))
}

#run rdadapt function with K=5
rdadapt_env <- rdadapt(rda_env, 5)

#if not re-running above analyses, can load in rdadapt output for K=4 below:
rdadapt_env <- readRDS("rdadapt_env_k5.gz")


#identify loci that are below qvalue threshold
outliers_qval <- data.frame(Loci = locus_ids[which(rdadapt_env$q.values<0.01)], q.value = rdadapt_env$q.values[which(rdadapt_env$q.values<0.01)])
outliers_qval_index <- data.frame(locus_index = which(rdadapt_env$q.values<0.01) ,outliers_qval)

#isolate just SNP index
rda_qval_outliers <- which(rdadapt_env$q.values<0.01)


##############################
#identify environmental factor with greatest correlation with each candidate SNP
##############################

#get axis loadings for all SNPs
load_rda_all <- data.frame(snp_index=1:3431, scores(rda_env, choices=c(1:5), display="species"))

#subset out axis loadings for outlier SNPs
load_rda_outliers <- load_rda_all[rda_qval_outliers,]

#subset env variables df to retain only those used in the final RDA function call
retained_variables <- subset.data.frame(variables, select=c("Bio1", "Bio3", "Bio4", "Bio12", "Bio15", "Def", "Mol", "Ngl", "Nef", "Raf", "Swl", "Scr", "Wef", "Gen_DFTD", "Gen10_Devil"))

###
#Pearson's correlation

#find Pearson's correlation of each outlier SNP with the fifteen environmental predictors
pearcor <- matrix(nrow=nrow(load_rda_outliers), ncol=15)  # 15 columns for 15 predictors
colnames(pearcor) <- c("Bio1", "Bio3", "Bio4", "Bio12", "Bio15", "Def", "Mol", "Ngl", "Nef", "Raf", "Swl", "Scr", "Wef", "Gen_DFTD", "Gen10_Devil")

#loop through outlier SNPs and estimate Pearson's correlation coefficient
#modified from https://popgen.nescent.org/2018-03-27_RDA_GEA.html
for (i in 1:nrow(load_rda_outliers)){
  snp_index <- load_rda_outliers$snp_index[i]
  snp_geno <- imputed_geno[,snp_index]
  pearcor[i,] <- apply(X=retained_variables, MARGIN=2, FUN=function(x) cor(x, snp_geno))
}


###
#Binomial correlation
binomcor <- matrix(nrow=nrow(load_rda_outliers), ncol=15)
colnames(binomcor) <- c("Bio1", "Bio3", "Bio4", "Bio12", "Bio15", "Def", "Mol", "Ngl", "Nef", "Raf", "Swl", "Scr", "Wef", "Gen_DFTD", "Gen10_Devil")

#define function for carrying out binomial regression and calculating McFadden's pseudo-R2 from binomial GLM (1 - deviance/null.deviance)
binom_r2_calc <- function(x, geno){
  snp_binom <- glm(snp_binom_res ~ x, family=binomial(link="logit"))
  r2 <- 1-snp_binom$deviance/snp_binom$null.deviance
}

#loop through outlier SNPs and estimate pseudo R-squared from binomial models
for (i in 1:nrow(load_rda_outliers)){
  snp_index <- load_rda_outliers$snp_index[i]
  snp_geno <- imputed_geno[,snp_index]
  snp_binom_res <- as.matrix(data.frame(snp_geno, 2)) #create response matrix for binomial glm
  binomcor[i,] <- apply(X=retained_variables, MARGIN=2, FUN=function(x, geno) binom_r2_calc(x, geno=snp_binom_res))
}

###
#identify env factor with maximum correlation

#pearson correlation
snp_max_cor_pearson <- data.frame(snp_index=load_rda_outliers[,1], env_var=rep(NA, nrow(load_rda_outliers)), cor_pearson=rep(NA, nrow(load_rda_outliers)))
for (i in 1:nrow(load_rda_outliers)){
  snp_cor <- pearcor[i,]
  snp_max_cor_pearson$env_var[i] <- colnames(pearcor)[which.max(abs(snp_cor))]
  snp_max_cor_pearson$cor_pearson[i] <- max(abs(snp_cor))
}

#binomial pseudo-r2
snp_max_cor_binom <- data.frame(snp_index=load_rda_outliers[,1], env_var=rep(NA, nrow(load_rda_outliers)), pseudo_r2=rep(NA, nrow(load_rda_outliers)))
for (i in 1:nrow(load_rda_outliers)){
  snp_cor <- binomcor[i,]
  snp_max_cor_binom$env_var[i] <- colnames(binomcor)[which.max(abs(snp_cor))]
  snp_max_cor_binom$pseudo_r2[i] <- max(abs(snp_cor))
}

#check consistency of env factor with best fit between pearson's r and binomial pseudo-R2
##calculates the percent of SNPs for which Pearson coefficient and Binomial pseudo-R2 identified the same environmental factor as showing the strongest association
sum(snp_max_cor_pearson$env_var==snp_max_cor_binom$env_var)/nrow(snp_max_cor_pearson)

#export outlier-env correlation data
rda_outlier_loadings_pearson <- data.frame(load_rda_outliers, pearcor)
rda_outlier_loadings_binom <- data.frame(load_rda_outliers, binomcor)

#write.table(x=rda_outlier_loadings_pearson, file="rda_outlier_loadings_pearson", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
#write.table(x=rda_outlier_loadings_binom, file="rda_outlier_loadings_binom", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
#write.table(x=snp_max_cor_pearson, file="rda_outliers_max_cor_pearson", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
#write.table(x=snp_max_cor_binom, file="rda_outliers_max_cor_binom", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
