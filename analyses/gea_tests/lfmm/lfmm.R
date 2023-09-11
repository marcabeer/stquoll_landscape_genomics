library(pegas)
library(ggplot2)
library(raster)
library(rgdal)
library(lfmm)
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
library(ggrepel)

##############################
#Load in genetic data and metadata
##############################

#load in imputed genotype matrix (individual genotypes, not population allele frequencies)
imputed_geno <- read.lfmm("stq.recode.lfmm_imputed.lfmm")

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
#reformat data
##############################

#make sample IDs the rownames of genotype matrix
row.names(imputed_geno)<-ids

#standardize environmental data
env_all_scaled<-scale(env_all, center=TRUE, scale=TRUE)
colnames(env_all_scaled)<-c("Bio1",
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

env_pruned <- env_all_scaled[,c(1,3,4,5,6,9,10,11,12,13,14,15,16,18,20)]


##############################
#run LFMM
##############################

#running lfmm
lfmm_res<-lfmm_ridge(Y=imputed_geno, X=env_pruned, K=3)
lfmm_res_test <- lfmm_test(Y = imputed_geno, X = env_pruned, lfmm = lfmm_res, calibrate = "median+MAD")

#save LFMM results to file
#saveRDS(lfmm_res, file="lfmm_res") #can be loaded with readRDS()
#saveRDS(lfmm_res_test, file="lfmm_res_test") #can be loaded with readRDS()

#if not re-running analysis, read in our results below:
lfmm_res <- readRDS("lfmm_res.gz")
lfmm_res_test <- readRDS("lfmm_res_test.gz")

##############################
#work through LFMM outputs
##############################

#obtain p-values
lfmm_pval <- lfmm_res_test$pvalue

#convert p-values to q-values
lfmm_qval <- qvalue(lfmm_pval)$qvalues

#identify SNPs significantly associated with each env factor at q<0.01
lfmm_qval_sig <- list()
for (i in 1:ncol(lfmm_pval)){
  lfmm_qval_sig[[i]] <- which(lfmm_qval[,i]<0.01)
}
names(lfmm_qval_sig) <- colnames(lfmm_qval)

#reformat sig SNPs into table
lfmm_qval_sig_tab <- data.frame(lapply(lfmm_qval_sig, 'length<-', max(lengths(lfmm_qval_sig))))

#calculate total number of unique significant SNPs
sum(!is.na(unique(unlist(lfmm_qval_sig_tab))))

#save significant associations (i.e., the SNP indices) for each environmental variable
#write.table(lfmm_qval_sig_tab, "lfmm_sig_index_q001", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
