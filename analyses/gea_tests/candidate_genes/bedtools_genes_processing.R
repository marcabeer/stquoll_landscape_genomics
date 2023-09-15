
library(vcfR)
library(tidyr)
library(tidyverse)

##############################
#Load in locus info
##############################

#load locus info for all SNPs
stq_locus_info <- read.table("locus_info", header=TRUE, sep="\t")

#put positional info into BED format for all SNPs
stq_chrom_pos <- stq_locus_info[,c(1,2)]
stq_chrom_pos_comp <- paste(stq_chrom_pos[,1], (as.numeric(stq_chrom_pos[,2])-1), sep="_")

stq_locus_info_mapped_end <- as.numeric(stq_chrom_pos[,2])
stq_locus_info_mapped_start <- stq_locus_info_mapped_end-1

stq_all_snps_bedtools <- data.frame(cbind(stq_locus_info[,1], stq_locus_info_mapped_start, stq_locus_info_mapped_end))

#output df of all SNPs in BED format
##used as input for bedtools closest
write.table(x=stq_all_snps_bedtools, file="stq_all_snps_bedtools", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)


##############################
#After running BEDTools "closest" command with the above table as input, process the output
##############################

#load bedtools closest output file
bt_out_raw <- read.table("closest_genes", sep="\t")

#split column containing genomic feature information
bt_out <- bt_out_raw %>% separate(col=V13,
                    sep=";",
                    into=c("ID", "Dbxref", "Name", "gbkey", "gene", "gene_biotype", "note1", "note2"),
                    fill="right")


#add column containing a composite SNP ID (chr_pos)
bt_out$comp_id <- paste(bt_out[,1], bt_out[,2], sep="_")

#add column containing SNP index
bt_out$snp_index <- NA
for (i in 1:nrow(bt_out)){
  bt_out$snp_index[i] <-  match(bt_out$comp_id[i], stq_chrom_pos_comp)
}

#output file of all unique gene symbols, which can be used for other things such as GO enrichment
write.table(unique(bt_out$gene), "background_genes", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)


##############################
#Load in indices of SNPs with significant GEAs
##############################

#RDA significant outliers
rda_sig <- read.table("rda_outliers_max_cor_binom", sep="\t", header=TRUE)

#LFMM significant associations
lfmm_sig <- read.table("lfmm_sig_index_q001", sep="\t", header=TRUE)

#RDA-LFMM overlapping SNPs
shared_sig <- read.table("lfmm_rda_overlap_index_sig_q001", sep="\t", header=TRUE)


##############################
#find ids of genes near SNPs with significant GEAs determined by RDA and LFMM
##############################

#create vectors of RDA and LFMM significant SNP indices
lfmm_all <- sort(na.omit(unlist(lfmm_sig)), decreasing=FALSE)
rda_all <- rda_sig[,1]

#create vector of SNPs found to be significant by either RDA or LFMM
both_all <- sort(unique(c(lfmm_all, rda_all)), decreasing=FALSE)

#for lfmm-identified SNPs, obtain environmental factor of association
lfmm_table <- data.frame(lfmm_bio1=rep(NA, nrow(bt_out)),
                       lfmm_bio3=rep(NA, nrow(bt_out)),
                       lfmm_bio4=rep(NA, nrow(bt_out)),
                       lfmm_bio12=rep(NA, nrow(bt_out)),
                       lfmm_bio15=rep(NA, nrow(bt_out)),
                       lfmm_def=rep(NA, nrow(bt_out)),
                       lfmm_mol=rep(NA, nrow(bt_out)),
                       lfmm_ngl=rep(NA, nrow(bt_out)),
                       lfmm_nef=rep(NA, nrow(bt_out)),
                       lfmm_raf=rep(NA, nrow(bt_out)),
                       lfmm_swl=rep(NA, nrow(bt_out)),
                       lfmm_scr=rep(NA, nrow(bt_out)),
                       lfmm_wef=rep(NA, nrow(bt_out)),
                       lfmm_dftd=rep(NA, nrow(bt_out)),
                       lfmm_devil=rep(NA, nrow(bt_out))
)

for (i in 1:length(lfmm_sig_split)){
  if (length(lfmm_sig_split[[i]])>0){
    for (j in 1:length(lfmm_sig_split[[i]])){
      lfmm_table[which(bt_out$snp_index==lfmm_sig_split[[i]][j]),i] <- "x"
    }
  }
}


#for rda-identified snps, obtain environmental factor of association
rda_outliers_env <- read.table("rda_outliers_max_cor_binom", sep="\t", header=TRUE)[,-3]

rda_outliers_split <- rda_outliers_env %>%
  split(.$env_var)

for (i in 1:length(rda_outliers_split)){
  rda_outliers_split[[i]] <- as.integer(unlist(rda_outliers_split[[i]][,1]))
}

rda_table <- data.frame(rda_bio1=rep(NA, nrow(bt_out)),
                      rda_bio12=rep(NA, nrow(bt_out)),
                      rda_bio15=rep(NA, nrow(bt_out)),
                      rda_bio3=rep(NA, nrow(bt_out)),
                      rda_bio4=rep(NA, nrow(bt_out)),
                      rda_dftd=rep(NA, nrow(bt_out)),
                      rda_devil=rep(NA, nrow(bt_out)),
                      rda_mol=rep(NA, nrow(bt_out)),
                      rda_nef=rep(NA, nrow(bt_out)),
                      rda_ngl=rep(NA, nrow(bt_out)),
                      rda_scr=rep(NA, nrow(bt_out)),
                      rda_swl=rep(NA, nrow(bt_out)),
                      rda_wef=rep(NA, nrow(bt_out))
)

for (i in 1:length(rda_outliers_split)){
  if (length(rda_outliers_split[[i]])>0){
    for (j in 1:length(rda_outliers_split[[i]])){
      rda_table[which(bt_out$snp_index==rda_outliers_split[[i]][j]),i] <- "x"
    }
  }
}

rda_table <- data.frame(rda_table[,c(1,4,5,2,3)], rda_def=rep(NA, nrow(bt_out)), rda_mol=rda_table$rda_mol, rda_ngl=rda_table$rda_ngl, rda_nef=rda_table$rda_nef, rda_raf=rep(NA, nrow(bt_out)), rda_swl=rda_table$rda_swl, rda_scr=rda_table$rda_scr, rda_wef=rda_table$rda_wef, rda_table[,c(6,7)])

#bind lfmm_table and rda_table to bedtools output
bt_out_bound <- data.frame(bt_out, both_geas=NA, lfmm_table, rda_table)

#isolate bed entries associated with significant SNPs
bt_sig_genes <- bt_out_bound[bt_out_bound$snp_index %in% both_all,]

#make sure all entries are marked as signficant for at least one test/env factor
##maximum number of NAs (one test/env factor is significant) is 29
isna_sig_genes <- list()
for (i in 1:nrow(bt_sig_genes)){
  isna_sig_genes[[i]] <- sum(is.na(bt_sig_genes[i,25:54]))
  
}
max(unlist(isna_sig_genes))

#populate both_geas column for SNPs identified by both RDA and LFMM
for (i in 1:nrow(bt_sig_genes)){
  if (sum(is.na(bt_sig_genes[i, 25:39]))<15 & sum(is.na(bt_sig_genes[i, 40:54]))<15){
    bt_sig_genes$both_geas[i] <- "x"
  }
}


##############################
#Format table for export
##############################

#isolate columns of interest
bt_sig_genes_trimmed <- bt_sig_genes[,-c(4,7,8,9,10,11,12,19,20)]

#reorder columns
bt_sig_genes_trimmed <- bt_sig_genes_trimmed[,c(13, 14, 1:12, 15:45)]
colnames(bt_sig_genes_trimmed)[1] <- "snp_id"
colnames(bt_sig_genes_trimmed)[2] <- "snp_index"
colnames(bt_sig_genes_trimmed)[3] <- "chr"
colnames(bt_sig_genes_trimmed)[4] <- "snp_start"
colnames(bt_sig_genes_trimmed)[5] <- "snp_end"
colnames(bt_sig_genes_trimmed)[6] <- "feature_start"
colnames(bt_sig_genes_trimmed)[7] <- "feature_end"
colnames(bt_sig_genes_trimmed)[14] <- "distance"

write.table(bt_sig_genes_trimmed, file="sig_genes_table", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
