
##############################
#Load in significant SNPs from LFMM and RDA
##############################

#load in LFMM, keeping only unique SNPs
lfmm_sig <- unlist(read.table("lfmm_sig_index_q001", header=TRUE, sep="\t"))
lfmm_sig <- unique(lfmm_sig[which(!is.na(lfmm_sig))])

#load in RDA
rda_sig <- read.table("rda_outliers_q001", header=TRUE, sep="\t")[,1]


#intersect LFMM and RDA results
rda_qval_outliers<-which(rdadapt_env$q.values<0.01)
lfmm_sig_all<-unlist(lfmm_qval_sig_tab)
lfmm_sig_all<-lfmm_sig_all[which(!is.na(lfmm_sig_all))]

##############################
#Intersect significant SNPs from LFMM and RDA
##############################

overlap_sig <- intersect(lfmm_sig, rda_sig)
write.table(data.frame(rda_lfmm_overlap=overlap_sig), "lfmm_rda_overlap_index_sig_q001", row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
