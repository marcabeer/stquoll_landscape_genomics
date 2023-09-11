library(ggplot2)
library(cowplot)
library(tidyr)
library(plotrix)
library(adegenet)
options(digits=15)


##############################
#Load in marginal likelihoods of replicate FastStructure runs
##############################

#read in list of raw marginal likelihoods
##first column has log file name, which contains replicate number followed by k number (r#.#)
mlikelihoods_raw <- read.table(file="mlikelihood_format.txt", sep="\t", header=FALSE)

#split filename and rename columns
file_name_split <- separate(mlikelihoods_raw, col=1, into=c("test", "rep", "k", "log"))

#save only columns containing replicate number, k value, and marginal likelihood
mlikelihoods_v1<-file_name_split[,c(2,3,5)]
colnames(mlikelihoods_v1)<-c("rep", "k", "mlikelihood")

#sort by value of K
mlikelihoods_v1_sorted <- mlikelihoods_v1[order(as.numeric(as.character(mlikelihoods_v1[,2]))),]


##############################
#Summarize replicate FastStructure runs for each value of K
##############################

#define function to summarize replicates of each value of K (mean and std.error of marginal likelihood)
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = std.error(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

#get mean and standard error of marginal likelihood for each value of k
k_summ <- data_summary(data=mlikelihoods_v1_sorted, varname="mlikelihood", 
                    groupnames="k")

#re-sort by K value
k_summ_sorted <- k_summ[order(as.numeric(as.character(k_summ[,1]))),]


##############################
#Plot marginal likelihood versus K
##############################

#plot marginal likelihood versus k
level_order<-c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20")

faststructure_marginal_likelihoods <- ggplot(k_summ_sorted, aes(x=factor(k, level=level_order), y=mlikelihood, group=1)) + 
  geom_point(size=1)+
  geom_line(linewidth=0.5)+
  geom_errorbar(aes(ymin=mlikelihood-sd, ymax=mlikelihood+sd), width=.2)+
  ylim(-0.75,-0.675)+
  xlab("Number of clusters")+
  ylab("Marginal likelihood")+
  theme(axis.title.x = element_text(margin=margin(t=10),color="black", size=14),
        axis.title.y = element_text(margin=margin(r=15), color="black", size=14))+
  theme_bw()+
  theme(axis.text.x = element_text(size=10, angle=0), axis.text.y = element_text(size=10, angle=0))+
  theme(legend.position="none")

faststructure_marginal_likelihoods

