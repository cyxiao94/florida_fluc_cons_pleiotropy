## ----setup, include=FALSE-----------------------------------------------------
#knitr::opts_chunk$set(echo = TRUE)

options(warn = -1)
## -----------------------------------------------------------------------------
library(ggplot2)
library(data.table)
library(poolSeq)
#library(lme4)




## -----------------------------------------------------------------------------
#this follows the order of the output allele frequency and coverage from poolSeq package
sample_list<-NULL
for (r in 1:10){
  #sample_list<-c(sample_list,paste0("anc.0.",r))
  #for (pop in c("fluc","cons")){
    for (gen in c(0,10,20,30,40,50,60)){
      sample_list<-c(sample_list,paste0("F",gen,".R",r))
    }
  #}
}

#c("pop","gen","rep")
split_list<- tstrsplit(sample_list, split="\\.")
gen <- as.numeric(sub("F","",split_list[[1]]))
rep <- as.numeric(sub("R","",split_list[[2]]))


## -----------------------------------------------------------------------------
list_gen<-c(0,10,20,30,40,50,60) 
list_repl<-1:10
freq<-NULL
cov<-NULL

sync<-read.sync(file = "./01genomic_data/Dsim_F0-F60_Q20_polymorphic.sync.gz",
                    gen = rep(list_gen,each=length(list_repl)),
                    repl = rep(list_repl,length(list_gen)), polarization = "rising")
freq<-af(sync,  
             gen = rep(c(0,10,20,30,40,50,60),each=length(list_repl)), #remove F160
             repl = rep(list_repl,length(list_gen)))
#colnames(freq)<-paste0(sample_list,".freq")

cov<-coverage(sync,  
              gen = rep(c(0,10,20,30,40,50,60),each=length(list_repl)),
              repl = rep(list_repl,length(list_gen)))
#colnames(cov)<-paste0(sample_list,".cov")
rm(sync)  


## -----------------------------------------------------------------------------
#pos <- splitLocusID(rownames(freq))
var <- apply(freq, 1, function(x) var(na.omit(x)))
ind<-which(var == 0)
#if (length(ind) !=0){
  freq_filtered <- freq[-ind,]
  cov_filtered <- cov[-ind,]
#}
#rm(freq, cov)

#save data

##HAFs
###afmat
  afmat <- freq
  colnames(afmat) <- sub("\\.freq$", "", colnames(afmat))
###samps
  samps <- data.table(sampID = colnames(afmat))
  samps[ , c("generation","replicate") := tstrsplit(sampID,"\\.")]
  samps[, tpt := as.numeric(sub("F","",generation))/10+1]
  samps[, cage := as.numeric(sub("R","",replicate))]
  samps <- samps[ , c("sampID","tpt","cage","generation","replicate")]
###sites
  sites <- as.data.frame(splitLocusID(rownames(afmat)))
  colnames(sites) <- c("chrom","pos")
save(afmat, samps, sites, file = "01genomic_data/Florida_1stwave_HAF.RData")

##cov
covmat <- cov
colnames(covmat) <- sub("\\.cov$", "", colnames(afmat))
saveRDS(covmat, file = "01genomic_data/Florida_1stwave_cov.RDS")

rm(afmat, samps, sites, covmat)
