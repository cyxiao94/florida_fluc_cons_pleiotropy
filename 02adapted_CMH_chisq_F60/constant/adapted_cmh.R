require(poolSeq)
require(data.table)
require(dplyr)
require(tidyr)
require(IHW)
require(ACER)
require(fdrtool)
library(gridExtra)
library(grid)
require(mclust)



##intial parameters-------
input_gen=c(0,60)
input_rep=5
#fet_order=c(2*input_rep+1,2*input_rep+2,seq(from=1, to=2*input_rep-1, by=2),seq(from=2, to=2*input_rep, by=2))
fet_order<-c("chr","pos",paste("r",1:input_rep,"_p",sep=""), paste("r",1:input_rep,"_adj_p",sep=""))
chromosomes=c("2L","2R","3L","3R","4","X")
setwd("~/my_project/fluctuation_constant/20final_dataset/01adapted_CMH/constant/")

#1. remove 0 variance site--------------
list_no_var <- NULL
  for (chr in chromosomes){
    sync1 <- read.sync(file = paste("00",chr, ".filtered.cons.q20.no_0cov.sync", sep = ""), 
                      gen = rep(c(0,10,20,30,40,50,60), each = input_rep), repl = rep(1:input_rep, 7), polarization = "rising")
    sync_todo<-tbl_df(fread(file = paste("00",chr, ".filtered.cons.q20.no_0cov.sync", sep = ""), h = F, stringsAsFactors = F, sep = "\t"))
    freqf <- af(sync, gen = rep(input_gen, each = input_rep), repl = rep(1:input_rep, 2)); col_freq <- colnames(freqf)
    pos <-  splitLocusID(rownames(freqf))
    var <- apply(freqf, 1, function(x) var(na.omit(x)))
    ind <- which(var == 0)
    if (length(ind) != 0) {pos_keep<-pos$pos[-ind]
    sync_todo<-subset(sync_todo, V2 %in% pos_keep); rm(pos_keep)}
    write.table(subset(sync_todo,select = paste("V",c(1:(input_rep+3),(6*input_rep+4):(7*input_rep+3)), sep="")), paste("01no_0var.F0_F60.",chr, ".filtered.cons.q20.no_0cov.sync", sep = ""),
                col.names = F, row.names = F, quote = F, sep= "\t")
    rm(sync);rm(sync_todo);rm(freqf);rm(pos);rm(ind);rm(var);rm(chr)
  }

#saveRDS(list_no_var, "var_freq_per_mapper_F0_F70_all_chr_subsnps.RData")

#2.remove site that contain NA
list_no_NA <- NULL
for (chr in chromosomes){
  sync <- read.sync(file = paste("01no_0var.F0_F60.",chr,".filtered.cons.q20.no_0cov.sync", sep = ""),
                    gen = rep(c(0,60), each = input_rep), repl = rep(1:input_rep, 2), polarization = "rising")
  sync_todo<-tbl_df(fread(file = paste("01no_0var.F0_F60.",chr,".filtered.cons.q20.no_0cov.sync", sep = ""), h = F, stringsAsFactors = F, sep = "\t"))
  freqf <- af(sync, gen = rep(input_gen, each = input_rep), repl = rep(1:input_rep, 2)); col_freq <- colnames(freqf)
  pos <-  splitLocusID(rownames(freqf))
  ind<-which(rowSums(is.na(freqf))>0)
  if (length(ind) != 0) {pos_keep<-pos$pos[-ind]
  sync_todo<-subset(sync_todo, V2 %in% pos_keep); rm(pos_keep)}
  write.table(sync_todo, paste("02noNA.no_0var.F0_F60.",chr, ".filtered.cons.q20.no_0cov.sync", sep = ""),
              col.names = F, row.names = F, quote = F, sep= "\t")
  rm(sync);rm(sync_todo);rm(freqf);rm(pos);rm(ind);rm(chr)
}

#2. count coverage for genome wide------------
#cov_all <- NULL
#for (chr in chromosomes){
#  sync <- read.sync(file = paste("01no_0var.F0_F60.", chr, ".filtered.cons.q20.no_0cov.sync", sep = ""),
#                           gen = rep(input_gen, each = input_rep), repl = rep(1:input_rep, 2), polarization = "rising")
#  covf <- poolSeq::coverage(sync, gen = rep(input_gen, each = input_rep), repl = rep(1:input_rep, 2)); 
#  cov_all <- rbind(cov_all, covf);
#  rm(sync);rm(covf)}
#summary_count<-cbind(data.frame(repl=rep(c(1:input_rep), 2), gen = rep(input_gen,each=input_rep)),
#                   as.data.frame(t(apply(subset(cov_all, select = paste("F", rep(input_gen,each=input_rep), ".R", rep(1:input_rep), ".cov", sep = "")), 2, function(x) 
#                   c(mean(x, na.rm = T), sd(x, na.rm = T), sd(x, na.rm = T)/mean(x, na.rm = T), median(x, na.rm = T), 
#                     quantile(x, c(0.001, 0.01, 0.02, 0.1, 0.25, 0.75, 0.98, 0.99, 0.999), na.rm = T))))))
#colnames(summary_count)<-c("repl", "gen", "mean", "sd", "sd_mean", "median", "0.1%", "1%", "2%", "10%", "25%", "75%", "98%", "99%", "99.9%")
#write.table(summary_count,"summary_count_F0_F60_nodup_q20.txt",sep="\t", quote = F, col.names = T, row.names = T)
#rm(cov_all);rm(summary_count)
###2.1remove extreme SNPs-------------
#summary_count<-read.table("summary_count_F0_F60_nodup_q20.txt", sep="\t", h=T)
#base
#thres_down_F0<-summary_count[which(summary_count$gen==0),"X1."]
#thres_up_F0<-summary_count[which(summary_count$gen==0),"X99."]
#F60
#thres_down_F60<-summary_count[which(summary_count$gen==60),"X1."]
#thres_up_F60<-summary_count[which(summary_count$gen==60),"X99."]
#for (chr in chromosomes){
#  sync_modify<-fread(file = paste("01no_0var.F0_F60.", chr, ".filtered.cons.q20.no_0cov.sync", sep = ""), h = F, stringsAsFactors = F, sep = "\t")
#  sync<-read.sync(file = paste("01no_0var.F0_F60.", chr, ".filtered.cons.q20.no_0cov.sync", sep = ""),
#                  gen = rep(input_gen, each = input_rep), repl = rep(1:input_rep, 2), polarization = "rising")
#  covf<-as.data.frame(poolSeq::coverage(sync, gen =rep(input_gen, each = input_rep), repl = rep(1:input_rep, 2)))
#  covf<-cbind(covf,splitLocusID(rownames(covf)))
#  discard<-which(covf$F0.R1.cov<thres_down_F0[1] | covf$F0.R2.cov<thres_down_F0[2] | covf$F0.R3.cov<thres_down_F0[3] | covf$F0.R4.cov<thres_down_F0[4] | covf$F0.R5.cov<thres_down_F0[5] | #covf$F0.R6.cov<thres_down_F0[6] #base_down
#                 covf$F60.R1.cov<thres_down_F60[1] | covf$F60.R2.cov<thres_down_F60[2] | covf$F60.R3.cov<thres_down_F60[3] | covf$F60.R4.cov<thres_down_F60[4] | covf$F60.R5.cov<thres_down_F60[5] | #covf$F60.R6.cov<thres_down_F60[6] | #F60_down
#                 covf$F0.R1.cov>thres_up_F0[1] | covf$F0.R2.cov>thres_up_F0[2] | covf$F0.R3.cov>thres_up_F0[3] | covf$F0.R4.cov>thres_up_F0[4] | covf$F0.R5.cov>thres_up_F0[5] | #covf$F0.R6.cov>thres_up_F0[6] | #base_up
#                 covf$F60.R1.cov>thres_up_F60[1] | covf$F60.R2.cov>thres_up_F60[2] | covf$F60.R3.cov>thres_up_F60[3] | covf$F60.R4.cov>thres_up_F60[4] | covf$F60.R5.cov>thres_up_F60[5] )#| covf$F60.R6.cov>thres_up_F60[6]) #F60_up
#  site_to_keep<-covf$pos[-discard]
#  sync_modify<-subset(sync_modify, V2 %in% site_to_keep)
#  write.table(sync_modify, paste("02coverage.no_0var.F0_F60.", chr, ".filtered.cons.q20.no_0cov.sync", sep = ""),
#              col.names = F, row.names = F, quote = F, sep= "\t")
#  rm(covf);rm(sync);rm(sync_modify);rm(discard);rm(site_to_keep)}
#rm(chr);rm(summary_count);rm(thres_down_F0);rm(thres_down_F60);rm(thres_up_F0);rm(thres_up_F60)


#3. estimate Ne-------------
###3.1 autosome---------
med_Ne_auto=NULL
for (i in 1:input_rep){
  ne=NULL;af=NULL;cov=NULL
  for (chr in c("2L","2R","3L","3R","4")){
    sync<-read.sync(file=paste("02noNA.no_0var.F0_F60.",chr,".filtered.cons.q20.no_0cov.sync",sep=""), gen = rep(input_gen, each = input_rep), repl = rep(1:input_rep, 2), polarization = "rising")
    af=rbind(af, as.data.frame(af(sync, gen = c(0,60), repl = i)))
    cov=rbind(cov, as.data.frame(poolSeq::coverage(sync, gen = c(0,60), repl = i)))}
    colnames(af)<-c("F0","F60"); colnames(cov)<-c("F0","F60")
    NB_trials <- 100
    off <- 1
    nb_obs <- dim(af)[1]
    #ne <- NULL
    for(j in 1:NB_trials){
      set.seed(off)
      ind <- sample(x = 1:nb_obs, size = 1000)
      ne <- rbind(ne, data.frame(trial = j,
                                 ne = estimateNe(p0 = af[ind,"F0"], pt = af[ind,"F60"], cov0 = cov[ind,"F0"], covt = cov[ind,"F60"],
                                                 t = 60, ploidy=2, truncAF=0.05, method=c("P.planI"), poolSize=rep(600, times=2), Ncensus=1250)))
      off <- off+1
    }
  Ne <- median(ne$ne)
  ind <- which(is.na(ne$ne) | ne$ne<0)
  if(length(ind)>0){ne <- ne[-ind, ]}
  med_Ne_auto<-c(med_Ne_auto,round(median(ne$ne)))
}
rm(af);rm(cov);rm(ne);rm(sync);rm(chr);rm(i);rm(ind);rm(j);rm(nb_obs);rm(NB_trials);rm(Ne);rm(off)
###3.2 X chromosome-----------------
med_Ne_X=NULL
for (i in 1:input_rep){
  ne=NULL
  sync<-read.sync(file="02noNA.no_0var.F0_F60.X.filtered.cons.q20.no_0cov.sync", gen = rep(input_gen, each = input_rep), repl = rep(1:input_rep, 2), polarization = "rising")
  af=as.data.frame(af(sync, gen = c(0,60), repl = i));colnames(af)<-c("F0","F60")
  cov=as.data.frame(poolSeq::coverage(sync, gen = c(0,60), repl = i)); colnames(cov)<-c("F0","F60")
  NB_trials <- 100
  off <- 1
  nb_obs <- dim(af)[1]
  #ne <- NULL
  for(j in 1:NB_trials){
    set.seed(off)
    ind <- sample(x = 1:nb_obs, size = 1000)
    ne <- rbind(ne, data.frame(trial = j,
                               ne = estimateNe(p0 = af[ind,"F0"], pt = af[ind,"F60"], cov0 = cov[ind,"F0"], covt = cov[ind,"F60"],
                                               t = 60, ploidy=2, truncAF=0.05, method=c("P.planI"), poolSize=rep(600, times=2), Ncensus=1250)))
    off <- off+1
    }
  Ne <- median(ne$ne)
  ind <- which(is.na(ne$ne) | ne$ne<0)
  if(length(ind)>0){ne <- ne[-ind, ]}
  med_Ne_X<-c(med_Ne_X,round(median(ne$ne)))
}
rm(af);rm(cov);rm(ne);rm(sync);rm(i);rm(ind);rm(j);rm(nb_obs);rm(NB_trials);rm(Ne);rm(off)
#med_Ne_R1 <- med_Ne[1]
#med_Ne_R2 <- med_Ne[2]
#med_Ne_R3 <- med_Ne[3]
write.table(data.frame(med_Ne_auto,med_Ne_X),file="Ne.log", quote=F, sep='\t')


#4. adapted CMH-----------
cmh_p<-NULL
p_list<-NULL
locus<-NULL
###4.1 autosome-----------
for (chr in c("2L", "2R", "3L", "3R","4")){
  sync<-read.sync(file=paste("02noNA.no_0var.F0_F60.",chr,".filtered.cons.q20.no_0cov.sync",sep=""), 
                  gen = rep(input_gen, each = input_rep), repl = rep(1:input_rep, 2), polarization = "rising")
  af=as.matrix(poolSeq::af(sync, gen = rep(input_gen,each=input_rep), repl = rep(1:input_rep,2)))
  cov=as.matrix(poolSeq::coverage(sync, gen = rep(input_gen,each=input_rep), repl = rep(1:input_rep,2)))
  p.values<-adapted.cmh.test(freq = af, coverage = cov, Ne = med_Ne_auto, gen = input_gen, repl = 1:input_rep, poolSize = rep(600,ncol(af)))
  p_list<-c(p_list, p.values)
  locus<-rbind(locus, splitLocusID(rownames(af)))}
  adj_p<-p.adjust(p_list,method = "fdr", n=length(p_list))
  cmh_p<-cbind(p_list,adj_p,locus)
  #cmh_p<-rbind(cmh_p,p.values)

###4.2 X chromosome------
sync<-read.sync(file="02noNA.no_0var.F0_F60.X.filtered.cons.q20.no_0cov.sync", 
                gen = rep(input_gen, each = input_rep), repl = rep(1:input_rep, 2), polarization = "rising")
af=as.matrix(poolSeq::af(sync, gen = rep(input_gen,each=input_rep), repl = rep(1:input_rep,2)))
cov=as.matrix(poolSeq::coverage(sync, gen = rep(input_gen,each=input_rep), repl = rep(1:input_rep,2)))
p_list<-adapted.cmh.test(freq = af, coverage = cov, Ne = med_Ne_X, gen = input_gen, repl = 1:input_rep, poolSize = rep(600,ncol(af)))
adj_p<-p.adjust(p_list,method = "fdr", n=length(p_list))
p_list<-cbind(p_list,adj_p,splitLocusID(rownames(af)))
cmh_p<-rbind(cmh_p,p_list)
colnames(cmh_p)<-c("p","adj_p","chr","pos")
#adjusted_cmh_p<-p.adjust(cmh_p$p.values, method = "fdr", n = length(cmh_p$p.values))
#result<-cbind(adjusted_cmh_p,cmh_p);colnames(result)<-c("adj_p", "p", "chr", "pos")
write.table(cmh_p, "CMH.result", quote = F, sep = "\t", row.names = F)
#hist(adjusted_cmh_p)
#d1<-subset(cmh_p,chr == "2L")
#plot(d1$pos,-log(d1$adj_p))

#5. FET test----------
fet_p<-NULL
locus<-NULL
###5.1 autosome----------
p.values_result<-NULL
for (t_rep in 1:input_rep){
  af<-NULL
  cov<-NULL
  #p.values_result<-NULL
  for (chr in c("2L", "2R", "3L", "3R","4")){
    sync<-read.sync(file=paste("02noNA.no_0var.F0_F60.",chr,".filtered.cons.q20.no_0cov.sync", sep=""),
                  gen = rep(input_gen, each = input_rep), repl = rep(1:input_rep, 2), polarization = "rising")
    #p.values_result<-NULL
#  for (t_rep in 1:input_rep){
    af=rbind(af, as.matrix(poolSeq::af(sync, gen = rep(input_gen,each=input_rep), repl = rep(t_rep,2))))
    cov=rbind(cov, as.matrix(poolSeq::coverage(sync, gen = rep(input_gen,each=input_rep), repl = rep(t_rep,2))))}
    p.values<-adapted.chisq.test(freq = af, coverage = cov, Ne = med_Ne_auto[t_rep], gen = input_gen, poolSize = rep(600,ncol(af)))
  adj_p<-p.adjust(p.values,method = "fdr", n=length(p.values))
  p.values_result<-cbind(p.values_result, p.values, adj_p)}
fet_p<-cbind(p.values_result, splitLocusID(rownames(af)))
#  p.values_result<-cbind(p.values_result,splitLocusID(rownames(af)))
#  fet_p<-p.values_result

###5.2 X chromsome---------
sync<-read.sync(file="02noNA.no_0var.F0_F60.X.filtered.cons.q20.no_0cov.sync", 
                gen = rep(input_gen, each = input_rep), repl = rep(1:input_rep, 2), polarization = "rising")
p.values_result<-NULL
for (t_rep in 1:input_rep){
  af=as.matrix(poolSeq::af(sync, gen = rep(input_gen,each=input_rep), repl = rep(t_rep,2)))
  cov=as.matrix(poolSeq::coverage(sync, gen = rep(input_gen,each=input_rep), repl = rep(t_rep,2)))
  p.values<-adapted.chisq.test(freq = af, coverage = cov, Ne = med_Ne_X[t_rep], gen = input_gen, poolSize = rep(600,ncol(af)))
  adj_p<-p.adjust(p.values,method = "fdr", n=length(p.values))
  p.values_result<-cbind(p.values_result, p.values, adj_p)}
p.values_result<-cbind(p.values_result,splitLocusID(rownames(af)))
#p.values<-cbind(p.values,adj_p,splitLocusID(rownames(af)))
fet_p<-rbind(fet_p,p.values_result)
colnames(fet_p)<-c("r1_p","r1_adj_p","r2_p","r2_adj_p","r3_p","r3_adj_p","r4_p","r4_adj_p","r5_p","r5_adj_p","chr","pos")
fet_p_neworder<-fet_p[,c("chr","pos","r1_p", "r2_p","r3_p","r4_p","r5_p","r1_adj_p","r2_adj_p","r3_adj_p","r4_adj_p","r5_adj_p")]
write.table(fet_p_neworder, "FET.result", quote = F, sep = "\t", row.names = F)


