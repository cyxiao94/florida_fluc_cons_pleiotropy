
library(data.table)
library(ACER)
library(poolSeq)
library(scales)
library(ggplot2)
library(ggsignif)

#function of cmh_test
#organize data
load_data=function(num_snp,pleio,set){
  ##load af
      rep1=fread(paste0("../03results/03results_",num_snp,"SNP_",pleio,"_set",set,"_p1.freq.gz"), select=c("Pos","F0.freq","F60.freq"));colnames(rep1)=c("pos","F0.r1.freq","F60.r1.freq")
      rep2=fread(paste0("../03results/03results_",num_snp,"SNP_",pleio,"_set",set,"_p2.freq.gz"), select=c("Pos","F0.freq","F60.freq"));colnames(rep2)=c("pos","F0.r2.freq","F60.r2.freq")
      rep3=fread(paste0("../03results/03results_",num_snp,"SNP_",pleio,"_set",set,"_p3.freq.gz"), select=c("Pos","F0.freq","F60.freq"));colnames(rep3)=c("pos","F0.r3.freq","F60.r3.freq")
      rep4=fread(paste0("../03results/03results_",num_snp,"SNP_",pleio,"_set",set,"_p4.freq.gz"), select=c("Pos","F0.freq","F60.freq"));colnames(rep4)=c("pos","F0.r4.freq","F60.r4.freq")
      rep5=fread(paste0("../03results/03results_",num_snp,"SNP_",pleio,"_set",set,"_p5.freq.gz"), select=c("Pos","F0.freq","F60.freq"));colnames(rep5)=c("pos","F0.r5.freq","F60.r5.freq")
      rep6=fread(paste0("../03results/03results_",num_snp,"SNP_",pleio,"_set",set,"_p6.freq.gz"), select=c("Pos","F0.freq","F60.freq"));colnames(rep6)=c("pos","F0.r6.freq","F60.r6.freq")
      rep7=fread(paste0("../03results/03results_",num_snp,"SNP_",pleio,"_set",set,"_p7.freq.gz"), select=c("Pos","F0.freq","F60.freq"));colnames(rep7)=c("pos","F0.r7.freq","F60.r7.freq")
      rep8=fread(paste0("../03results/03results_",num_snp,"SNP_",pleio,"_set",set,"_p8.freq.gz"), select=c("Pos","F0.freq","F60.freq"));colnames(rep8)=c("pos","F0.r8.freq","F60.r8.freq")
      rep9=fread(paste0("../03results/03results_",num_snp,"SNP_",pleio,"_set",set,"_p9.freq.gz"), select=c("Pos","F0.freq","F60.freq"));colnames(rep9)=c("pos","F0.r9.freq","F60.r9.freq")
      rep10=fread(paste0("../03results/03results_",num_snp,"SNP_",pleio,"_set",set,"_p10.freq.gz"), select=c("Pos","F0.freq","F60.freq"));colnames(rep10)=c("pos","F0.r10.freq","F60.r10.freq")
      #combine data into a single matrix
      list_af=list(rep1,rep2,rep3,rep4,rep5,rep6,rep7,rep8,rep9,rep10)
      af<-Reduce(function(...) merge(...,by="pos",all=TRUE), list_af)
      pos=af$pos
      af=as.matrix(af[,2:ncol(af)]);rownames(af)=pos
      #af=af/2
      return(af)
}

#estimate Ne
cal_ne=function(af,cov){

  ##estimate Ne
  NE=NULL
  for (i in 1:10){
    NB_trials=100
    off=1
    nb_obs=dim(af)[1]
    ne=NULL
    for (j in 1:NB_trials){
      set.seed(off)
      ind=sample(x=1:nb_obs, size=1000)
      ne=rbind(ne, data.frame(trail=j,
                              ne=estimateNe(p0 =af[ind,(2*i-1)] , pt =af[ind,2*i] , cov0 =cov[ind,1] , covt =cov[ind,2] ,
                              t = 60, ploidy = 2, truncAF = 0.001, method = c("P.planI"), poolSize = c(1000,1000), Ncensus = 1000)))
      off=off+1
    }
    #Ne=median(ne$ne)
    ind=which(is.na(ne$ne) | ne$ne<0)
    if (length(ind)>0){ne=ne[-ind,]}
    NE=c(NE,round(median(ne$ne)))
  }
  
  
  ##filled NA estimate to median
  median_ne_fluc=median(na.omit(NE[1:5]))
  median_ne_cons=median(na.omit(NE[6:10]))
  NE[which(is.na(NE[1:5]))]=median_ne_fluc
  NE[5+which(is.na(NE[6:10]))]=median_ne_cons
  
  return(NE)
}

#perform CMH-test
cal_cmh = function(af,cov,NE){
  ##cmh test
  cmh_fluc=adapted.cmh.test(freq = af[,1:10], coverage = cov, order = 0, Ne = NE[1:5], gen = c(0,60), repl = 1:5)
  adj_cmh_fluc=p.adjust(cmh_fluc, method = "BH")

  cmh_cons=adapted.cmh.test(freq = af[,11:20], coverage = cov, order = 0, Ne = NE[6:10], gen = c(0,60), repl = 1:5)
  adj_cmh_cons=p.adjust(cmh_cons, method = "BH")


  cmh_result=data.table(pos=pos,
                        p_fluc=adj_cmh_fluc,
                        p_cons=adj_cmh_cons)


  return(cmh_result)
}



#function to calculate number of shared SNPs

num_sig = function(input_matrix){
  fluc=length(which(input_matrix$p_fluc<0.05))
  cons=length(which(input_matrix$p_cons<0.05))
  both=length(which(input_matrix$p_fluc<0.05 & input_matrix$p_cons<0.05))
  or=length(which(input_matrix$p_fluc<0.05 | input_matrix$p_cons<0.05))
  return(c(fluc,cons,both,or))
}


#perform CMH-test
suppressWarnings({
sum_info=data.frame()
NE_list=data.frame()
for (pleio in c("pleio")){
  for (num_snp in c("10000")){
    for (set in 1:100){
      af=load_data(num_snp,pleio,set)
      pos=rownames(af)
      #coverage matrix, were set to 1000 for every SNP
      cov=matrix(1000,nrow=nrow(af), ncol=ncol(af)/2);rownames(cov)=pos;#colnames(cov)=colnames(af)
      NE=cal_ne(af,cov);
      NE_list=rbind(NE_list,c(set,num_snp,pleio,NE))
      cmh_result=cal_cmh(af,cov,NE)
      fwrite(cmh_result,file = paste0("cmh_result/",num_snp,"_",pleio,"_set",set,".cmh"), quote = FALSE, sep="\t", na = "NA")
      sum_info=rbind(sum_info,c(set,num_snp,pleio,num_sig(cmh_result)))
    }
  }
}
colnames(sum_info)=c("set","num_SNP","pleiotropy","fluc","cons","both","total_significant")
fwrite(sum_info, file="Number_significant.txt",quote = FALSE, sep="\t")

colnames(NE_list)=c("set","num_SNP","pleiotropy","fluc_r1","fluc_r2","fluc_r3","fluc_r4","fluc_r5","cons_r1","cons_r2","cons_r3","cons_r4","cons_r5")
fwrite(NE_list, file="Estimated_Ne.txt", quote = FALSE, sep="\t", na="NA")

})


