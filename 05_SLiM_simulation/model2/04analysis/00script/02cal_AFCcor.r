library(ggplot2)
library(data.table)
library(corrplot)
#function of cmh_test
#organize data
cal_afc=function(num_snp,pleio,set){
    #filename "03results_10000SNP_nopleio_0.3333lab_set9_p9.freq"
  ##load af
      rep1=fread(paste0("../03results/03results_",num_snp,"SNP_",pleio,"_0.3333lab_set",set,"_p1.freq.gz"), select=c("Pos","F0.freq","F60.freq"));afc_r1=rep1[,3]-rep1[,2];
      rep2=fread(paste0("../03results/03results_",num_snp,"SNP_",pleio,"_0.3333lab_set",set,"_p2.freq.gz"), select=c("Pos","F0.freq","F60.freq"));afc_r2=rep2[,3]-rep2[,2]
      rep3=fread(paste0("../03results/03results_",num_snp,"SNP_",pleio,"_0.3333lab_set",set,"_p3.freq.gz"), select=c("Pos","F0.freq","F60.freq"));afc_r3=rep3[,3]-rep3[,2]
      rep4=fread(paste0("../03results/03results_",num_snp,"SNP_",pleio,"_0.3333lab_set",set,"_p4.freq.gz"), select=c("Pos","F0.freq","F60.freq"));afc_r4=rep4[,3]-rep4[,2]
      rep5=fread(paste0("../03results/03results_",num_snp,"SNP_",pleio,"_0.3333lab_set",set,"_p5.freq.gz"), select=c("Pos","F0.freq","F60.freq"));afc_r5=rep5[,3]-rep5[,2]
      rep6=fread(paste0("../03results/03results_",num_snp,"SNP_",pleio,"_0.3333lab_set",set,"_p6.freq.gz"), select=c("Pos","F0.freq","F60.freq"));afc_r6=rep6[,3]-rep6[,2]
      rep7=fread(paste0("../03results/03results_",num_snp,"SNP_",pleio,"_0.3333lab_set",set,"_p7.freq.gz"), select=c("Pos","F0.freq","F60.freq"));afc_r7=rep7[,3]-rep7[,2]
      rep8=fread(paste0("../03results/03results_",num_snp,"SNP_",pleio,"_0.3333lab_set",set,"_p8.freq.gz"), select=c("Pos","F0.freq","F60.freq"));afc_r8=rep8[,3]-rep8[,2]
      rep9=fread(paste0("../03results/03results_",num_snp,"SNP_",pleio,"_0.3333lab_set",set,"_p9.freq.gz"), select=c("Pos","F0.freq","F60.freq"));afc_r9=rep9[,3]-rep9[,2]
      rep10=fread(paste0("../03results/03results_",num_snp,"SNP_",pleio,"_0.3333lab_set",set,"_p10.freq.gz"), select=c("Pos","F0.freq","F60.freq"));afc_r10=rep10[,3]-rep10[,2]
      #combine data into a single matrix
      dt_afc=data.table(cbind(rep1[,1],afc_r1,afc_r2,afc_r3,afc_r4,afc_r5,afc_r6,afc_r7,afc_r8,afc_r9,afc_r10))
      colnames(dt_afc)=c("pos",paste0("r",1:10))
      #list_af=list(rep1,rep2,rep3,rep4,rep5,rep6,rep7,rep8,rep9,rep10)
      #af<-Reduce(function(...) merge(...,by="pos",all=TRUE), list_af)
      #pos=af$pos
      #af=as.matrix(af[,2:ncol(af)]);rownames(af)=pos
      #af=af/2
      return(dt_afc)
}


summary_cor=function(cor_matrix){
    sum_cor=NULL
    for (i in 1:9){
        for (j in (i+1):10){
            if (i %in% 1:5){env1="fluc"}else{env1="cons"}
            if (j %in% 1:5){env2="fluc"}else{env2="cons"}
            if (env1==env2){comparison="within"}else{comparison="between"}
            pair1=paste0("r",i);
            pair2=paste0("r",j)
            cor1=cor_matrix[i,j]
            sum_cor=rbind(sum_cor,c(cor1,pair1,pair2,env1,env2,comparison))
        }
    }
    colnames(sum_cor)=c("correlation","pair1","pair2","env1","env2","comparison")
    #sum_cor$correlation=as.numeric(sum_cor$correlation)
    return(data.frame(sum_cor))
}


summary_cor_all_SNP=data.frame()
summary_cor_sig_SNP=data.frame()

for (num_snp1 in c(10000)){
    for (pleio1 in c("nopleio")){
        for (set1 in 1:100){
            afc=cal_afc(num_snp1,pleio1,set1)
            #all SNPs
            cor_matrix=cor(afc[,2:11], use="pairwise.complete.obs",method="spearman")
            sum_cor=summary_cor(cor_matrix)
            sum_cor$num_snp=num_snp1
            sum_cor$pleio=pleio1
            sum_cor$set=set1
            summary_cor_all_SNP=rbind(summary_cor_all_SNP, sum_cor)

            #significant SNPs only
            #cmhfile "10000_nopleio_0.3333lab_set85.cmh"
            cmh=fread(paste0('cmh_result/',num_snp1,'_',pleio1,'_0.3333lab_set',set1,'.cmh'))
            sig_snp=subset(cmh,p_fluc<=0.05 | p_cons<=0.05)
            cor_matrix_sig=cor(afc[pos %in% sig_snp$pos, .SD, .SDcols=2:11], use="pairwise.complete.obs", method="spearman")
            sum_cor_sig=summary_cor(cor_matrix_sig)
            sum_cor_sig$num_snp=num_snp1
            sum_cor_sig$pleio=pleio1
            sum_cor_sig$set=set1
            summary_cor_sig_SNP=rbind(summary_cor_sig_SNP, sum_cor_sig)
        }
    }
}
fwrite(summary_cor_all_SNP,"summary_correlation_all_SNP.txt",quote=FALSE,sep="\t")
fwrite(summary_cor_sig_SNP,"summary_correlation_sig_SNP.txt",quote=FALSE,sep="\t")
            #corrplot(cor_matrix_sig,diag=FALSE,is.corr=FALSE)
