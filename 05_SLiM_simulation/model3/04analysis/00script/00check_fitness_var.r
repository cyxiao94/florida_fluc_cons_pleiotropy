library(data.table)
library(Rmisc)

cal_sd=function(num_snp,ispleiotropy,per_lab,set,pop){
    #filename1="../03results/03results_10000SNP_nopleio_0.001lab_set1_p1.pheno"
    filename1=paste0("../03results/03results_",num_snp,"SNP_",ispleiotropy,"_",per_lab,"lab_set",set,"_",pop,".pheno")
    f1=fread(filename1)

    #pheno_lab
    pheno_f1_lab=summarySE(f1, measurevar= "phenotype_lab",, groupvars = "generation")[,1:4]
    pheno_f1_lab$F_pheno_lab=(pheno_f1_lab$sd/(pheno_f1_lab$sd[1]))^2
    colnames(pheno_f1_lab)=c("generation","N","mean_pheno_lab","sd_pheno_lab","F_pheno_lab")

    #pheno_23
    pheno_f1_23=summarySE(f1, measurevar= "phenotype_23",, groupvars = "generation")[,1:4]
    pheno_f1_23$F_pheno_23=(pheno_f1_23$sd/(pheno_f1_23$sd[1]))^2
    colnames(pheno_f1_23)=c("generation","N","mean_pheno_23","sd_pheno_23","F_pheno_23")

    #pheno_2818
    pheno_f1_2818=summarySE(f1, measurevar= "phenotype_2818",, groupvars = "generation")[,1:4]
    pheno_f1_2818$F_pheno_2818=(pheno_f1_2818$sd/(pheno_f1_2818$sd[1]))^2
    colnames(pheno_f1_2818)=c("generation","N","mean_pheno_2818","sd_pheno_2818","F_pheno_2818")

    #fit_lab
    fit_f1_lab=summarySE(f1, measurevar= "effects_lab",, groupvars = "generation")[,1:4]
    fit_f1_lab$F_fit_lab=(fit_f1_lab$sd/(fit_f1_lab$sd[1]))^2
    colnames(fit_f1_lab)=c("generation","N","mean_fit_lab","sd_fit_lab","F_fit_lab")

    #fit_23
    fit_f1_23=summarySE(f1, measurevar= "effects_23",, groupvars = "generation")[,1:4]
    fit_f1_23$F_fit_23=(fit_f1_23$sd/(fit_f1_23$sd[1]))^2
    colnames(fit_f1_23)=c("generation","N","mean_fit_23","sd_fit_23","F_fit_23")

    #fit_2818
    fit_f1_2818=summarySE(f1, measurevar= "effects_2818",, groupvars = "generation")[,1:4]
    fit_f1_2818$F_fit_2818=(fit_f1_2818$sd/(fit_f1_2818$sd[1]))^2
    colnames(fit_f1_2818)=c("generation","N","mean_fit_2818","sd_fit_2818","F_fit_2818")

    #merge
    merged_dt <- Reduce(function(x, y) merge(x, y, by = c("generation","N"), all = TRUE), list(pheno_f1_lab, pheno_f1_23, pheno_f1_2818, fit_f1_lab, fit_f1_23, fit_f1_2818))
    merged_dt$num_snp=num_snp
    merged_dt$ispleiotropy=ispleiotropy
    merged_dt$percent_lab=per_lab
    merged_dt$set=set
    merged_dt$pop=pop
    
    return(merged_dt)
}


#num_snp,ispleiotropy,linkage,set,pop
#num_snp=10000
#ispleiotropy='pleio'
#linkage="LD"
#set="1"
#pop="p1"
#a=cal_sd(num_snp,ispleiotropy,linkage,set,pop)

list_pheno_sum=data.table()
num_snp=10000
for (ispleiotropy in c('nopleio')){
    for (per_lab in c('0.001')){
        for (set in seq(1,100)){
            for (pop in c(paste0("p",seq(1,10)))){
                temp_f=cal_sd(num_snp,ispleiotropy,per_lab,set,pop)
                list_pheno_sum=rbind(list_pheno_sum,temp_f)
            }
        }
    }
}
fwrite(list_pheno_sum, file="summary_phenotypic_variation.txt",sep="\t",quote=FALSE)
#list_pheno_sum=rbind(list_pheno_sum,a)
