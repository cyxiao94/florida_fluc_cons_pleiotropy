library(data.table)
library(ggplot2)
library(dplyr)
library(ggpubr)

for (model in c("model1","model2","model3","model4")){
f=fread(paste0('../',model,'/04analysis/summary_phenotypic_variation.txt'))

f$generation=f$generation-1
f$generation=factor(f$generation, levels=c(0,1,2,3,4,5,10,20,30,40,50,60))
f$ispleiotropy=factor(f$ispleiotropy)
f$model="model1"
f[, population := ifelse(pop %in% c("p1", "p2", "p3", "p4", "p5"), "fluctuating", "constant")]

f$population=factor(f$population)

#measurement error
#f[, sampled_mean_fit_lab := sapply(1:.N, function(i) mean(rnorm(100, mean_fit_lab[i], sd_fit_lab[i])))]
#f[, sampled_mean_fit_23 := sapply(1:.N, function(i) mean(rnorm(100, mean_fit_23[i], sd_fit_23[i])))]
#f[, sampled_mean_fit_2818 := sapply(1:.N, function(i) mean(rnorm(100, mean_fit_2818[i], sd_fit_2818[i])))]

#f$percent_lab=factor(f$percent_lab)

#plot mean fitness in 23 and 2818 environment
#
p1=ggplot(f,aes(x=generation,y=mean_fit_lab*mean_fit_23, fill=population))+
geom_boxplot()+
scale_fill_manual(values=c("#4fc0c6","#e8a11c"))+
#stat_compare_means(aes(group=population), label="p.signif", method="wilcox.test")+
ylim(0,3)+
labs(title="23 CGE",y="Fitness")+
theme_bw()
p2=ggplot(f,aes(x=generation,y=mean_fit_lab*mean_fit_2818, fill=population))+
geom_boxplot()+
scale_fill_manual(values=c("#4fc0c6","#e8a11c"))+
ylim(0,3)+
labs(title="2818 CGE",y="Fitness")+
theme_bw()
p=ggarrange(p1,p2, nrow=2,common.legend=TRUE)
ggsave(paste0("00_",model,"_fitness.pdf"),plot=p, width=4, heigh=7)

}


