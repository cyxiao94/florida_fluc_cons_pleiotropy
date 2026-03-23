
library(data.table)
#library(ACER)
#library(poolSeq)
library(scales)
library(ggplot2)
library(ggsignif)

#plot
#significant SNPs
#setwd("04analysis")
info=NULL
for (model in c(1,2,3,4)){
sum_info=fread(paste0("../model",model,"/04analysis/Number_significant.txt"))
sum_info[,"model"]=paste0("model",model)
sum_info$fluc=as.numeric(sum_info$fluc)
sum_info$cons=as.numeric(sum_info$cons)
sum_info$both=as.numeric(sum_info$both)
sum_info$total_significant=as.numeric(sum_info$total_significant)
#sum_info$cate=paste(sum_info$pleiotropy, sum_info$num_SNP,sep="_")
#sum_info$cate=factor(sum_info$cate)
#sum_info$cate=factor(sum_info$cate, levels = c("pleio_LD","nopleio_LD","pleio_LE","nopleio_LE"))
sum_info$proportion=sum_info$both/sum_info$total_significant
temp_info=subset(sum_info, select=c("set","model","proportion"))
info=rbind(info,temp_info)
}
info$model=factor(info$model)

p1=ggplot(info, aes(x=model,y=proportion))+
    geom_boxplot(outlier.size = 0.5)+
    geom_signif(comparisons=list(c("model2","model3"),
                                 c("model2","model4"),
                                 c("model3","model4")),
                map_signif_level=TRUE,
                color="darkred",
                #y_position=c(0.15,0.19,0.125),
                y_position=c(0.1,0.11,0.09),
                tip_length=0.005)+
    #ylim(0,0.25)+
    ylim(0,0.15)+
    labs(y="Proportion of SNPs shared",x="Model")+
    #annotate("text",x=1,y=-0.5,size=2,vjust=1,
    #        label="Model1\nEqual architecture\nEqual selection\nNon-pleiotropy")+
    scale_x_discrete(labels=c("Model1\nEqual architecture\nEqual selection\nNon-pleiotropy",
                              "Model2\nEqual architecture\nLab selection\nNon-pleiotropy",
                              "Model3\nUnequal architecture\nLab selection\nNon-pleiotropy",
                              "Model4\nEqual architecture\nLab selection\nPleiotropy"))+
    coord_cartesian(clip = "off") +
    theme_bw()+
    theme(legend.position="none",
          axis.ticks.x=element_blank())
ggsave("01Proportion_shared_SNP.pdf", width=7, height=4)






