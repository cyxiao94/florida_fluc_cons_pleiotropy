library(ggplot2)
library(data.table)
library(ggsignif)
library(ggpubr)




sum_f_all=NULL
sum_f_sig=NULL
for (model in c("model1","model2","model3","model4")){

    #significant SNPs only
    f_sig=fread(paste0('../',model,'/04analysis/summary_correlation_sig_SNP.txt'))
    f_sig$model=model
    sum_f_sig=rbind(sum_f_sig,f_sig)
}

sum_f_sig$model=factor(sum_f_sig$model)


my_comparisons1=list(c("model2_between","model3_between"),
                    c("model2_between","model4_between"),
                    c("model3_between","model4_between"))

my_comparisons2=list(c("model2_within","model3_within"),
                    c("model2_within","model4_within"),
                    c("model3_within","model4_within"))




spaces=strrep(" ",24)




##plot AFC correlation
sum_f_sig$cate=paste(sum_f_sig$model,sum_f_sig$comparison,sep="_")
p2=ggplot(sum_f_sig,aes(x=cate,y=correlation,color=comparison))+
            geom_boxplot(outlier.size = 0.2)+
            ylim(-1,1)+
           geom_signif(comparisons = my_comparisons1, map_signif_level = TRUE,  # Use significance stars or symbols
                y_position = c(0.66, 0.78, 0.9),  textsize = 3,  tip_length=0.03,test = "wilcox.test",
                vjust = 0.2) +  
            geom_signif(comparisons = my_comparisons2, map_signif_level = TRUE,  # Use significance stars or symbols
                y_position = c(-0.95, -0.85, -0.75),  textsize = 3,  tip_length=-0.03,test = "wilcox.test",
                vjust = 3.0) +  
            xlab("")+
                        scale_x_discrete(labels=c(paste0(spaces,"Model1\n",spaces,"Equal architecture\n",spaces,"Equal selection\n",spaces,"Non-pleiotropy"),"",
                                      paste0(spaces,"Model2\n",spaces,"Equal architecture\n",spaces,"Lab selection\n",spaces,"Non-pleiotropy"),"",
                                      paste0(spaces,"Model3\n",spaces,"Unequal architecture\n",spaces,"Lab selection\n",spaces,"Non-pleiotropy"),"",
                                      paste0(spaces,"Model4\n",spaces,"Equal architecture\n",spaces,"Lab selection\n",spaces,"Pleiotropy"),""))+
            #facet_wrap(~per_lab, nrow=1)+
            theme_bw()+
            theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0),
                  axis.ticks.x = element_blank(),
                  legend.position=c(0.07,0.88))
#plot(p1)

ggsave("02summary_correlation_sig_SNP.pdf",plot=p2,width=8,height=5)


#check AFC correlation
##model2 vs model3
#wilcox.test(sum_f_sig[model =="model3" & comparison =="within",correlation],
#            sum_f_sig[model =="model2" & comparison =="within",correlation], alternative = "less")
#p-value=0.23

#wilcox.test(sum_f_sig[model =="model3" & comparison =="between",correlation],
#              sum_f_sig[model =="model2" & comparison =="between",correlation], alternative = "greater") 
#p-value=0.0005495

#compared with model2, models 3 have similar within treatment correlation and lower between treatment correlation, suggesting model3 show less divergence, opposite from what we observed in the experimental evolution


#model2 vs model 4
#wilcox.test(sum_f_sig[model =="model4" & comparison =="within",correlation],
#            sum_f_sig[model =="model3" & comparison =="within",correlation], alternative = "greater")
#p-value<2.2e-16

#wilcox.test(sum_f_sig[model =="model4" & comparison =="between",correlation],
#              sum_f_sig[model =="model3" & comparison =="between",correlation], alternative = "less") 
#p-value<2.2e-16
#compare with model2, model 4 have higher within treatment correlation and lower between treatment correlation, more divergence, consistent with experimental evolution observation

