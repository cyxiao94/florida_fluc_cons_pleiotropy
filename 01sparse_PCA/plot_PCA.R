library(poolSeq)
library(ggplot2)
library(factoextra)
library(reshape2)
library(sparsepca)
library(data.table)
#setwd("~/my_project/fluctuation_constant/20final_dataset/05sparse_PCA/sparse_PCA_new/centered_scaled/")
#mergefile="merge.sync"
spca<-readRDS("spca_score.rds")
dat_spca<-as.data.frame(spca)
suminfo<-fread("suminfo.centered.scaled.txt")

samples<-rownames(dat_spca)
environment<-rep(c("Ancestral", rep("Constant",6), rep("Fluctuating",6)),5)
generation<-paste("F",rep(c(0,10,20,30,40,50,60,10,20,30,40,50,60),5),sep="")
replicate<-as.factor(rep(1:5, each=13))
environment_generation=as.factor(paste(environment, generation, sep="_"))
t1<-data.frame(samples, environment, generation, replicate)

p1<-ggplot(dat_spca, aes(x=V1,y=V2,col=environment_generation, shape=replicate))+
  geom_point(cex=3)+
  geom_vline(aes(xintercept=0), color='grey', linetype='dashed')+
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid = element_blank())+
  scale_color_manual(values = c("#6e326e",
                                "#48cae4","#00b4d8","#0096c7","#0077b6","#023e8a","#03045e",
                                "#f7b267","#f79d65","#f4845f","#f27059","#f25c54","#800f2f"))+
  #scale_color_manual(values = c("#ff8636",
  #                              "#a2ce69","#45d38a","#06d6a0","#05b184","#059b74","#056e52",
  #                              "#e37593","#ef476f","#962d46","#782438","#5a1b2a","#3c121c"))+
  scale_shape_manual(values = c(15,16,17,18,3))+
  xlab(paste("PC1(", suminfo[3,PC1]*100,"%)",sep=""))+
  ylab(paste("PC2(", suminfo[3,PC2]*100,"%)",sep=""))+
  guides(shape=guide_legend(order=1), col=guide_legend(order=2))

#tiff('sparsePCA_allsample_allSNP.tiff', width=6, height=5, units="cm", res=1000)
pdf('sparsePCA_allsample_allSNP.pdf', width=7.5, height=6)
plot(p1)
dev.off()
print("finished, congratulations!")

