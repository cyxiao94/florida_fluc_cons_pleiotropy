#population adapting to the hot fluctuating and hot constant environment,
#F0-F60

library(poolSeq)
library(ggplot2)
library(factoextra)
library(reshape2)
library(data.table)
library(sparsepca)

#setwd("/Users/cxiao/my_project/fluctuation_constant/13PCA_contributingmostSNP/sparse_PCA/test/")
syncfile="fluc_cons.sync.gz"

##5replicate
print("PCA analysis for merged sync file")
my_sync<-read.sync(syncfile, polarization = "rising",
                   gen = rep(c(0,10,20,30,40,50,60,61,71,81,91,101,111), each=5), 
                   repl = rep(1:5, 15))

my_freqf<-as.data.frame(af(my_sync, gen = rep(c(0,10,20,30,40,50,60,61,71,81,91,101,111), each=5),repl = rep(1:5),13))
my_freqf<-my_freqf[which(apply(my_freqf,1,var)!=0),]
colnames(my_freqf)<-paste(rep(c("E_base",rep("E_cons",6),rep("E_fluc",6)),5),
                          paste("F",rep(c(0,10,20,30,40,50,60,10,20,30,40,50,60),5),sep="_"),
                          paste("R",rep(1:5,each=13),sep=""),"freq",sep=".")
#transform frequencies and organize data
dat = 2*asin(sqrt(my_freqf))
clmnames<-colnames(dat)
sampleIDs <- sub(".freq","",clmnames)
splits <- unlist(strsplit(sampleIDs, "[.]"))
#generation <- as.factor(splits[grep("F",splits)])
#rep<- as.factor(splits[grep("R",splits)])
environment<-as.factor(splits[grep("E",splits)])
t_dat<-as.data.frame(t(dat))
row.names(t_dat)<-sampleIDs


pcadata<-na.omit(t_dat)#remove missing data
rm(my_sync,my_freqf,dat,clmnames, sampleIDs,splits,t_dat)
##1. sparse PCA_centered_scaled----
spca.res<-spca(pcadata, center = TRUE, scale = TRUE)
saveRDS(spca.res$scores, file="spca_score.rds")
dat_spca<-spca.res$scores
suminfo<-summary(spca.res)
p1<-ggplot(as.data.frame(spca.res$scores), aes(x=V1,y=V2, col=environment, shape = environment))+
  geom_point()+
  xlab(paste("PC1(", suminfo[3,1]*100,"%)",sep=""))+
  ylab(paste("PC2(", suminfo[3,2]*100,"%)",sep=""))
pdf("PCA.anc_cons_fluc.centered.scaled.pdf")
plot(p1)
dev.off()

data_loading<-spca.res$loadings[,1:10]
rownames(data_loading)<-colnames(pcadata)
colnames(data_loading)<-paste("PC",1:10,sep="")
write.table(data_loading,"loadings.centered.scaled.txt", quote = F, row.names = T, col.names = T, sep="\t")
write.table(suminfo,"suminfo.centered.scaled.txt", quote = F, row.names = T, col.names = T, sep="\t")
rm(spca.res, dat_spca, suminfo, p1, data_loading)

