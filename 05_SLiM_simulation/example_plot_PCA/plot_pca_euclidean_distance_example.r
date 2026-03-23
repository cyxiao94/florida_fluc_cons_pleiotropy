
library(data.table)
#library(ACER)
library(factoextra)
library(reshape2)
#library(poolSeq)
library(scales)
library(ggplot2)
library(ggsignif)
library(dplyr)

#function of cmh_test
#organize data
load_data=function(model, set){
  ##load af
      rep1 <- fread(Sys.glob(paste0("../",model,"/03results/03results_*_set",set,"_p1.freq")))
      colnames(rep1)=c("pos","F0.r1.freq","F10.r1.freq","F20.r1.freq","F30.r1.freq","F40.r1.freq","F50.r1.freq","F60.r1.freq") # nolint: line_length_linter, line_length_linter.
      rep2 <- fread(Sys.glob(paste0("../",model,"/03results/03results_*_set",set,"_p2.freq"))) # nolint: line_length_linter.
      colnames(rep2)=c("pos","F0.r2.freq","F10.r2.freq","F20.r2.freq","F30.r2.freq","F40.r2.freq","F50.r2.freq","F60.r2.freq")
      rep3 <- fread(Sys.glob(paste0("../",model,"/03results/03results_*_set",set,"_p3.freq")))
      colnames(rep3)=c("pos","F0.r3.freq","F10.r3.freq","F20.r3.freq","F30.r3.freq","F40.r3.freq","F50.r3.freq","F60.r3.freq")
      rep4 <- fread(Sys.glob(paste0("../",model,"/03results/03results_*_set",set,"_p4.freq")))
      colnames(rep4)=c("pos","F0.r4.freq","F10.r4.freq","F20.r4.freq","F30.r4.freq","F40.r4.freq","F50.r4.freq","F60.r4.freq")
      rep5 <- fread(Sys.glob(paste0("../",model,"/03results/03results_*_set",set,"_p5.freq")))
      colnames(rep5)=c("pos","F0.r5.freq","F10.r5.freq","F20.r5.freq","F30.r5.freq","F40.r5.freq","F50.r5.freq","F60.r5.freq")
      rep6 <- fread(Sys.glob(paste0("../",model,"/03results/03results_*_set",set,"_p6.freq")))
      colnames(rep6)=c("pos","F0.r6.freq","F10.r6.freq","F20.r6.freq","F30.r6.freq","F40.r6.freq","F50.r6.freq","F60.r6.freq")
      rep7 <- fread(Sys.glob(paste0("../",model,"/03results/03results_*_set",set,"_p7.freq")))
      colnames(rep7)=c("pos","F0.r7.freq","F10.r7.freq","F20.r7.freq","F30.r7.freq","F40.r7.freq","F50.r7.freq","F60.r7.freq")
      rep8 <- fread(Sys.glob(paste0("../",model,"/03results/03results_*_set",set,"_p8.freq")))
      colnames(rep8)=c("pos","F0.r8.freq","F10.r8.freq","F20.r8.freq","F30.r8.freq","F40.r8.freq","F50.r8.freq","F60.r8.freq")
      rep9 <- fread(Sys.glob(paste0("../",model,"/03results/03results_*_set",set,"_p9.freq")))
      colnames(rep9)=c("pos","F0.r9.freq","F10.r9.freq","F20.r9.freq","F30.r9.freq","F40.r9.freq","F50.r9.freq","F60.r9.freq")
      rep10 <- fread(Sys.glob(paste0("../",model,"/03results/03results_*_set",set,"_p10.freq")))
      colnames(rep10)=c("pos","F0.r10.freq","F10.r10.freq","F20.r10.freq","F30.r10.freq","F40.r10.freq","F50.r10.freq","F60.r10.freq")
      #combine data into a single matrix
      list_af=list(rep1,rep2,rep3,rep4,rep5,rep6,rep7,rep8,rep9,rep10)
      af<-Reduce(function(...) merge(...,by="pos",all=TRUE), list_af)
      pos=af$pos
      af=as.matrix(af[,2:ncol(af)]);rownames(af)=pos
      #af=af/2
      return(af)
}



num_snp=10000
pleio="pleio"
set=1
sum_pca=data.frame()

#model3_set7
model="model3"
num_snp=10000
pleio="nopleio"
set=7

my_freqf=as.data.frame(load_data(model,set))
dat=2*asin(sqrt(my_freqf))
clmnames<-colnames(dat)
sampleIDs <- sub(".freq","",clmnames)
splits <- unlist(strsplit(sampleIDs, "[.]"))
generation <- as.factor(splits[grep("F",splits)])
rep<- factor(splits[grep("r",splits)], levels=c("r1","r2","r3","r4","r5",
                                                "r6","r7","r8","r9","r10"))
environment=factor(c(rep(c("Ancestral","Fluctuating","Fluctuating","Fluctuating","Fluctuating","Fluctuating","Fluctuating"),5),
                  rep(c("Ancestral","Constant","Constant","Constant","Constant","Constant","Constant"),5)))
t_dat<-as.data.frame(t(dat))
row.names(t_dat)<-sampleIDs

pcadata=na.omit(t_dat)
pca.res=prcomp(pcadata, retx=TRUE, center=TRUE, scale. = TRUE)
pca.var=get_eigenvalue(pca.res)$variance.percent
saveRDS(pca.res$x,file=paste0(model,"_",num_snp,"targets_",pleio,"_0.001lab_set",set,"_loadings.rds"))
saveRDS(pca.var,file=paste0(model,"_",num_snp,"targets_",pleio,"_0.001lab_set",set,"_var.rds"))


pca_loading1=as.data.frame(pca.res$x[paste0(rep(c("F0.r","F60.r"),10),rep(1:10,each=2)),1:2])
pca_loading1$group=c(rep(c("Ancestral","Fluctuating"),5),
                      rep(c("Ancestral","Constant"),5))
# Calculate centroids for each group， only use F60 data
centroids <- pca_loading1 %>%
  group_by(group) %>%
  summarize(PC1 = mean(PC1), PC2 = mean(PC2))

# Calculate Euclidean distance between groups
distances = as.matrix(dist(centroids[, c("PC1", "PC2")]))
colnames(distances) = centroids$group
rownames(distances) = centroids$group
sum_pca=rbind(sum_pca,c(paste0("set",set),distances["Ancestral","Fluctuating"],distances["Ancestral","Constant"],distances["Fluctuating","Constant"]))
      
p1=ggplot()+
  #geom_point(data=as.data.frame(pca.res$x),aes(x=PC1, y=PC2,col=environment,shape=rep),size=1,alpha=0.5)+
  geom_point(data=as.data.frame(pca.res$x),aes(x=PC1, y=PC2,col=environment),size=1,alpha=0.5, shape=2)+
  #stat_ellipse()
  #scale_shape_manual(values=c(15,16,17,18,1,15,16,17,18,1))+
  xlab(paste0("PC1(",round(pca.var[1],2),"%)"))+
  ylab(paste0("PC2(",round(pca.var[2],2),"%)"))+
  stat_ellipse(data=subset(pca_loading1,group!="Ancestral"),aes(x=PC1,y=PC2,color=group),type="norm", linetype="dashed")+
  geom_point(data=centroids, aes(x=PC1, y=PC2, color=group), shape=17,size=2)+
  geom_segment(data=centroids, aes(x=PC1,y=PC2,xend=PC1[c(2,3,1)],yend=PC2[c(2,3,1)]),color="black",size=0.5)+
  theme(legend.position = "none")

ggsave(paste0(model,"plot_",num_snp,"targets_",pleio,"_0.001lab_set",set,".pdf"),plot=p1, width = 5, height = 5)

 
colnames(sum_pca)=c("set","anc_fluc","anc_cons","fluc_cons")
fwrite(sum_pca,paste0(model,"_Euclidean_distant_PC1_PC2.txt"),quote=FALSE,sep="\t")


#model4_set16
model="model4"
num_snp=10000
pleio="pleio"
set=16

my_freqf=as.data.frame(load_data(model,set))
dat=2*asin(sqrt(my_freqf))
clmnames<-colnames(dat)
sampleIDs <- sub(".freq","",clmnames)
splits <- unlist(strsplit(sampleIDs, "[.]"))
generation <- as.factor(splits[grep("F",splits)])
rep<- factor(splits[grep("r",splits)], levels=c("r1","r2","r3","r4","r5",
                                                "r6","r7","r8","r9","r10"))
environment=factor(c(rep(c("Ancestral","Fluctuating","Fluctuating","Fluctuating","Fluctuating","Fluctuating","Fluctuating"),5),
                  rep(c("Ancestral","Constant","Constant","Constant","Constant","Constant","Constant"),5)))
t_dat<-as.data.frame(t(dat))
row.names(t_dat)<-sampleIDs

pcadata=na.omit(t_dat)
pca.res=prcomp(pcadata, retx=TRUE, center=TRUE, scale. = TRUE)
pca.var=get_eigenvalue(pca.res)$variance.percent
saveRDS(pca.res$x,file=paste0(model,"_",num_snp,"targets_",pleio,"_0.001lab_set",set,"_loadings.rds"))
saveRDS(pca.var,file=paste0(model,"_",num_snp,"targets_",pleio,"_0.001lab_set",set,"_var.rds"))


pca_loading1=as.data.frame(pca.res$x[paste0(rep(c("F0.r","F60.r"),10),rep(1:10,each=2)),1:2])
pca_loading1$group=c(rep(c("Ancestral","Fluctuating"),5),
                      rep(c("Ancestral","Constant"),5))
# Calculate centroids for each group， only use F60 data
centroids <- pca_loading1 %>%
  group_by(group) %>%
  summarize(PC1 = mean(PC1), PC2 = mean(PC2))

# Calculate Euclidean distance between groups
distances = as.matrix(dist(centroids[, c("PC1", "PC2")]))
colnames(distances) = centroids$group
rownames(distances) = centroids$group
sum_pca=rbind(sum_pca,c(paste0("set",set),distances["Ancestral","Fluctuating"],distances["Ancestral","Constant"],distances["Fluctuating","Constant"]))
      
p1=ggplot()+
  #geom_point(data=as.data.frame(pca.res$x),aes(x=PC1, y=PC2,col=environment,shape=rep),size=1,alpha=0.5)+
  geom_point(data=as.data.frame(pca.res$x),aes(x=PC1, y=PC2,col=environment),size=1,alpha=0.5, shape=2)+
  #stat_ellipse()
  #scale_shape_manual(values=c(15,16,17,18,1,15,16,17,18,1))+
  xlab(paste0("PC1(",round(pca.var[1],2),"%)"))+
  ylab(paste0("PC2(",round(pca.var[2],2),"%)"))+
  stat_ellipse(data=subset(pca_loading1,group!="Ancestral"),aes(x=PC1,y=PC2,color=group),type="norm", linetype="dashed")+
  geom_point(data=centroids, aes(x=PC1, y=PC2, color=group), shape=17,size=2)+
  geom_segment(data=centroids, aes(x=PC1,y=PC2,xend=PC1[c(2,3,1)],yend=PC2[c(2,3,1)]),color="black",size=0.5)+
  theme(legend.position = "none")

ggsave(paste0(model,"plot_",num_snp,"targets_",pleio,"_0.001lab_set",set,".pdf"),plot=p1, width = 5, height = 5)

 
colnames(sum_pca)=c("set","anc_fluc","anc_cons","fluc_cons")
fwrite(sum_pca,paste0(model,"_Euclidean_distant_PC1_PC2.txt"),quote=FALSE,sep="\t")
