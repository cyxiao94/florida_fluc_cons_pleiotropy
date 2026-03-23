library(data.table)
library(dplyr)



sum_pca_pc1=NULL
sum_pca_pc1_pc2=NULL
sum_pca_pc1_20=NULL
for (num_snp in c("10000")){
  for (pleio in c("pleio")){
    for (set in seq(1,100)){
      pca_loading=readRDS(paste0("plot_pca_new/",num_snp,"targets_",pleio,"_set",set,"_loadings.rds"))
      pca_var=readRDS(paste0("plot_pca_new/",num_snp,"targets_",pleio,"_set",set,"_var.rds"))
      pca_var=pca_var/100


      clmnames<-colnames(pca_loading)
      sampleIDs <- rownames(pca_loading)
      splits <- unlist(strsplit(sampleIDs, "[.]"))
      generation <- as.factor(splits[grep("F",splits)])
      rep<- factor(splits[grep("r",splits)], levels=c("r1","r2","r3","r4","r5",
                                                "r6","r7","r8","r9","r10"))
      environment=factor(c(rep(c("Ancestral","Fluctuating","Fluctuating","Fluctuating","Fluctuating","Fluctuating","Fluctuating"),5),
                        rep(c("Ancestral","Constant","Constant","Constant","Constant","Constant","Constant"),5)))
      


      pca_loading1=as.matrix(pca_loading[paste0(rep(c("F0.r","F60.r"),10),rep(1:10,each=2)),])
      #pca_loading1$group=c(rep(c("Ancestral","Fluctuating"),5),
      #                     rep(c("Ancestral","Constant"),5))
      weighted_loadings=as.data.frame(sweep(pca_loading1,2,pca_var,`*`))
      weighted_loadings$group=c(rep(c("Ancestral","Fluctuating"),5),
                           rep(c("Ancestral","Constant"),5))
      # Calculate centroids for each group
      centroids <- weighted_loadings %>%
        group_by(group) %>%
        summarize(across(everything(), \(x) mean(x, na.rm = TRUE)))

      ##distant for pc1 only
      # Calculate Euclidean distance between groups
      distances = as.matrix(dist(centroids[, c("PC1")]))
      colnames(distances) = centroids$group
      rownames(distances) = centroids$group
      sum_pca_pc1=rbind(sum_pca_pc1,c(paste0("set",set),distances["Ancestral","Fluctuating"],distances["Ancestral","Constant"],distances["Fluctuating","Constant"]))
      
      ##distant for pc1 and pc2
      # Calculate Euclidean distance between groups
      distances = as.matrix(dist(centroids[, c("PC1","PC2")]))
      colnames(distances) = centroids$group
      rownames(distances) = centroids$group
      sum_pca_pc1_pc2=rbind(sum_pca_pc1_pc2,c(paste0("set",set),distances["Ancestral","Fluctuating"],distances["Ancestral","Constant"],distances["Fluctuating","Constant"]))
      
      ##distant for pc1 to pc20
      # Calculate Euclidean distance between groups
      distances = as.matrix(dist(centroids[, c(paste0("PC",1:70)),method="euclidean"]))
      colnames(distances) = centroids$group
      rownames(distances) = centroids$group
      sum_pca_pc1_20=rbind(sum_pca_pc1_20,c(paste0("set",set),distances["Ancestral","Fluctuating"],distances["Ancestral","Constant"],distances["Fluctuating","Constant"]))

    }
  }
}
colnames(sum_pca_pc1)=c("set","anc_fluc","anc_cons","fluc_cons")
fwrite(sum_pca_pc1,"weighted_Euclidean_distant_PC1.txt",quote=FALSE,sep="\t")

colnames(sum_pca_pc1_pc2)=c("set","anc_fluc","anc_cons","fluc_cons")
fwrite(sum_pca_pc1_pc2,"weighted_Euclidean_distant_PC1_PC2.txt",quote=FALSE,sep="\t")

colnames(sum_pca_pc1_20)=c("set","anc_fluc","anc_cons","fluc_cons")
fwrite(sum_pca_pc1_20,"weighted_Euclidean_distant_PC1_to_PC20.txt",quote=FALSE,sep="\t")


