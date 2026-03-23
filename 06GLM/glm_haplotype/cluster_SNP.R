knitr::opts_chunk$set(echo = TRUE)

#library(dplyr)
source("source/workflow_functions.R")
source("source/helper_functions.R")
source("source/load_packages.R")

#df.sig<-read.csv("df.glm_sig.result.1_7.csv") 
df.sig_all <- read.csv("df.glm_sig_new_score.result.1_7.csv")
#df.sig <- subset(df.sig, chrom=="2L")
df.sig_all$comparison <- factor(df.sig_all$comparison) #make sure comparison is factor, required by score_wins function
load("01genomic_data/Florida_1stwave_HAF.RData")

winSize=500
winShift=100
maxClusterBreak=100
### score windows
cat("scoring windows..\n");flush.console();Sys.sleep(1)
df.wins <- suppressMessages(score_wins(df.sig_all ,sites,winSize,winShift))
df.sig.shuff <- suppressMessages(df.sig_all %>% group_by(comparison) %>% mutate(ix=sample(1:nrow(sites),n())) %>%
  ungroup() %>% mutate(chrom=sites$chrom[ix],pos=sites$pos[ix]))
df.wins.shuff <- suppressMessages(score_wins(df.sig.shuff,sites,winSize,winShift))
df.winfdr <- suppressMessages(get_win_fdr(df.wins,df.wins.shuff))

### cluster windows and merge by linkage
cat("finding initial clusters..\n");flush.console();Sys.sleep(1)
df.clust <- suppressMessages(cluster_wins(df.wins,df.winfdr,maxBreak=maxClusterBreak))
if(is.null(df.clust)){
  cat("**No clusters found** Ending early and returning results so far.\n");flush.console();Sys.sleep(1)
  results=list("sigSites"=df.sig_all,"wins"=df.wins,"params"=params)
  return(results)
  }
df.clust_all <- suppressMessages(df.clust %>% mutate(startPos=sites$pos[startSNP],endPos=sites$pos[endSNP]))
nClusters_init=nrow(df.clust)
cat(nClusters_init,"initial clusters found..\n");flush.console();Sys.sleep(1)

##save cluster file
fwrite(df.clust_all, file = "Clustered_SNPs_step1_poolseq.txt.gz", quote = FALSE, sep="\t")

maxSNPPairDist<-3000000
ncores=20
#snpFile="01genomic_data/snpFile_CHROM"
linkedClusterThresh<-0.03

df.clust_all <- as.data.frame(df.clust_all)
for (chr in c("2L","2R","3L","3R","X")){
  df.sig <- df.sig_all %>% filter(chrom==chr)
  snpFile <- paste0("01genomic_data/snpFile_",chr)
  snpSite <- fread(snpFile, select = 1, col.names = "pos")
  df.clust <- df.clust_all %>% filter(chrom==chr)

  nClusters_init=nrow(df.clust)

  if(nClusters_init>1){
    cat("merging linked clusters for chr", chr, " ..\n");flush.console();Sys.sleep(1)
    df.clust <- df.sig %>% filter(sigLevel>1) %>% 
      filter(pos %in% snpSite$pos) %>% ##only keep SNPs that exist in the haplotype sequencing file
      associate_snps_to_clusters(df.clust) %>%
      find_snp_pairs(maxDist=maxSNPPairDist, max_snps_per_cluster = 10000) %>%
      filter(pairType=="inter")  %>% ##only cal r2 for the SNPs across clusters
      calc_Rsq_for_snp_pairs(ncores=ncores,snpFile) %>%
      dplyr::select(-snp1.cl,-snp2.cl) %>%
      merge_linked_clusters(df.clust,df.sig,Rsq.thresh = linkedClusterThresh)
    }

  fwrite(df.clust, file = paste0("Clustered_SNPs_step2_founderhap_chr",chr,".txt.gz"), quote = FALSE, sep="\t")
}
