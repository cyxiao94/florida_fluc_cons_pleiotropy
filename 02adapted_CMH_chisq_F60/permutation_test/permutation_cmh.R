library(optparse)
library(data.table)
library(poolSeq)
library(ACER)

#functions
permu_estimate_ne<-function(mat_af,mat_cov){
  med_Ne=NULL
  input_rep <- 5
  NB_trials <- 100
  off <- 1
  nb_obs <- dim(mat_af)[1]
  #ne <- NULL
  for (i in 1:input_rep){
    ne=NULL;
    for(j in 1:NB_trials){
      set.seed(off)
      ind <- sample(x = 1:nb_obs, size = 1000)
      ne <- rbind(ne, data.frame(trial = j,
                                 ne = estimateNe(p0 = mat_af[ind, i], pt = mat_af[ind, i+5], 
                                                 cov0 = mat_cov[ind, i], covt = mat_cov[ind, i+5],
                                                 t = 60, ploidy=2, truncAF=0.05, method=c("P.planI"), 
                                                 poolSize=rep(600, times=2), Ncensus=1250)))
      off <- off+1
    }
    
    Ne <- median(ne$ne)
    ind <- which(is.na(ne$ne) | ne$ne<0)
    if(length(ind)>0){ne <- ne[-ind, ]}
    med_Ne<-c(med_Ne,round(median(ne$ne)))
    
  }
  return(med_Ne)
}

permu_cmh_test <- function(pop_af, pop_cov){
  ##autosome
  pop_af_auto <- pop_af[!grepl("^X\\.", rownames(pop_af)),]
  pop_cov_auto <- pop_cov[!grepl("^X\\.", rownames(pop_cov)),]
  pop_ne_auto<-permu_estimate_ne(pop_af_auto, pop_cov_auto)
  pop_cmh_auto<-adapted.cmh.test(freq = pop_af_auto, 
                                 coverage = pop_cov_auto, 
                                 Ne = pop_ne_auto, order=1, gen = c(0,60), repl = 1:5, poolSize = rep(600,10))
  ##X
  pop_af_x <- pop_af[grepl("^X\\.", rownames(pop_af)),]
  pop_cov_x <- pop_cov[grepl("^X\\.", rownames(pop_cov)),]
  pop_ne_x<-permu_estimate_ne(pop_af_x, pop_cov_x)
  pop_cmh_x<-adapted.cmh.test(freq = pop_af_x, 
                              coverage = pop_cov_x, 
                              Ne = pop_ne_auto, order=1, gen = c(0,60), repl = 1:5, poolSize = rep(600,10))
  
  pop_adj_p<-p.adjust(c(pop_cmh_auto,pop_cmh_x), method = "fdr")
  return(pop_adj_p)
}




option_list <- list(
  make_option(c("-p", "--permu"), type = "integer")
)


opt <- parse_args(OptionParser(option_list = option_list))

permu <- opt$permu

set.seed(as.integer(permu))

freqf<-readRDS("Anc_ConsF60_FlucF60_AF.rds")
covf<-readRDS("Anc_ConsF60_FlucF60_COV.rds")

anc_cols <- grepl("^Anc",colnames(freqf))
evo_cols <- grepl("^(Cons|Fluc)",colnames(freqf));

anc_freqf <- freqf[, anc_cols]
anc_covf <- covf[, anc_cols]

evo_freqf<-freqf[,evo_cols]
evo_covf<-covf[,evo_cols]

#permute evolved replicates
cols <- 1:10

perm_idx <- t(replicate(nrow(evo_freqf), sample(cols)))

af_perm <- evo_freqf
cov_perm <- evo_covf

for(i in seq_len(nrow(evo_freqf))){
  af_perm[i, cols]  <- evo_freqf[i, perm_idx[i, ]]
  cov_perm[i, cols] <- evo_covf[i, perm_idx[i, ]]
}

#cons
pop_af<-cbind(anc_freqf,af_perm[,c(1,3,5,7,9)])
pop_cov<-cbind(anc_covf,cov_perm[,c(1,3,5,7,9)])
cons_adj_p<-permu_cmh_test(pop_af, pop_cov)


#fluc
pop_af<-cbind(anc_freqf,af_perm[,c(2,4,6,8,10)])
pop_cov<-cbind(anc_covf,cov_perm[,c(2,4,6,8,10)])
fluc_adj_p<-permu_cmh_test(pop_af, pop_cov)

rm(pop_af, pop_cov)

#compare between pop1 and pop2
dt_cmh<-data.table(loci=rownames(anc_freqf),
                   p_cons=cons_adj_p,
                   p_fluc=fluc_adj_p)
dt_cmh<-dt_cmh[p_cons<0.05 | p_fluc < 0.05,]

#extract the AF information
dt_freqf <- as.data.table(freqf, keep.rownames = "loci")
dt_combined <- merge(dt_cmh, dt_freqf, by = "loci", all.x = TRUE)
dt_combined[, `:=`(
  AFC_Cons_R1 = Cons_F60.R1.freq-Anc_F0.R1.freq,
  AFC_Cons_R2 = Cons_F60.R2.freq-Anc_F0.R2.freq,
  AFC_Cons_R3 = Cons_F60.R3.freq-Anc_F0.R3.freq,
  AFC_Cons_R4 = Cons_F60.R4.freq-Anc_F0.R4.freq,
  AFC_Cons_R5 = Cons_F60.R5.freq-Anc_F0.R5.freq,
  AFC_Fluc_R1 = Fluc_F60.R1.freq-Anc_F0.R1.freq,
  AFC_Fluc_R2 = Fluc_F60.R2.freq-Anc_F0.R2.freq,
  AFC_Fluc_R3 = Fluc_F60.R3.freq-Anc_F0.R3.freq,
  AFC_Fluc_R4 = Fluc_F60.R4.freq-Anc_F0.R4.freq,
  AFC_Fluc_R5 = Fluc_F60.R5.freq-Anc_F0.R5.freq
)]

dt_combined[, `:=`(
  median_AFC_Cons = apply(.SD, 1, median, na.rm =TRUE)
), .SDcols = paste0("AFC_Cons_R",1:5)]
dt_combined[, `:=`(
  median_AFC_Fluc = apply(.SD, 1, median, na.rm =TRUE)
), .SDcols = paste0("AFC_Fluc_R",1:5)]

#subset the columns
sel_cols <- c("loci", "p_cons", "p_fluc", "median_AFC_Cons", "median_AFC_Fluc")
dt_sub <- dt_combined[, ..sel_cols]
rm(dt_combined, dt_freqf)

fwrite(dt_sub, 
       file = paste0("permutation_result/Permutation", permu, "_CMH_result.txt.gz"),
       quote = FALSE, sep="\t")

#write summary information
num_sig_cons<-sum(dt_sub$p_cons<0.05)
num_sig_fluc<-sum(dt_sub$p_fluc<0.05)
num_sig_total<-nrow(dt_sub[p_cons<0.05 | p_fluc<0.05,])
num_sig_share<-nrow(dt_sub[p_cons<0.05 & p_fluc<0.05,])
percent_share<-num_sig_share/num_sig_total

#AFC in the Cons population
#SNPs significant in the constant population
#polarize to increase in the cons
dt_sub[, polarize_median_AFC_Cons := abs(median_AFC_Cons)]
dt_sub[, polarize_median_AFC_Fluc := fifelse(median_AFC_Cons>0, median_AFC_Fluc, (-1)*median_AFC_Fluc)]
#average AFC across the genome
meanAFC_cons_sig_in_cons <- dt_sub[p_cons <0.05, mean(polarize_median_AFC_Cons, na.rm =TRUE)]
meanAFC_cons_sig_in_fluc <- dt_sub[p_cons <0.05, mean(polarize_median_AFC_Fluc, na.rm =TRUE)]

#AFC in the Fluc population
#SNPs significant in the constant population
#polarize to increase in the fluc
dt_sub[, polarize_median_AFC_Cons := fifelse(median_AFC_Fluc>0, median_AFC_Cons, (-1)*median_AFC_Cons)]
dt_sub[, polarize_median_AFC_Fluc := abs(median_AFC_Fluc)]
#average AFC across the genome
meanAFC_fluc_sig_in_cons <- dt_sub[p_fluc <0.05, mean(polarize_median_AFC_Cons, na.rm =TRUE)]
meanAFC_fluc_sig_in_fluc <- dt_sub[p_fluc <0.05, mean(polarize_median_AFC_Fluc, na.rm =TRUE)]

cat(
  paste(
    permu,
    num_sig_cons,
    num_sig_fluc,
    num_sig_total,
    num_sig_share,
    percent_share,
    meanAFC_cons_sig_in_cons,
    meanAFC_cons_sig_in_fluc,
    meanAFC_fluc_sig_in_fluc,
    meanAFC_fluc_sig_in_cons,
    sep="\t"),
  "\n",
  file = paste0("permutation_result/Permutation", permu, "_summary_result.txt")
)

