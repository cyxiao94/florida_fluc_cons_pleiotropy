library('dplyr')
library('tidyr')
library('data.table')
library('doMC')

source("source/workflow_functions.R")
source("source/load_packages.R")
#source("source/")


###FUNCTIONS
fit_GLM_ContinuousTime = function (afMatrix, rdMatrix, sampleData, vec, model.vars, poolCt = 100, 
                                   ncores) 
{
  registerDoMC(ncores)
  df = as.data.frame(sampleData[, colnames(sampleData) %in% 
                                  model.vars])
  colnames(df) <- model.vars
  formulaString = paste0(colnames(df), collapse = " + ")
  formulaString = paste0("cts ~ ", formulaString)
  cat("Model Formula is: \n", formulaString, "\n")
  Neff = calc_Neff(rdMatrix, poolCt)
  do.call(rbind, mclapply(1:nrow(afMatrix), function(ix) {
    if (ix%%10000 == 0) {
      cat("working on site ", ix, "\n")
    }
    cts = cbind(t(round(Neff[ix, ] * afMatrix[ix, ])), t(round(Neff[ix, 
    ] * (1 - afMatrix[ix, ]))))
    df$cts = cts
    df = as.data.frame(df)
    model = glm(formulaString, family = "quasibinomial", 
                data = df)
    cp = summary(model)$coefficients[-1, c(1, 4), drop = FALSE]
    results = c(cp[, 1], cp[, 2])
    names(results) = c(paste0("coef.", vec[1], "_", vec[length(vec)]), 
                       paste0("p.", vec[1], "_", vec[length(vec)]))
    return(results)
  }, mc.cores = ncores))
}

#calucluate N effective
calc_Neff = function (rd, poolCt) 
{
  ((poolCt * 2 * rd) - 1)/(poolCt * 2 + rd)
}


###Loop through treatments and run glm - save results as RData file
for (pop in c("fluc","cons")){
  load(paste0("./01genomic_data/",pop,"_HAF.Rdata"))
  vec = c(1:7)
  samps = samps %>% mutate(tpt = as.numeric(tpt))
  df.glm = sites
  eec = readRDS(paste0("./01genomic_data/",pop,"_cov.RDS"))
  afmat = as.data.frame(afmat)
  eec = as.data.frame(eec)
  samps = as.data.frame(samps)
  res = fit_GLM_ContinuousTime(afMatrix = afmat ,rdMatrix = eec, vec = vec, sampleData = samps, model.vars = 'tpt', poolCt=600, ncores = 16)
  df.glm = cbind(df.glm, res) #%>% dplyr::select(-chrom, -pos)
  names(df.glm) = c("chrom","pos",paste0('coef.1_7'),paste0('p.1_7'))
  save(df.glm, file = paste0("glm_",pop,".RData"))
}

##Getting significant sites of each treatment through time:
get.sig.sites = function(glm.file, rdata, comps, fdrThreshs, esThreshs){
  chroms = c('2L', '2R', '3L', '3R', 'X')
  load(rdata)
  af.shifts = get_af_shifts(afmat, samps, cage_set = NULL, comparisons = comps)
  #Get fdr 
  load(glm.file)
  FDR = get_glm_FDR(df.glm = df.glm, comparisons = comps)
  #load(glm.file)
  #df.glm = cbind(sites, df.glm)
  df.sig = get_sig_sites(df.glm, comparisons = comps, FDR  , af.shifts, fdrThreshs,
                         esThreshs)
  return(df.sig)
}


df.sig.fluc = get.sig.sites(glm.file = 'glm_fluc.RData', rdata = './01genomic_data/fluc_HAF.Rdata',
                         comps = c('1_7'), fdrThreshs= c(.2, .05, 0.01), esThreshs= c(0.005, 0.02, 0.02))
write.csv(df.sig.fluc, './df.sig.fluc.1_7.csv', row.names = FALSE)
df.sig.cons = get.sig.sites(glm.file = 'glm_cons.RData', rdata = './01genomic_data/cons_HAF.Rdata',
                            comps = c('1_7'), fdrThreshs= c(.2, .05, 0.01), esThreshs= c(0.005, 0.02, 0.02)) #here we keep all the information, also non-signifiacnt ones
write.csv(df.sig.cons, './df.sig.cons.1_7.csv', row.names = FALSE)




#neutral simulation----
##plot the expected AFC in our experimental scale
library(poolSeq)
dt<-data.table(
  af0 = seq(0,0.5, length=101)
)
#af0<-seq(0,1, length=101)

for (i in 1:100){ #100 simulations 
  
  af_t<-wf.traj(dt$af0, Ne = 300, t = 60, h = 0.5)
  dt[, paste0("AFC60_",i) := abs(af_t-af0)]
  
}

dt[, `:=` (
  mean_AFC = rowMeans(.SD),
  sd_AFC = apply(.SD, 1, sd),
  CI_lower = apply(.SD, 1, quantile, probs = 0.025),
  CI_upper = apply(.SD, 1, quantile, probs = 0.975)
), .SDcols = paste0("AFC60_",1:100)]


p<-ggplot(dt, aes(x=af0, y = mean_AFC))+
  geom_line()+
  geom_ribbon(aes(ymin = CI_lower,
                  ymax = CI_upper),
              alpha=0.2, fill = "blue")+
  #ylim(0,1)+
  geom_hline(yintercept = c(0.02), linetype="dashed",color="red")+
  labs(x = "Initial Frequecny", y = "Allele frequency change", title = "Expected AFC under Neutrality_Ne300_Gen60")+
  theme_minimal()
ggsave(file="plots/Expected_AFC_Neutral_Ne300_Gen60.pdf", plot=p)
#ggsave("")



