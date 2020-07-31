library(BayesMendel) # 2.1.6.1
packageVersion("BayesMendel")
library(foreach) # 1.4.4
library(doMC) # 1.3.5
library(parallel) # 3.5.2


source('../../penetrance_modification/brcapro_bcrat.R')

registerDoMC(cores=12)

script_num=1
prop_aj = 1
# get script_num from command line
args = commandArgs(trailingOnly=T)
if (length(args)>=1) {
  script_num = as.numeric(args[1])
  
  if (length(args)>=2) {
    prop_aj = as.numeric(args[2])/100
  }
}


load(paste0("../families/simFam_", prop_aj, "_", script_num, ".RData"))

# set proband's baseline breast cancer status to 0 so BRCAPRO doesn't calculate contralateral BC risk (will eventually exclude the probands affected at baseline)
fam$AffectedBreast[which(fam$ID==50)] = 0
fam$AgeBreast[which(fam$ID==50)] = fam$agecur[which(fam$ID==50)]
# change diagnosis age 1 to 2 so BRCAPRO doesn't treat it as unknown
fam$AgeOvary[which(fam$AffectedOvary==1 & fam$AgeOvary==1)] = 2


ids = unique(fam$FamID)
probands = fam[which(fam$ID==50),]
probands$T1 = probands$agecur
probands$T2 = probands$T1 + 5
probands$Race = 4
rr = relative.risk(probands)

bparams = brcaparams(age.by=5, age.to=94)
bparams$comprisk = 0*bparams$comprisk + 1

# get carrier probabilities, 5-year risk and 10-year risk
run.penmod = function(family, rr1, rr2, seed=1, impute=F) {
  set.seed(seed)
  out = penmod(family[, c("ID", "Gender", "MotherID", "FatherID", "ethnic", "Twins", "AffectedBreast", "AgeBreast", "AffectedOvary", "AgeOvary", "AgeBreastContralateral")], counselee.id=50, params=bparams,
               rr.gail=c(rr1, rr2), rr.pop=rr.ref,
               print=F, imputeAges=impute, imputeRelatives=impute, net=F)
  return(tryCatch(unlist(c(out@probs[1, c(1, 2, 3, 4)], out@predictions[1:2,2])),
                  error = function(e) rep(NA, 6)))
}

ind.start = 1
ind.end = length(ids)
results = foreach (i = ind.start:ind.end, .combine=rbind) %dopar% {
  c(ids[i], run.penmod(fam[which(fam$FamID==ids[i]), ], rr1=rr$RR_Star1[i], rr2=rr$RR_Star2[i]))
}

penmod.results = data.frame(results)
names(penmod.results) = c("FamID", "Prob.Carrier", "Prob.BRCA1", "Prob.BRCA2", "Prob.BRCA12",
                           "penmod5", "penmod10")


save(penmod.results, file=paste0("../model_results/run_penmod_sim_", prop_aj, "_", script_num, ".RData"))


