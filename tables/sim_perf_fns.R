
source("../../../ModelEvaluation/ModelEvaluation.R")

library(BayesMendel)
library(BCRA)
library(mc2d)
library(data.table)
source("../sim_code/pheno.gen.R")
load("../../penetrance_modification/rr_pop_nhis2015.RData")



##### functions for simulating times to breast cancer based on BRCAPRO+BCRAT (M) penetrances

# function to convert hazard to penetrance
# haz.bc.mat and haz.death.mat are matrices where each row corresponds to a proband
haz2dens = function(haz.bc.mat, haz.death.mat) {
  
  dens = 0*haz.bc.mat
  surv.all <- 1 - haz.bc.mat[, 1] - haz.death.mat[, 1]
  dens[, 1] = haz.bc.mat[, 1]
  
  
  for(i in 2:ncol(dens)){
    dens[, i] = haz.bc.mat[, i]*surv.all
    surv.all <- surv.all * (1 - haz.bc.mat[, i] - haz.death.mat[, i])
  }
  return(dens)
}

# function to generate times to breast cancer based on BRCAPRO+BCRAT (M) penetrances
generate.outcomes = function(probands, penetrance = penet.brca.crude, haz.death = haz.death.crude$femaleBC, times=NULL) {
  
  probands$T1 = probands$agecur
  probands$T2 = probands$T1 + 5
  probands$Race = 4
  rr = relative.risk(probands)
  rr1 = rr$RR_Star1/rr.ref$noncarrier[1]
  rr2 = rr$RR_Star2/rr.ref$noncarrier[2]
  
  ind0 = which(probands$BRCA1==0 & probands$BRCA2==0)
  ind1 = which(probands$BRCA1==1 & probands$BRCA2==0)
  ind2 = which(probands$BRCA1==0 & probands$BRCA2==1)
  ind12 = which(probands$BRCA1==1 & probands$BRCA2==1)
  
  penet.mat = matrix(NA, nrow=nrow(probands), ncol=94)
  penet.mat = data.frame(cbind(probands$FamID, penet.mat))
  names(penet.mat)[1] = "FamID"
  names(penet.mat)[2:95] = paste0("p", 1:94)
  
  ### get hazards for carriers 
  haz10 = dens2haz.cs(penet.brca.crude$fFX[, "B10"], haz.death)
  haz01 = dens2haz.cs(penet.brca.crude$fFX[, "B01"], haz.death)
  haz11 = dens2haz.cs(penet.brca.crude$fFX[, "B11"], haz.death)
  
  ### get hazards for non-carriers
  # baseline non-carrier hazard
  haz0 = dens2haz.cs(penet.brca.crude$fFX[, "B00"], haz.death)
  haz.mat = matrix(rep(haz0, nrow(probands)), nrow=nrow(probands), ncol=94, byrow=T)
  haz.mat = data.frame(cbind(probands$FamID, haz.mat))
  names(haz.mat)[1] = "FamID"
  names(haz.mat)[2:95] = paste0("h", 1:94)
  ## modify hazard using BCRAT relative risk
  # 1-(1-haz0)^rr
  rr1.mat = matrix(rep(rr1, length(1:49)), nrow=nrow(probands))
  rr2.mat = matrix(rep(rr2, length(50:94)), nrow=nrow(probands))
  rr.mat = cbind(haz.mat$FamID, rr1.mat, rr2.mat)
  haz.mat.mod = haz.mat
  haz.mat.mod[, 2:95] = 1-(1-haz.mat[, 2:95])^rr.mat[, 2:95]
  #i = 200
  #sum(haz.mat.mod[i, 2:95] - c((1-(1-haz0)^rr1[i])[1:49], (1-(1-haz0)^rr2[i])[50:94]))
  
  ### combine hazards for carriers and non-carriers
  haz.all = haz.mat.mod
  haz.all[ind1, 2:95] = matrix(rep(haz10, length(ind1)), nrow=length(ind1), ncol=94, byrow=T)
  haz.all[ind2, 2:95] = matrix(rep(haz01, length(ind2)), nrow=length(ind2), ncol=94, byrow=T)
  haz.all[ind12, 2:95] = matrix(rep(haz11, length(ind12)), nrow=length(ind12), ncol=94, byrow=T)
  
  ### get death hazard
  haz.death.mat = matrix(rep(haz.death, nrow(probands)), nrow=nrow(probands), ncol=94, byrow=T)
  
  ### set breast cancer and death hazards to 0 for ages <= current age
  unique.ages = unique(probands$agecur) 
  for (i in 1:length(unique.ages)) {
    age.i = unique.ages[i]
    ind.i = which(probands$agecur==age.i)
    haz.all[ind.i, 2:(age.i+1)] = 0
    haz.death.mat[ind.i, 1:age.i] = 0
  }
  
  ### convert hazards to penetrances
  penet.mat[, 2:95] = haz2dens(haz.all[, 2:95], haz.death.mat)
  penet.mat = cbind(penet.mat, 1-rowSums(penet.mat[, 2:95]))
  
  ### sample age to breast cancer 
  AgeBreast.lt = rmultinomial(nrow(probands), 1, as.matrix(penet.mat[, 2:96]))
  probands$AgeBreast.lt = NA
  for (i in 1:95) {
    probands$AgeBreast.lt[which(AgeBreast.lt[, i]==1)] = i
  }
  
  probands$AffectedBreast.lt = probands$AgeBreast.lt<95
  
  if (!is.null(times)) {
    for (i in times) {
      probands[, paste0("BC.", i)] = probands$AffectedBreast.lt==1 & probands$AgeBreast.lt <= probands$agecur+i
    }
  }
  
  return(probands)
  
}





###### helper functions for BRCAPRO+BCRAT (E2)


# function to format input for BRCAPRO+BCRAT (E2)
convert.to.long = function(df, times=1:5) {
  df.list = vector("list", length = length(times))
  for (i in times) {
    df[, paste0("BC.", i)] = df$AffectedBreast.lt==1 & df$AgeBreast.lt <= df$agecur+i
    
    df.list[[i]] = df[, c("id", paste0("BC.", i), paste0("gail", i), paste0("brcapro", i))]
    df.list[[i]][, "BC"] = as.numeric(df.list[[i]][, paste0("BC.", i)])
    df.list[[i]][, "brcapro"] = df.list[[i]][, paste0("brcapro", i)]
    df.list[[i]][, "gail"] = df.list[[i]][, paste0("gail", i)]
    df.list[[i]][, "time"] = i
    df.list[[i]][, c(paste0("BC.", i), paste0("gail", i), paste0("brcapro", i))] = NULL
  }
  
  df.long = rbindlist(df.list)
  return(df.long)
}


expit = function(x) 1/(1+exp(-x))

# function to get GEE predictions
predict.gee = function(fit.geeglm, df) {
  df$sqrt.brcapro = sqrt(df$brcapro)
  df$sqrt.gail = sqrt(df$gail)
  df$sqrt.brcapro.gail = sqrt(df$brcapro) * sqrt(df$gail)
  df$sqrt.brcapro.time = sqrt(df$brcapro) * df$time
  df$sqrt.gail.time = sqrt(df$gail) * df$time
  df$sqrt.brcapro.gail.time = sqrt(df$brcapro) * sqrt(df$gail) * df$time
  pred.gee = expit(as.matrix(cbind(rep(1, nrow(df)), df[, c("sqrt.brcapro", "sqrt.gail", "time", "sqrt.brcapro.gail", "sqrt.brcapro.time", "sqrt.gail.time", "sqrt.brcapro.gail.time")])) %*% matrix(coef(fit.geeglm), nrow=length(coef(fit.geeglm))))
  return(pred.gee)
}




##### NCCN criteria

# function for checking NCCN criteria
nccn = function(probands, fam) {
  firstdeg = c("Mother", "Father", "Sister", "Brother", "Daughter", "Son")
  maternal = paste0("Maternal", c("Grandmother", "Grandfather", "Aunt", "Uncle"))
  paternal = paste0("Paternal", c("Grandmother", "Grandfather", "Aunt", "Uncle"))
  
  oc = data.frame(table(fam$FamID[which(fam$AffectedOvary==1)]))
  male.bc = data.frame(table(fam$FamID[which(fam$AffectedBreast==1 & fam$Gender==1)]))
  early.bc = data.frame(table(fam$FamID[which(fam$AffectedBreast==1 & fam$AgeBreast<=45)]))
  bc = data.frame(table(fam$FamID[which(fam$AffectedBreast==1)]))
  bc.mat = data.frame(table(fam$FamID[which(fam$AffectedBreast==1 & fam$relationship %in% c(firstdeg, maternal))]))
  bc.pat = data.frame(table(fam$FamID[which(fam$AffectedBreast==1 & fam$relationship %in% c(firstdeg, paternal))]))
  bc.mat.under50 = data.frame(table(fam$FamID[which(fam$AffectedBreast==1 & fam$AgeBreast<=50 & fam$relationship %in% c(firstdeg, maternal))]))
  bc.pat.under50 = data.frame(table(fam$FamID[which(fam$AffectedBreast==1 & fam$AgeBreast<=50 & fam$relationship %in% c(firstdeg, paternal))]))
  
  probands$nccn = probands$FamID %in% as.numeric(as.character(c(oc$Var1[which(oc$Freq>0)], male.bc$Var1[which(male.bc$Freq>0)]))) |
    (probands$FamID %in%  as.numeric(as.character(bc.mat$Var1[which(bc.mat$Freq>1)])) & probands$FamID %in%  as.numeric(as.character(bc.mat.under50$Var1[which(bc.mat.under50$Freq>0)]))) | 
    (probands$FamID %in%  as.numeric(as.character(bc.pat$Var1[which(bc.pat$Freq>1)])) & probands$FamID %in%  as.numeric(as.character(bc.pat.under50$Var1[which(bc.pat.under50$Freq>0)]))) | 
    probands$FamID %in%  as.numeric(as.character(bc.mat$Var1[which(bc.mat$Freq>2)])) |
    probands$FamID %in%  as.numeric(as.character(bc.pat$Var1[which(bc.pat$Freq>2)])) 
  
  return(probands)
}




##### performance measures

# root brier score (lower is better)
calc.sqrt.brier = function(scores, outcomes, weights=NULL) {
  if (is.null(weights)) {
    weights = rep(1, length(scores))
  }
  ind = which(outcomes>=0 & scores>=0)
  return(sqrt(mean((outcomes[ind]-scores[ind])^2 * weights[ind])))
}

# log scoring rule (higher is better)
calc.lsr = function(scores, outcomes, weights=NULL) {
  if (is.null(weights)) {
    weights = rep(1, length(scores))
  }
  ind = which(outcomes>=0 & scores>0)
  return(sum(outcomes*log(scores)* weights[ind] + (1-outcomes)*(log(1 - scores)* weights[ind])))
}

# standardized net benefit
# r: risk threshold
# p: outcome prevalence (if NULL, then estimated from observed outcomes)
calc.snb = function(scores, outcomes, weights=NULL, r=0.0167, p=NULL) {
  if (is.null(weights)) {
    weights = rep(1, length(scores))
  }
  if (is.null(p)) {
    p = sum(outcomes*weights)/length(outcomes)
  }
  tpr = sum((outcomes==1 & scores>=r)*weights )/sum( (outcomes==1)*weights )
  fpr = sum((outcomes==0 & scores>=r)*weights)/sum( (outcomes==0)*weights )
  SNB = tpr - fpr*(r/(1-r))*(1-p)/p
  return(SNB)
}


# IPCW weights based on Kaplan-Meier stratified by center and carrier probability>0.1
km.weights.2 = function(dataset, t) {
  # estimate censoring distribution using Kaplan-Meier
  km = summary(survfit(Surv(Time, y==0) ~ CENTER, data=dataset[which(dataset$Prob.Carrier>0.1),]))
  km2 = summary(survfit(Surv(Time, y==0) ~ CENTER, data=dataset[which(dataset$Prob.Carrier<=0.1),]))
  
  
  km.strata = gsub("CENTER=", "", km$strata)
  unique.strata = unique(km.strata)
  ind.1 = which(km.strata == unique.strata[1])
  survest = stepfun(km$time[ind.1], c(1, km$surv[ind.1]))
  # set each person's weight to their inverse probability of not being censored by min(Time, t)
  weights = 1/survest(pmin(dataset$Time, t))
  for (i in 2:length(unique.strata)) {
    ind.i = which(km.strata == unique.strata[i])
    survest.i = stepfun(km$time[ind.i], c(1, km$surv[ind.i]))
    weights.i = 1/survest.i(pmin(dataset$Time, t))
    weights[which(dataset$CENTER==km.centers[i])] = weights.i[which(dataset$CENTER==km.centers[i])]
  }
  
  # if censored before time t, set weight to 0
  weights[which(dataset$y==0 & dataset$Time<t | is.na(dataset$Time))] = 0
  
  
  
  km2.strata = gsub("CENTER=", "", km2$strata)
  unique.strata = unique(km2.strata)
  ind.1 = which(km2.strata == unique.strata[1])
  survest = stepfun(km2$time[ind.1], c(1, km2$surv[ind.1]))
  weights2 = 1/survest(pmin(dataset$Time, t))
  
  for (i in 2:length(unique.strata)) {
    ind.i = which(km2.strata == unique.strata[i])
    survest.i = stepfun(km2$time[ind.i], c(1, km2$surv[ind.i]))
    weights2.i = 1/survest.i(pmin(dataset$Time, t))
    weights2[which(dataset$CENTER==km.centers[i])] = weights2.i[which(dataset$CENTER==km.centers[i])]
  }
  
  # if censored before time t, set weight to 0
  weights2[which(dataset$y==0 & dataset$Time<t | is.na(dataset$Time))] = 0
  
  weights[which(dataset$Prob.Carrier<=0.1)] = weights2[which(dataset$Prob.Carrier<=0.1)]
  
  return(weights)
}

# IPCW weights based on Kaplan-Meier stratified by center
km.weights.3 = function(dataset, t) {
  # estimate censoring distribution using Kaplan-Meier
  km = summary(survfit(Surv(Time, y==0) ~ CENTER, data=dataset))
  km.strata = gsub("CENTER=", "", km$strata)
  unique.strata = unique(km.strata)
  ind.1 = which(km.strata == unique.strata[1])
  survest = stepfun(km$time[ind.1], c(1, km$surv[ind.1]))
  # set each person's weight to their inverse probability of not being censored by min(Time, t)
  weights = 1/survest(pmin(dataset$Time, t))
  for (i in 2:length(unique.strata)) {
    ind.i = which(km.strata == unique.strata[i])
    survest.i = stepfun(km$time[ind.i], c(1, km$surv[ind.i]))
    weights.i = 1/survest.i(pmin(dataset$Time, t))
    weights[which(dataset$CENTER==unique.strata[i])] = weights.i[which(dataset$CENTER==unique.strata[i])]
  }
  
  # if censored before time t, set weight to 0
  weights[which(dataset$y==0 & dataset$Time<t | is.na(dataset$Time))] = 0
  return(weights)
}


perf.boot.2 = function(score.matrix, outcomes, cens.dist="none", dataset=NULL, covars=NULL, t=NULL, 
                       model.names=NULL, nboot=200, seed=1) {
  set.seed(seed)
  score.matrix = as.matrix(score.matrix)
  n.models = ncol(score.matrix)
  
  # matrices for storing bootstrap results
  oe = matrix(ncol=n.models, nrow=nboot+1)
  auc = matrix(ncol=n.models, nrow=nboot+1)
  brier = matrix(ncol=n.models, nrow=nboot+1)
  
  # get weights if using IPCW
  weights = NULL
  if (cens.dist=="km") {
    if (is.null(dataset) | is.null(t)) {
      stop("dataset and t cannot be NULL.")
    }
    weights = km.weights(dataset, t)
  } else if (cens.dist=="cox") {
    if (is.null(dataset) | is.null(t) | is.null(covars)) {
      stop("dataset, covars, and t cannot be NULL.")
    }
    weights = cox.weights(dataset, covars, t)
  }
  
  # evaluate performance in observed sample
  oe[1, ] = apply(score.matrix, 2, function(x) calc.oe(x, outcomes, weights))
  auc[1, ] = apply(score.matrix, 2, function(x) calc.auc(x, outcomes, weights))
  brier[1, ] = apply(score.matrix, 2, function(x) calc.sqrt.brier(x, outcomes, weights))
  
  for (i in 2:(nboot+1)) {
    # generate bootstrap sample
    ind.boot = sample(1:length(outcomes), length(outcomes), replace=T)
    score.matrix.boot = as.matrix(score.matrix[ind.boot, ])
    outcomes.boot = outcomes[ind.boot]
    # re-calculate weights for bootstrap sample if using IPCW
    weights.boot = NULL
    if (cens.dist=="km") {
      weights.boot = km.weights(dataset[ind.boot,], t)
    } else if (cens.dist=="cox") {
      weights.boot = cox.weights(dataset[ind.boot,], covars, t)
    }
    # evaluate performance in bootstrap sample
    oe[i, ] = apply(score.matrix.boot, 2, function(x) calc.oe(x, outcomes.boot, weights.boot))
    auc[i, ] = apply(score.matrix.boot, 2, function(x) calc.auc(x, outcomes.boot, weights.boot))
    brier[i, ] = apply(score.matrix.boot, 2, function(x) calc.sqrt.brier(x, outcomes.boot, weights.boot))
  }
  
  # calculate 95% CIs for performance measures using 2.5th and 97.5th percentiles of bootstrap estimates
  oe.summary = apply(oe, 2, function(x) c(x[1], quantile(x[-1], 0.025), quantile(x[-1], 0.975)))
  auc.summary = apply(auc, 2, function(x) c(x[1], quantile(x[-1], 0.025), quantile(x[-1], 0.975)))
  brier.summary = apply(brier, 2, function(x) c(x[1], quantile(x[-1], 0.025), quantile(x[-1], 0.975)))
  perf.matrix = as.matrix(rbind(oe.summary, auc.summary, brier.summary))
  rownames(perf.matrix) = c("oe", "oe.lo", "oe.hi",
                            "auc", "auc.lo", "auc.hi",
                            "brier", "brier.lo", "brier.hi")
  colnames(perf.matrix) = model.names
  if (is.null(colnames(perf.matrix)) & !is.null(colnames(score.matrix))) {
    colnames(perf.matrix) = colnames(score.matrix)
  } else if (is.null(colnames(perf.matrix))) {
    colnames(perf.matrix) = paste("model", 1:n.models)
  }
  return(perf.matrix)
}


perf.boot.3 = function(score.matrix, outcomes, cens.dist="none", dataset=NULL, covars=NULL, t=NULL, 
                       model.names=NULL, nboot=200, seed=1) {
  set.seed(seed)
  score.matrix = as.matrix(score.matrix)
  n.models = ncol(score.matrix)
  
  # matrices for storing bootstrap results
  oe = matrix(ncol=n.models, nrow=nboot+1)
  auc = matrix(ncol=n.models, nrow=nboot+1)
  brier = matrix(ncol=n.models, nrow=nboot+1)
  
  # get weights if using IPCW
  weights = NULL
  if (cens.dist=="km") {
    if (is.null(dataset) | is.null(t)) {
      stop("dataset and t cannot be NULL.")
    }
    weights = km.weights.3(dataset, t)
  } else if (cens.dist=="cox") {
    if (is.null(dataset) | is.null(t) | is.null(covars)) {
      stop("dataset, covars, and t cannot be NULL.")
    }
    weights = cox.weights(dataset, covars, t)
  }
  
  # evaluate performance in observed sample
  oe[1, ] = apply(score.matrix, 2, function(x) calc.oe(x, outcomes, weights))
  auc[1, ] = apply(score.matrix, 2, function(x) calc.auc(x, outcomes, weights))
  brier[1, ] = apply(score.matrix, 2, function(x) calc.sqrt.brier(x, outcomes, weights))
  
  for (i in 2:(nboot+1)) {
    # generate bootstrap sample
    ind.boot = sample(1:length(outcomes), length(outcomes), replace=T)
    score.matrix.boot = as.matrix(score.matrix[ind.boot, ])
    outcomes.boot = outcomes[ind.boot]
    # re-calculate weights for bootstrap sample if using IPCW
    weights.boot = NULL
    if (cens.dist=="km") {
      weights.boot = km.weights(dataset[ind.boot,], t)
    } else if (cens.dist=="cox") {
      weights.boot = cox.weights(dataset[ind.boot,], covars, t)
    }
    # evaluate performance in bootstrap sample
    oe[i, ] = apply(score.matrix.boot, 2, function(x) calc.oe(x, outcomes.boot, weights.boot))
    auc[i, ] = apply(score.matrix.boot, 2, function(x) calc.auc(x, outcomes.boot, weights.boot))
    brier[i, ] = apply(score.matrix.boot, 2, function(x) calc.sqrt.brier(x, outcomes.boot, weights.boot))
  }
  
  # calculate 95% CIs for performance measures using 2.5th and 97.5th percentiles of bootstrap estimates
  oe.summary = apply(oe, 2, function(x) c(x[1], quantile(x[-1], 0.025), quantile(x[-1], 0.975)))
  auc.summary = apply(auc, 2, function(x) c(x[1], quantile(x[-1], 0.025), quantile(x[-1], 0.975)))
  brier.summary = apply(brier, 2, function(x) c(x[1], quantile(x[-1], 0.025), quantile(x[-1], 0.975)))
  perf.matrix = as.matrix(rbind(oe.summary, auc.summary, brier.summary))
  rownames(perf.matrix) = c("oe", "oe.lo", "oe.hi",
                            "auc", "auc.lo", "auc.hi",
                            "brier", "brier.lo", "brier.hi")
  colnames(perf.matrix) = model.names
  if (is.null(colnames(perf.matrix)) & !is.null(colnames(score.matrix))) {
    colnames(perf.matrix) = colnames(score.matrix)
  } else if (is.null(colnames(perf.matrix))) {
    colnames(perf.matrix) = paste("model", 1:n.models)
  }
  return(perf.matrix)
}



# modified version of perf.boot() from ModelEvaluation.R that includes standardized net benefit and log scoring rule, as well as relative improvement in brier score / log scoring rule compared to a reference model 
### input
# score.matrix: matrix of predicted probabilities where each column corresponds to a different model
# outcomes: vector of outcomes (1 for case, 0 for control)
# cens.dist: character string specifying IPCW type 
# "none" if there is no censoring/IPCW should not be used
# "km" if Kaplan-Meier (if cens.dist="km", then dataset and t cannot be NULL)
# "cox" if Cox (if cens.dist="cox", then dataset, covars, and t cannot be NULL)
# "km3" if Kaplan-Meier stratified by CENTER (for CGN data)
# dataset: if cens.dist="km" or "cox", dataframe with observed event statuses and times 
# must have columns Time (minimum of censoring time and event time) and y (event status: 1 if event observed, 0 otherwise)
# if cens.dist="cox", must also have columns corresponding to the names in covars
# covars: if cens.dist="cox", vector of names of covariates to be included in Cox model for censoring
# t: if cens.dist="km" or "cox", cutoff time used for dichotomizing the time-to-event outcome
# model.names: vector of model names corresponding to the columns of score.matrix 
# nboot: number of bootstrap samples
# seed: seed for random number generator
# ref.col: column in score.matrix corresponding to reference model for calculating difference in brier score and log scoring rule
# r: risk threshold for standardized net benefit
# p: outcome prevalence for standardized net benefit (if NULL, then estimated from observed outcomes)
### output
# list where the first element is the performance matrix output of perf.boot(), the second element is the matrix of bootstrap O/E values, the third element is the matrix of bootstrap AUC values, the fourth element is the matrix of brier score bootstrap values 
perf.boot.2 = function(score.matrix, outcomes, cens.dist="none", dataset=NULL, covars=NULL, t=NULL, 
                       model.names=NULL, nboot=200, seed=1, ref.col=3, print=T) {
  set.seed(seed)
  score.matrix = as.matrix(score.matrix)
  n.models = ncol(score.matrix)
  
  # matrices for storing bootstrap results
  oe = matrix(ncol=n.models, nrow=nboot+1)
  auc = oe
  brier = oe
  brier.change = oe
  nb = oe
  lsr = oe
  lsr.change = oe
  
  # get weights if using IPCW
  weights = NULL
  if (cens.dist=="km") {
    if (is.null(dataset) | is.null(t)) {
      stop("dataset and t cannot be NULL.")
    }
    weights = km.weights(dataset, t)
  } else if (cens.dist=="cox") {
    if (is.null(dataset) | is.null(t) | is.null(covars)) {
      stop("dataset, covars, and t cannot be NULL.")
    }
    weights = cox.weights(dataset, covars, t)
  } else if (cens.dist=="km3") {
    if (is.null(dataset) | is.null(t)) {
      stop("dataset and t cannot be NULL.")
    }
    weights = km.weights.3(dataset, t)
  }
  
  # evaluate performance in observed sample
  oe[1, ] = apply(score.matrix, 2, function(x) calc.oe(x, outcomes, weights))
  auc[1, ] = apply(score.matrix, 2, function(x) calc.auc(x, outcomes, weights))
  brier[1, ] = apply(score.matrix, 2, function(x) calc.brier(x, outcomes, weights))
  brier.ref = brier[1, ref.col]
  brier.change[1, ] = 100*(brier.ref - brier[1, ])/brier.ref
  nb[1, ] = apply(score.matrix, 2, function(x) calc.snb(x, outcomes, weights))
  lsr[1,] = apply(score.matrix, 2, function(x) calc.lsr(x, outcomes, weights))
  lsr.ref = lsr[1, ref.col]
  lsr.change[1, ] = 100*(lsr[1, ] - lsr.ref)/abs(lsr.ref)
  
  for (i in 2:(nboot+1)) {
    if (print) {
      print(i)
    }
    
    # generate bootstrap sample
    ind.boot = sample(1:length(outcomes), length(outcomes), replace=T)
    score.matrix.boot = as.matrix(score.matrix[ind.boot, ])
    outcomes.boot = outcomes[ind.boot]
    # re-calculate weights for bootstrap sample if using IPCW
    weights.boot = NULL
    if (cens.dist=="km") {
      weights.boot = km.weights(dataset[ind.boot,], t)
    } else if (cens.dist=="cox") {
      weights.boot = cox.weights(dataset[ind.boot,], covars, t)
    } else if (cens.dist=="km3") {
      weights.boot = km.weights.3(dataset[ind.boot,], t)
    }
    # evaluate performance in bootstrap sample
    oe[i, ] = apply(score.matrix.boot, 2, function(x) calc.oe(x, outcomes.boot, weights.boot))
    auc[i, ] = apply(score.matrix.boot, 2, function(x) calc.auc(x, outcomes.boot, weights.boot))
    brier[i, ] = apply(score.matrix.boot, 2, function(x) calc.brier(x, outcomes.boot, weights.boot))
    brier.ref = brier[i, ref.col]
    brier.change[i, ] = 100*(brier.ref - brier[i, ])/brier.ref
    nb[i, ] = apply(score.matrix.boot, 2, function(x) calc.snb(x, outcomes.boot, weights.boot))
    lsr[i,] = apply(score.matrix.boot, 2, function(x) calc.lsr(x, outcomes.boot, weights.boot))
    lsr.ref = lsr[i, ref.col]
    lsr.change[i, ] = 100*(lsr[i, ] - lsr.ref)/abs(lsr.ref)
    
  }
  
  # calculate 95% CIs for performance measures using 2.5th and 97.5th percentiles of bootstrap estimates
  oe.summary = apply(oe, 2, function(x) c(x[1], quantile(x[-1], 0.025), quantile(x[-1], 0.975)))
  auc.summary = apply(auc, 2, function(x) c(x[1], quantile(x[-1], 0.025), quantile(x[-1], 0.975)))
  brier.summary = apply(brier, 2, function(x) c(x[1], quantile(x[-1], 0.025), quantile(x[-1], 0.975)))
  brier.change.summary = apply(brier.change, 2, function(x) c(x[1], quantile(x[-1], 0.025), quantile(x[-1], 0.975)))
  nb.summary = apply(nb, 2, function(x) c(x[1], quantile(x[-1], 0.025), quantile(x[-1], 0.975)))
  lsr.summary = apply(lsr, 2, function(x) c(x[1], quantile(x[-1], 0.025), quantile(x[-1], 0.975)))
  lsr.change.summary = apply(lsr.change, 2, function(x) c(x[1], quantile(x[-1], 0.025), quantile(x[-1], 0.975)))
  
  perf.matrix = as.matrix(rbind(oe.summary, auc.summary, brier.summary, brier.change.summary, nb.summary, lsr.summary, lsr.change.summary))
  rownames(perf.matrix) = c("oe", "oe.lo", "oe.hi",
                            "auc", "auc.lo", "auc.hi",
                            "brier", "brier.lo", "brier.hi",
                            "brier.c", "brier.c.lo", "brier.c.hi",
                            "snb", "snb.lo", "snb.hi",
                            "lsr", "lsr.lo", "lsr.hi",
                            "lsr.c", "lsr.c.lo", "lsr.c.hi")
  colnames(perf.matrix) = model.names
  if (is.null(colnames(perf.matrix)) & !is.null(colnames(score.matrix))) {
    colnames(perf.matrix) = colnames(score.matrix)
  } else if (is.null(colnames(perf.matrix))) {
    colnames(perf.matrix) = paste("model", 1:n.models)
  }
  return(list(perf.matrix, oe, auc, brier, nb, lsr))
}


# modified version of format.perf.table() from ModelEvaluation.R that formats first element of output of perf.boot.2()
format.perf.table.2 = function(perf.matrix, digits=c(2, 2, 2, 2, 2), model.names=NULL) {
  oe.matrix = matrix(round(perf.matrix[1:3, ], digits[1]), ncol=ncol(perf.matrix))
  auc.matrix = matrix(round(perf.matrix[4:6, ], digits[2]), ncol=ncol(perf.matrix))
  brier.matrix = matrix(round(perf.matrix[10:12, ], digits[3]), ncol=ncol(perf.matrix))
  nb.matrix = matrix(round(perf.matrix[13:15, ], digits[4]), ncol=ncol(perf.matrix))
  lsr.matrix = matrix(round(perf.matrix[19:21, ], digits[5]), ncol=ncol(perf.matrix))
  
  oe.matrix = formatC(oe.matrix, digits[1], format="f")
  auc.matrix = formatC(auc.matrix, digits[2], format="f")
  brier.matrix = formatC(brier.matrix, digits[3], format="f")
  nb.matrix = formatC(nb.matrix, digits[4], format="f")
  lsr.matrix = formatC(lsr.matrix, digits[5], format="f")
  
  formatted.matrix = cbind(apply(oe.matrix, 2, function(x) paste0(x[1], " (", x[2], ", ",  x[3], ")") ),
                           apply(auc.matrix, 2, function(x) paste0(x[1], " (", x[2], ", ",  x[3], ")") ),
                           apply(nb.matrix, 2, function(x) paste0(x[1], " (", x[2], ", ",  x[3], ")") ),
                           apply(brier.matrix, 2, function(x) paste0(x[1], " (", x[2], ", ",  x[3], ")") ),
                           apply(lsr.matrix, 2, function(x) paste0(x[1], " (", x[2], ", ",  x[3], ")") ))
  colnames(formatted.matrix) = c("O/E", "AUC", "SNB", "diff.BS", "diff.LS")
  rownames(formatted.matrix) = model.names
  return(formatted.matrix)
}


