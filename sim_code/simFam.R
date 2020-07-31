### simulate baseline covariates
# simulate non-AJ and AJ families based on family structures sampled from CGN
# sample BCRAT covariates from CGN

## simulate under penetrance modification
# carrier probabilities are calculated in the same way as in BRCAPRO, based on net penetrances -> use net penetrances to generate baseline family history (including proband's breast/ovarian cancer status) and exclude probands with breast cancer at baseline (note: baseline family history only contributes to future risk through the carrier probabilities)
# future risk is calculated using modified crude penetrance -> generate baseline BCRAT covariates, use them to modify crude hazard, convert modified crude hazard to modified crude penetrance, sample age to breast cancer from modified crude penetrance (between current age and maximum age)

library(BayesMendel) # 2.1.6.1
library(dplyr) # 0.8.5
source("fam.gen.R")
source("geno.gen.R")
source("pheno.gen.R")
load("cgn.rel.counts.RData")

script_num = 1
prop_aj = 1
# get script_num from command line
args = commandArgs(trailingOnly=T)
if (length(args)>=1) {
  script_num = as.numeric(args[1])
  
  if (length(args)>=2) {
    prop_aj = as.numeric(args[2])/100
  }
}

# number of families to generate
N = 10000

# allow up to 10 relatives of a given type
max.count = 10
cgn.rel.counts[, 2:11] = apply(cgn.rel.counts[, 2:11], 2, function(c) pmin(c, max.count))
cgn.rel.counts[, c("HALF_SIS", "HALF_BRO")] = NULL
cgn.rel.counts = na.omit(cgn.rel.counts)

# allow proband ages between 20 and 85
cgn.ages = cgn.ages[which(cgn.ages>=20 & cgn.ages<=85)]

# sample CGN families
set.seed(script_num)
ind.sample = sample(1:nrow(cgn.rel.counts), N, replace=T)
ind.sample.ages = sample(1:length(cgn.ages), N, replace=T)
probands = data.frame(FamID=1:N, agecur=cgn.ages[ind.sample.ages])
count.name = c("FULL_BRO", "FULL_SIS", "TOT_SONS", "TOT_DAUGHTERS", 
               "TOT_PAT_UNCLES", "TOT_PAT_AUNTS", "TOT_MAT_UNCLES", "TOT_MAT_AUNTS")
probands = cbind(probands, cgn.rel.counts[ind.sample, count.name])

# simulate N families with the maximal structure allowed
# use default BRCAPRO AJ prevalences (v2.1-6) for 7% of families (AJ proportion in CGN)
# use default BRCAPRO non-AJ prevalences for 93% of families
# mean.age.diff=27, sd.age.diff=6 based on distribution of age differences between probands and their mothers in CGN
# mean.age.cens=80, sd.age.cens=15 based on U.S. life expectancy statistics
N.aj = ceiling(prop_aj*N)
set.seed(script_num)
fam.aj = fam.gen(n=pmax(1, N.aj),
              n.unc.pat=max.count, n.aunt.pat=max.count, n.unc.mat=max.count, n.aunt.mat=max.count, 
              n.bro=max.count, n.sis=max.count,
              n.son=max.count, n.daught=max.count,
              n.neph.bro=0, n.nie.bro=0, n.neph.sis=0, n.nie.sis=0,
              ethnic="AJ",
              age.pro=cgn.ages[ind.sample.ages][1:N.aj], 
              mean.age.diff=27, sd.age.diff=6, min.age.diff=14, max.age.diff=40,
              max.age=94,
              mean.age.cens=80, sd.age.cens=15, bl.yr=2001)
fam.aj = geno.gen(fam.aj, pbrca1=0.01366243, pbrca2=0.01168798)
fam.aj = pheno.gen(fam.aj, penetrance = penet.brca.net)



set.seed(script_num)
fam = fam.gen(n=pmax(1, N-N.aj),
                 n.unc.pat=max.count, n.aunt.pat=max.count, n.unc.mat=max.count, n.aunt.mat=max.count, 
                 n.bro=max.count, n.sis=max.count,
                 n.son=max.count, n.daught=max.count,
                 n.neph.bro=0, n.nie.bro=0, n.neph.sis=0, n.nie.sis=0,
                 ethnic="nonAJ",
                 age.pro=cgn.ages[ind.sample.ages][pmin(N.aj+1, N):N], 
                 mean.age.diff=27, sd.age.diff=6, min.age.diff=14, max.age.diff=40,
                 max.age=94,
                 mean.age.cens=80, sd.age.cens=15, bl.yr=2001)
fam = geno.gen(fam, pbrca1=0.0005829, pbrca2=0.000676)
fam = pheno.gen(fam, penetrance = penet.brca.net)

fam$FamID = fam$FamID+N.aj
fam$unique.id = paste0(fam$FamID, "_", fam$ID)
fam$unique.mid = paste0(fam$FamID, "_", fam$MotherID)
fam$unique.mid[which(fam$MotherID==0)] = NA
fam$unique.fid = paste0(fam$FamID, "_", fam$FatherID)
fam$unique.fid[which(fam$FatherID==0)] = NA

if (N.aj==N) {
  fam = fam.aj
} else {
  fam = rbind(fam.aj, fam)
}



# prune the simulated families so that they have the same structure as the sampled CGN families
rel.df = data.frame(firstID=c(60, 70, 80, 90, 10, 20, 30, 40), 
           relationship=c("Brother", "Sister", "Son", "Daughter", 
                          "Paternal Uncle", "Paternal Aunt", "Maternal Uncle", "Maternal Aunt"), 
           count.name=count.name)


fam = left_join(fam, probands[, c("FamID", count.name)], by="FamID")

fam$remove = 0
for (i in 1:length(count.name)) {
  fam$remove[which(fam$relationship %in% rel.df$relationship[i] & fam$ID >= rel.df$firstID[i] + fam[, count.name[i]])] = 1
}


fam = fam[which(fam$remove==0 & fam$agecur>0),]
fam[ , c("remove", count.name, "birth.yr.m", "birth.yr.f")] = NULL

fam$FamID = fam$FamID + (script_num-1)*N
probands$FamID = probands$FamID + (script_num-1)*N

### make proband alive at baseline 
ind.dead = which(fam$relationship=="Self" & fam$Death==1)
fam$agecur[ind.dead] = fam$bl.yr[ind.dead] - fam$birth.yr[ind.dead]
fam$Death[ind.dead] = 0
fam$AgeDeath[ind.dead] = NA
fam$AffectedBreast = as.numeric(fam$AgeBreast.lt <= fam$agecur)
fam$AffectedOvary = as.numeric(fam$AgeOvary.lt <= fam$agecur)
fam$AgeBreast = fam$agecur
fam$AgeOvary = fam$agecur
fam$AgeBreast[which(fam$AffectedBreast==1)] = fam$AgeBreast.lt[which(fam$AffectedBreast==1)]
fam$AgeOvary[which(fam$AffectedOvary==1)] = fam$AgeOvary.lt[which(fam$AffectedOvary==1)]


### sample BCRAT covariates (AgeMen, N_Biop; calculate N_Rels and Age1st from the simulated family history)
load("cgn.gail.vars.RData")
probands = fam[which(fam$relationship=="Self"),]
probands[, c("AgeMen", "N_Biop", "N_Rels", "Age1st", "HypPlas")] = NA
# sample AgeMen and N_Biop
set.seed(script_num)
ind.sample.AgeMen = sample(which(gail.vars$AgeMen>0), N, replace=T)
ind.sample.N_Biop = sample(which(gail.vars$N_Biop>=0), N, replace=T)
probands$AgeMen = gail.vars$AgeMen[ind.sample.AgeMen]
probands$N_Biop = gail.vars$N_Biop[ind.sample.N_Biop]
# calculate N_Rels
aff.rels = data.frame(table(fam$FamID[which(fam$relationship %in% c("Mother", "Sister", "Daughter") & fam$AffectedBreast>0)]))
aff.rels$FamID = as.numeric(as.character(aff.rels$Var1))
probands = left_join(probands, aff.rels[, c("FamID", "Freq")])
probands$N_Rels = probands$Freq
probands$N_Rels[which(probands$N_Rels %in% NA)] = 0
# calculate Age1st
children = fam[which(fam$relationship %in% c("Son", "Daughter")),]
children = children[order(children$FamID, children$birth.yr),]
children = children[which(!duplicated(children$FamID)),]
probands = left_join(probands, children[, c("FamID", "birth.yr")], by="FamID", suffix=c("", ".child"))
probands$Age1st = probands$birth.yr.child - probands$birth.yr
probands$Age1st[which(probands$Age1st %in% NA)] = 98
# check that Age1st >= AgeMen
probands$AgeMen[which(probands$Age1st < probands$AgeMen)] = probands$Age1st[which(probands$Age1st < probands$AgeMen)]
# atypical hyperplasia in 10% of biopsies 
ind.ah = sample(which(probands$N_Biop>0), floor(0.1*sum(probands$N_Biop)), replace=F)
probands$HypPlas[ind.ah] = 1
probands$HypPlas[which(probands$HypPlas %in% NA)] = 99


fam = left_join(fam, probands[, c("unique.id", "AgeMen", "N_Biop", "N_Rels", "Age1st", "HypPlas")])


save(fam, file=paste0("../families/simFam_", prop_aj, "_", script_num, ".RData"))


