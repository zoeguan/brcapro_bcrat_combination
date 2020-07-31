library(BCRA) # 2.1

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


# function to get t-year Gail risk 
gail.risk.tyear = function(input, t) {
  input$T2 = input$T1 + t
  return(absolute.risk(input))
}



probands = fam[which(fam$ID==50),]
probands$T1 = probands$agecur
probands$T2 = probands$T1 + 5
probands$Race = 4
probands$ID = probands$FamID

for (t in 1:5) {
  probands[, paste0("gail", t)] = gail.risk.tyear(probands, t)
}
probands[, paste0("gail", 1:5)] = probands[, paste0("gail", 1:5)]/100

gail.results = probands[, c("FamID", "T1", "T2", "AgeMen", "Age1st", "N_Biop", "HypPlas", "N_Rels", "Race", paste0("gail", 1:5))]


save(gail.results, file=paste0("../model_results/run_gail_sim_", prop_aj, "_", script_num, ".RData"))


