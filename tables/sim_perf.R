library(BayesMendel)
library(geepack)
library(xtable)
source('sim_perf_fns.R')


# combine predictions from BRCAPRO+BCRAT (M), BRCAPRO, BCRAT
script_num=1
prop_aj=1
load(paste0("../families/simFam_", prop_aj, "_", script_num, ".RData"))
probands = fam[which(fam$ID==50),]
probands = probands[which(probands$AffectedBreast==0),]
probands = nccn(probands, fam)
load(paste0("../model_results/run_penmod_sim_", prop_aj, "_", script_num, ".RData"))
penmod.results = penmod.results[which(penmod.results$FamID %in% probands$FamID),]
load(paste0("../model_results/run_brcapro_sim_", prop_aj, "_", script_num, ".RData"))
brcapro.results = brcapro.results[which(brcapro.results$FamID %in% probands$FamID),]
load(paste0("../model_results/run_gail_sim_", prop_aj, "_", script_num, ".RData"))
gail.results = gail.results[which(gail.results$FamID %in% probands$FamID),]

probands = left_join(probands, penmod.results)
probands = left_join(probands, brcapro.results[, c("FamID", paste0("brcapro", 1:5))])
probands = left_join(probands, gail.results[, c("FamID", paste0("gail", 1:5))])

probands.all = probands

  
for (script_num in 2:10) {
  load(paste0("../families/simFam_", prop_aj, "_", script_num, ".RData"))
  probands = fam[which(fam$ID==50),]
  probands = probands[which(probands$AffectedBreast==0),]
  probands = nccn(probands, fam)
  load(paste0("../model_results/run_penmod_sim_", prop_aj, "_", script_num, ".RData"))
  penmod.results = penmod.results[which(penmod.results$FamID %in% probands$FamID),]
  load(paste0("../model_results/run_brcapro_sim_", prop_aj, "_", script_num, ".RData"))
  brcapro.results = brcapro.results[which(brcapro.results$FamID %in% probands$FamID),]
  load(paste0("../model_results/run_gail_sim_", prop_aj, "_", script_num, ".RData"))
  gail.results = gail.results[which(gail.results$FamID %in% probands$FamID),]
  
  probands = left_join(probands, penmod.results)
  probands = left_join(probands, brcapro.results[, c("FamID", paste0("brcapro", 1:5))])
  probands = left_join(probands, gail.results[, c("FamID", paste0("gail", 1:5))])
  
  probands.all = rbind(probands.all, probands)
}

save(probands.all, file="combine_results.RData")



# generate outcomes
set.seed(12345)
probands2 = generate.outcomes(probands.all, times=1:5)

# split into training and test sets
ind.test = setdiff(1:nrow(probands2), 1:50000)
pro.test = probands2[ind.test, ]
pro.train = probands2[setdiff(1:nrow(probands2), ind.test),]

table(pro.test$BC.5)
table(pro.train$BC.5)

save(pro.test, pro.train, file="pro.sim.RData")

# fit logistic regression ensemble
lr.fit = glm(BC.5 ~ I(sqrt(gail5))*I(sqrt(brcapro5)), data=probands2[setdiff(1:nrow(probands2), ind.test),], family=binomial)
pro.test$lr5 = predict(lr.fit, newdata=pro.test, type="response")


# fit longitudinal logistic regression ensemble
pro.train$id = pro.train$FamID
pro.train.long = convert.to.long(pro.train)
fit.geeglm = geeglm(BC ~ I(sqrt(brcapro))*I(sqrt(gail))*time, id=id, family=binomial(link="logit"), corstr="ar1", data=pro.train.long)
# get predictions
pro.test$id = pro.test$FamID
pro.test$long5 = predict.gee(fit.geeglm, convert.to.long(pro.test, times=5))



# takes about 10 minutes to do 1000 bootstrap replicates 
perf.sim = perf.boot.2(pro.test[, c("penmod5", "lr5", "long5", "brcapro5", "gail5")], pro.test$BC.5, cens.dist="none", dataset=pro.test, covars=NULL, t=5, model.names=c("B+B (M)", "B+B (E)", "B+B (E2)", "BRCAPRO", "BCRAT"), nboot=1000, ref.col=4)

# pairwise comparisons across bootstrap replicates
oe = data.frame(perf.sim[[2]][2:1001,])
auc = data.frame(perf.sim[[3]][2:1001,])
bs = data.frame(perf.sim[[4]][2:1001,])
snb = data.frame(perf.sim[[5]][2:1001,])
lsr = data.frame(perf.sim[[6]][2:1001,])

boot.table = data.frame(comparison=c("B+B(M)>B+B(E)", "B+B(M)>B+B(E2)", "B+B(M)>BRCAPRO", "B+B(M)>BCRAT", "B+B(E)>B+B(E2)", "B+B(E)>BRCAPRO", "B+B(E)>BCRAT", "B+B(E2)>BRCAPRO", "B+B(E2)>BCRAT"), OE=NA, AUC=NA, SNB=NA, BS=NA, LSR=NA)
boot.table[, "AUC"] = c(sapply(2:5, function(x) length(which(auc[,1] > auc[,x]))), sapply(3:5, function(x) length(which(auc[,2] > auc[,x]))),  sapply(4:5, function(x) length(which(auc[,3] > auc[,x]))))
boot.table[, "BS"] = c(sapply(2:5, function(x) length(which(bs[,1] < bs[,x]))), sapply(3:5, function(x) length(which(bs[,2] < bs[,x]))), sapply(4:5, function(x) length(which(bs[,3] < bs[,x]))))
boot.table[, "OE"] = c(sapply(2:5, function(x) length(which( abs(oe[,1]-1) < abs(oe[,x]-1) ))), sapply(3:5, function(x) length(which(abs(oe[,2]-1) < abs(oe[,x]-1) ))), sapply(4:5, function(x) length(which(abs(oe[,3]-1) < abs(oe[,x]-1) ))))
boot.table[, "SNB"] = c(sapply(2:5, function(x) length(which(snb[,1] > snb[,x]))), sapply(3:5, function(x) length(which(snb[,2] > snb[,x]))),  sapply(4:5, function(x) length(which(snb[,3] > snb[,x]))))
boot.table[, "LSR"] = c(sapply(2:5, function(x) length(which(lsr[,1] > lsr[,x]))), sapply(3:5, function(x) length(which(lsr[,2] > lsr[,x]))),  sapply(4:5, function(x) length(which(lsr[,3] > lsr[,x]))))

boot.table[, 2:6] = boot.table[, 2:6]/1000

save(perf.sim, boot.table, file="sim_results.RData")

### generate latex table

# first section: performance measures
table.perf = format.perf.table.2(perf.sim[[1]], model.names = c("B+B (M)", "B+B (E)", "B+B (E2)", "BRCAPRO", "BCRAT"))

# second section: comparisons across bootstrap replicates
# apply shading to cells
table.boot = boot.table
table.boot[, 2:6] = formatC(as.matrix(table.boot[, 2:6]), 3, format="f")
table.boot[, 2:6] = apply(table.boot[, 2:6], 2, function(x) paste0("\\mycolor{", x, "}"))

# add sub-table headers
addtorow = list()
addtorow$pos = as.list(c(0, 6))
addtorow$command = as.vector(c('\\hline \\textbf{Performance Metrics} & \\\\ \n', paste0(paste0('\\hline \\multicolumn{3}{l}{\\textbf{Comparisons Across Bootstrap Replicates}} ', collapse=''), '& & & \\\\ \n')))

# combine sub-tables
rownames(table.boot) = table.boot$comparison
table.boot$comparison = NULL
colnames(table.boot) = colnames(table.perf)
table.overall = rbind(table.perf, table.boot)
colnames(table.overall)[4:5] = c("$\\Delta$BS", "$\\Delta$LS")


print(xtable(table.overall, align='llllll', caption='5-year performance in a simulated dataset with 45,557 probands (717 cases). B+B: BRCAPRO+BCRAT. $\\Delta$BS: \\% relative improvement in Brier Score compared to BRCAPRO. $\\Delta$LS: \\% relative improvement in logarithmic score compared to BRCAPRO. The ``Comparisons Across Bootstrap Replicates" section shows pairwise comparisons involving the combination models across 1000 bootstrap replicates of the validation dataset; the row for $A>B$ shows the proportion of bootstrap replicates where model A outperformed model B with respect to each metric. Proportions $>0.5$ are highlighted in blue (with darker shades of blue for higher proportions) and proportions $\\leq 0.5$ are highlighted in red (with darker shades of red for lower proportions).', label='table:table_sim'), add.to.row=addtorow, sanitize.text.function = identity, hline.after=c(-1, nrow(rbind(table.perf, table.boot))), size="\\footnotesize", table.placement="!htbp", file="table_sim.tex")



### calibration plots

png('oe_plots_sim.png', width=1000, height=1000, res=150)
oe.plot = plot.oe(pro.test[, c("penmod5", "lr5", "long5", "brcapro5", "gail5")], pro.test$BC.5, n.bins=10, model.names=c("BRCAPRO+BCRAT (M)", "BRCAPRO+BCRAT (E)", "BRCAPRO+BCRAT (E2)", "BRCAPRO", "BCRAT"),
                  xlimits=c(0, 0.06), ylimits=c(0, 0.06), bar.width=0.0025, threshold=0.00001,
                  cols=2)
dev.off()


