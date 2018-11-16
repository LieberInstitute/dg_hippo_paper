##### libraries
library(SummarizedExperiment)
library(jaffelab)
library(MatrixEQTL)
library(sva)
library(limma)
library(edgeR)

######################
## load data #########
######################

## filtered counts
load("count_data/merged_dg_hippo_allSamples_n596.rda")
colnames(rse_deg_joint) = colnames(rse_gene_joint)

##### get MDS #####
mds_dg = read.table("/dcl01/ajaffe/data/lab/dg_hippo/genotype_data/Astellas_DG_Genotypes_n263_maf05_geno10_hwe1e6.mds",
	header=TRUE,as.is=TRUE,row.names=1)[,-(1:2)]
mds_dg = mds_dg[rse_gene_joint$BrNum[rse_gene_joint$Dataset == "DG"],1:5]

load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/genotype_data/mds_extracted_from_BrainSeq_Phase2_RiboZero_Genotypes_n551.Rdata",ver=TRUE)
mds_hippo = mds[rse_gene_joint$BrNum[rse_gene_joint$Dataset == "Hippo"],1:5]
colnames(mds_dg) = colnames(mds_hippo)
mds = rbind(mds_dg, mds_hippo)

# relevel
rse_gene_joint$Dx = factor(rse_gene_joint$Dx, levels = c("Control", "Schizo", "Bipolar", "MDD"))

#####################
## get qSVs #########
#####################

mod = model.matrix(~Dx + Region + totalAssignedGene + Sex + Age + 
				mitoRate + as.matrix(mds[,1:5]),
				data = colData(rse_gene_joint))
# do qSVA
covCount = assays(rse_deg_joint[rowData(rse_deg_joint)$bonfSig,])$counts
k = num.sv(log2(covCount+1), mod)
qSVs = prcomp(t(log2(covCount+1)))$x[,1:k]
colnames(qSVs) = paste0("qSV", 1:k)

###############################
######### analysis ############
###############################

###### GENE LEVEL ###########

#### dg only #####
rse_gene_dg = rse_gene_joint[,rse_gene_joint$Dataset == "DG"]
qSVs_dg = as.matrix(qSVs[colnames(rse_gene_dg),])
mds_dg = as.matrix(mds[rse_gene_dg$BrNum,])

mod_dg = model.matrix(~ Age + Dx + totalAssignedGene + Sex + 
				mitoRate + mds_dg + qSVs_dg,
				data = colData(rse_gene_dg))

dge_dg <- DGEList(counts = assays(rse_gene_dg)$counts)
dge_dg <- calcNormFactors(dge_dg)
vGene_dg <- voom(dge_dg, mod_dg, plot = TRUE)
geneFit_dg = lmFit(vGene_dg, mod_dg)
geneFit_dg = eBayes(geneFit_dg)
geneAgeStats_dg = topTable(geneFit_dg, coef = 2, n = nrow(dge_dg), sort="none")
geneSzStats_dg = topTable(geneFit_dg, coef = 3, n = nrow(dge_dg), sort="none")
## other dx
geneBpdStats_dg = topTable(geneFit_dg, coef = 4, n = nrow(dge_dg), sort="none")
geneMddStats_dg = topTable(geneFit_dg, coef = 5, n = nrow(dge_dg), sort="none")

### just dg for all 3 disorders
vars = c(1,3,4,5)

tmp = cbind(geneSzStats_dg[,vars], geneBpdStats_dg[,vars], geneMddStats_dg[,vars])
colnames(tmp) = paste0(rep(c("SZ","BP", "MD"), each=length(vars)),"_", colnames(tmp))
geneDxStats_dg = rowRanges(rse_gene_joint)
mcols(geneDxStats_dg) = cbind(tmp, mcols(geneDxStats_dg))
save(geneDxStats_dg, file = "rdas/geneLevel_dxEffects_dg.rda")

#### hippo only #####
rse_gene_hippo = rse_gene_joint[,rse_gene_joint$Dataset == "Hippo"]
qSVs_hippo = as.matrix(qSVs[colnames(rse_gene_hippo),])
mds_hippo = as.matrix(mds[rse_gene_hippo$BrNum,])

mod_hippo = model.matrix(~ Age + Dx + totalAssignedGene + Sex + 
				mitoRate + mds_hippo + qSVs_hippo,
				data = colData(rse_gene_hippo))

dge_hippo <- DGEList(counts = assays(rse_gene_hippo)$counts)
dge_hippo <- calcNormFactors(dge_hippo)
vGene_hippo <- voom(dge_hippo, mod_hippo, plot = TRUE)
geneFit_hippo = lmFit(vGene_hippo, mod_hippo)
geneFit_hippo = eBayes(geneFit_hippo)
geneAgeStats_hippo = topTable(geneFit_hippo, coef = 2, n = nrow(dge_hippo), sort="none")
geneSzStats_hippo = topTable(geneFit_hippo, coef = 3, n = nrow(dge_hippo), sort="none")

plot(geneAgeStats_dg$t, geneAgeStats_hippo$t)
plot(geneSzStats_dg$t, geneSzStats_hippo$t)

#### joint interaction #####
mod_joint_age = model.matrix(~ Age*Region + Dx + totalAssignedGene + Sex + 
				mitoRate + as.matrix(mds) + qSVs,	data = colData(rse_gene_joint))

dge_joint <- DGEList(counts = assays(rse_gene_joint)$counts)
dge_joint <- calcNormFactors(dge_joint)
vGene_joint <- voom(dge_joint, mod_joint_age, plot = TRUE)

## duplicate correlation
corfit_gene <- duplicateCorrelation(vGene_joint$E, mod_joint_age, block=rse_gene_joint$BrNum)
dir.create("duplCorr_objects")
save(corfit_gene, file = "duplCorr_objects/geneLevel_dgPlusHippo_age17.rda")

## dup corr
geneFit_joint = lmFit(vGene_joint, mod_joint_age, block=rse_gene_joint$BrNum,
        correlation = corfit_gene$consensus.correlation)

geneFit_joint = eBayes(geneFit_joint)
geneAgeStats_int = topTable(geneFit_joint, coef = ncol(mod_joint_age), 
	n = nrow(dge_joint), sort="none")

## dup corr
mod_joint_dx = model.matrix(~ Age + Dx*Region + totalAssignedGene + Sex + 
				mitoRate + as.matrix(mds) + qSVs,	data = colData(rse_gene_joint))

geneFit_joint_dx = lmFit(vGene_joint, mod_joint_dx, block=rse_gene_joint$BrNum,
        correlation = corfit_gene$consensus.correlation)

geneFit_joint_dx = eBayes(geneFit_joint_dx)
geneSzStats_int = topTable(geneFit_joint_dx, coef = 24, 
	n = nrow(dge_joint), sort="none")

### merge
tmp = cbind(geneSzStats_dg[,vars], geneSzStats_hippo[,vars], geneSzStats_int[,vars])
colnames(tmp) = paste0(colnames(tmp), "_SZ_", rep(c("DG", "Hippo", "Inter"), each=4))
geneSzStats = rowRanges(rse_gene_joint)
mcols(geneSzStats) = cbind(tmp, mcols(geneSzStats))

tmp = cbind(geneAgeStats_dg[,vars], geneAgeStats_hippo[,vars], geneAgeStats_int[,vars])
colnames(tmp) = paste0(colnames(tmp), "_Age_", rep(c("DG", "Hippo", "Inter"), each=4))
geneAgeStats = rowRanges(rse_gene_joint)
mcols(geneAgeStats) = cbind(tmp, mcols(geneAgeStats))

save(geneSzStats,geneAgeStats,
	file = "rdas/geneLevel_ageAndSzInteraction.rda")

###### EXON LEVEL ###########

#### dg only #####
rse_exon_dg = rse_exon_joint[,rse_exon_joint$Dataset == "DG"]
dee_dg <- DGEList(counts = assays(rse_exon_dg)$counts)
dee_dg <- calcNormFactors(dee_dg)
vExon_dg <- voom(dee_dg, mod_dg, plot = TRUE)
exonFit_dg = lmFit(vExon_dg, mod_dg)
exonFit_dg = eBayes(exonFit_dg)
exonAgeStats_dg = topTable(exonFit_dg, coef = 2, n = nrow(dee_dg), sort="none")
exonSzStats_dg = topTable(exonFit_dg, coef = 3, n = nrow(dee_dg), sort="none")

#### hippo only #####
rse_exon_hippo = rse_exon_joint[,rse_exon_joint$Dataset == "Hippo"]
dee_hippo <- DGEList(counts = assays(rse_exon_hippo)$counts)
dee_hippo <- calcNormFactors(dee_hippo)
vExon_hippo <- voom(dee_hippo, mod_hippo, plot = TRUE)
exonFit_hippo = lmFit(vExon_hippo, mod_hippo)
exonFit_hippo = eBayes(exonFit_hippo)
exonAgeStats_hippo = topTable(exonFit_hippo, coef = 2, n = nrow(dee_hippo), sort="none")
exonSzStats_hippo = topTable(exonFit_hippo, coef = 3, n = nrow(dee_hippo), sort="none")

plot(exonAgeStats_dg$t, exonAgeStats_hippo$t)
plot(exonSzStats_dg$t, exonSzStats_hippo$t)

#### joint interaction #####
dee_joint <- DGEList(counts = assays(rse_exon_joint)$counts)
dee_joint <- calcNormFactors(dee_joint)
vExon_joint <- voom(dee_joint, mod_joint_age, plot = TRUE)

## duplicate correlation
corfit_exon <- duplicateCorrelation(vExon_joint$E, mod_joint_age, block=rse_exon_joint$BrNum)
save(corfit_exon, file = "duplCorr_objects/exonLevel_dgPlusHippo_age17.rda")

## dup corr
exonFit_joint = lmFit(vExon_joint, mod_joint_age, block=rse_exon_joint$BrNum,
        correlation = corfit_exon$consensus.correlation)

exonFit_joint = eBayes(exonFit_joint)
exonAgeStats_int = topTable(exonFit_joint, coef = ncol(mod_joint_age), 
	n = nrow(dee_joint), sort="none")

## dup corr
exonFit_joint_dx = lmFit(vExon_joint, mod_joint_dx, block=rse_exon_joint$BrNum,
        correlation = corfit_exon$consensus.correlation)

exonFit_joint_dx = eBayes(exonFit_joint_dx)
exonSzStats_int = topTable(exonFit_joint_dx, coef = 24, 
	n = nrow(dee_joint), sort="none")

### merge
vars = c(1,3,4,5)
tmp = cbind(exonSzStats_dg[,vars], exonSzStats_hippo[,vars], exonSzStats_int[,vars])
colnames(tmp) = paste0(colnames(tmp), "_SZ_", rep(c("DG", "Hippo", "Inter"), each=4))
exonSzStats = rowRanges(rse_exon_joint)
mcols(exonSzStats) = cbind(tmp, mcols(exonSzStats))

tmp = cbind(exonAgeStats_dg[,vars], exonAgeStats_hippo[,vars], exonAgeStats_int[,vars])
colnames(tmp) = paste0(colnames(tmp), "_Age_", rep(c("DG", "Hippo", "Inter"), each=4))
exonAgeStats = rowRanges(rse_exon_joint)
mcols(exonAgeStats) = cbind(tmp, mcols(exonAgeStats))

save(exonSzStats,exonAgeStats,
	file = "rdas/exonLevel_ageAndSzInteraction.rda")

#################################
###### Junction LEVEL ###########
#################################

#### dg only #####
rse_jx_dg = rse_jxn_joint[,rse_jxn_joint$Dataset == "DG"]
dje_dg <- DGEList(counts = assays(rse_jx_dg)$counts)
dje_dg <- calcNormFactors(dje_dg)
vJxn_dg <- voom(dje_dg, mod_dg, plot = TRUE)
jxFit_dg = lmFit(vJxn_dg, mod_dg)
jxFit_dg = eBayes(jxFit_dg)
jxAgeStats_dg = topTable(jxFit_dg, coef = 2, n = nrow(dje_dg), sort="none")
jxSzStats_dg = topTable(jxFit_dg, coef = 3, n = nrow(dje_dg), sort="none")

#### hippo only #####
rse_jx_hippo = rse_jxn_joint[,rse_jxn_joint$Dataset == "Hippo"]
dje_hippo <- DGEList(counts = assays(rse_jx_hippo)$counts)
dje_hippo <- calcNormFactors(dje_hippo)
vJxn_hippo <- voom(dje_hippo, mod_hippo, plot = TRUE)
jxFit_hippo = lmFit(vJxn_hippo, mod_hippo)
jxFit_hippo = eBayes(jxFit_hippo)
jxAgeStats_hippo = topTable(jxFit_hippo, coef = 2, n = nrow(dje_hippo), sort="none")
jxSzStats_hippo = topTable(jxFit_hippo, coef = 3, n = nrow(dje_hippo), sort="none")

plot(jxAgeStats_dg$t, jxAgeStats_hippo$t)
plot(jxSzStats_dg$t, jxSzStats_hippo$t)

#### joint interaction #####
dje_joint <- DGEList(counts = assays(rse_jxn_joint)$counts)
dje_joint <- calcNormFactors(dje_joint)
vJxn_joint <- voom(dje_joint, mod_joint_age, plot = TRUE)

## duplicate correlation
corfit_jx <- duplicateCorrelation(vJxn_joint$E, mod_joint_age, block=rse_jxn_joint$BrNum)
save(corfit_jx, file = "duplCorr_objects/jxLevel_dgPlusHippo_age17.rda")

## dup corr
jxFit_joint = lmFit(vJxn_joint, mod_joint_age, block=rse_jxn_joint$BrNum,
        correlation = corfit_jx$consensus.correlation)

jxFit_joint = eBayes(jxFit_joint)
jxAgeStats_int = topTable(jxFit_joint, coef = ncol(mod_joint_age), 
	n = nrow(dje_joint), sort="none")

## dup corr
jxFit_joint_dx = lmFit(vJxn_joint, mod_joint_dx, block=rse_jxn_joint$BrNum,
        correlation = corfit_jx$consensus.correlation)

jxFit_joint_dx = eBayes(jxFit_joint_dx)
jxSzStats_int = topTable(jxFit_joint_dx, coef = 24, 
	n = nrow(dje_joint), sort="none")

### merge
vars = c(1,3,4,5)
tmp = cbind(jxSzStats_dg[,vars], jxSzStats_hippo[,vars], jxSzStats_int[,vars])
colnames(tmp) = paste0(colnames(tmp), "_SZ_", rep(c("DG", "Hippo", "Inter"), each=4))
jxSzStats = rowRanges(rse_jxn_joint)
mcols(jxSzStats) = cbind(tmp, mcols(jxSzStats))

tmp = cbind(jxAgeStats_dg[,vars], jxAgeStats_hippo[,vars], jxAgeStats_int[,vars])
colnames(tmp) = paste0(colnames(tmp), "_Age_", rep(c("DG", "Hippo", "Inter"), each=4))
jxAgeStats = rowRanges(rse_jxn_joint)
mcols(jxAgeStats) = cbind(tmp, mcols(jxAgeStats))

save(jxSzStats,jxAgeStats,
	file = "rdas/jxLevel_ageAndSzInteraction.rda")

#################################
###### Tx LEVEL ###########
#################################

#### dg only #####
rse_tx_dg = rse_tx_joint[,rse_tx_joint$Dataset == "DG"]
txExprsDg = log2(assays(rse_tx_dg)$tpm+1)
txFit_dg = lmFit(txExprsDg, mod_dg)
txFit_dg = eBayes(txFit_dg)
txAgeStats_dg = topTable(txFit_dg, coef = 2, n = nrow(txExprsDg), sort="none")
txSzStats_dg = topTable(txFit_dg, coef = 3, n = nrow(txExprsDg), sort="none")

#### hippo only #####
rse_tx_hippo = rse_tx_joint[,rse_tx_joint$Dataset == "Hippo"]
txExprsHippo = log2(assays(rse_tx_hippo)$tpm+1)
txFit_hippo = lmFit(txExprsHippo, mod_hippo)
txFit_hippo = eBayes(txFit_hippo)
txAgeStats_hippo = topTable(txFit_hippo, coef = 2, n = nrow(txExprsHippo), sort="none")
txSzStats_hippo = topTable(txFit_hippo, coef = 3, n = nrow(txExprsHippo), sort="none")

plot(txAgeStats_dg$t, txAgeStats_hippo$t)
plot(txSzStats_dg$t, txSzStats_hippo$t)

#### joint interaction #####
txExprsJoint= log2(assays(rse_tx_joint)$tpm+1)

## duplicate correlation
corfit_tx <- duplicateCorrelation(txExprsJoint, mod_joint_age, block=rse_tx_joint$BrNum)
save(corfit_tx, file = "duplCorr_objects/txLevel_dgPlusHippo_age17.rda")

## dup corr
txFit_joint = lmFit(vTx_joint, mod_joint_age, block=rse_tx_joint$BrNum,
        correlation = corfit_tx$consensus.correlation)

txFit_joint = eBayes(txFit_joint)
txAgeStats_int = topTable(txFit_joint, coef = ncol(mod_joint_age), 
	n = nrow(vTx_joint), sort="none")

## dup corr
txFit_joint_dx = lmFit(vTx_joint, mod_joint_dx, block=rse_tx_joint$BrNum,
        correlation = corfit_tx$consensus.correlation)

txFit_joint_dx = eBayes(txFit_joint_dx)
txSzStats_int = topTable(txFit_joint_dx, coef = 24, 
	n = nrow(vTx_joint), sort="none")

### merge
vars = c(1,3,4,5)
tmp = cbind(txSzStats_dg[,vars], txSzStats_hippo[,vars], txSzStats_int[,vars])
colnames(tmp) = paste0(colnames(tmp), "_SZ_", rep(c("DG", "Hippo", "Inter"), each=4))
txSzStats = rowRanges(rse_tx_joint)
mcols(txSzStats) = cbind(tmp, mcols(txSzStats))

tmp = cbind(txAgeStats_dg[,vars], txAgeStats_hippo[,vars], txAgeStats_int[,vars])
colnames(tmp) = paste0(colnames(tmp), "_Age_", rep(c("DG", "Hippo", "Inter"), each=4))
txAgeStats = rowRanges(rse_tx_joint)
mcols(txAgeStats) = cbind(tmp, mcols(txAgeStats))

save(txSzStats,txAgeStats,
	file = "rdas/txLevel_ageAndSzInteraction.rda")
