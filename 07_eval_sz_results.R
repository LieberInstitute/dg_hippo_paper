###
library(GenomicRanges)
library(clusterProfiler)
library(limma)
library(RColorBrewer)

## load results
load("rdas/geneLevel_dxEffects_dg.rda")
load("rdas/geneLevel_ageAndSzInteraction.rda")
rm(geneAgeStats)

## load DEqual plot stats
load("rdas/Hippocampus_geneLevel_degradationStats_forDEqual_hg38.rda")
degradeStats = degradeStats[names(geneSzStats),]

#######################
## focus on dx in dg ##
#######################

sum(geneDxStats_dg$SZ_adj.P.Val < 0.1)
sum(geneDxStats_dg$BP_adj.P.Val < 0.1)
sum(geneDxStats_dg$MD_adj.P.Val < 0.1)

geneDxStats_dg$SZ_sig = geneDxStats_dg$SZ_adj.P.Val < 0.1
geneDxStats_dg$BP_sig = geneDxStats_dg$BP_adj.P.Val < 0.1
geneDxStats_dg$MD_sig = geneDxStats_dg$MD_adj.P.Val < 0.1
vc = vennCounts(as.data.frame(geneDxStats_dg)[,27:29])
vennDiagram(vc)  ## little overlap

suppTab = geneDxStats_dg
suppTab = as.data.frame(suppTab)
suppTab$gencodeTx = NULL
colnames(suppTab)[1:3] = c("chr_hg38", "start_hg38", "end_hg38")
suppTab = suppTab[,c(19,22, 26:28, 6:18,20:21,1:5, 23:25)]
write.csv(suppTab, file = "tables/caseControl_geneLevel_stats_noFilter.csv",row.names=FALSE)

## significane
suppTabSig = suppTab[suppTab$SZ_sig |suppTab$BP_sig | suppTab$MD_sig ,]
write.csv(suppTabSig, file = "tables/caseControl_geneLevel_stats_anySig.csv",row.names=FALSE)

## which shared
geneDxStats_dg[geneDxStats_dg$SZ_sig & geneDxStats_dg$BP_sig ,]

fcs = mcols(geneDxStats_dg)[,grep("logFC", colnames(mcols(geneDxStats_dg)))]
cor(as.matrix(fcs))

pdf("plots/dxPairs.pdf", h=10,w=10)
pairs(fcs,pch=21,bg="grey", cex.axis=2,cex.main=1.4,
	ylim=c(-1.5,1.5),xlim=c(-1.5,1.5),
	main = "Cross-Disorders Overlap: Expressed Genes")
dev.off()

## plots
## what genes?
geneDxStats_dg$Symbol[geneDxStats_dg$SZ_sig]
geneDxStats_dg$Symbol[geneDxStats_dg$BP_sig]
geneDxStats_dg$Symbol[geneDxStats_dg$MD_sig]

## compare to hippo
head(geneSzStats)
table(geneSzStats$adj.P.Val_SZ_DG < 0.1, 
	geneSzStats$adj.P.Val_SZ_Hippo < 0.1,
	dnn = c("DG", "HIPPO"))

# which geneSzStats
geneSzStats[geneSzStats$adj.P.Val_SZ_DG < 0.1 & geneSzStats$adj.P.Val_SZ_Hippo < 0.1 ,]
cor(geneSzStats$t_SZ_DG, geneSzStats$t_SZ_Hippo) 
cor(geneSzStats$logFC_SZ_DG, geneSzStats$logFC_SZ_Hippo) 

## colors
cols = rep(1, length(geneSzStats))
cols[geneSzStats$adj.P.Val_SZ_DG < 0.1 & geneSzStats$adj.P.Val_SZ_Hippo > 0.1] = 3
cols[geneSzStats$adj.P.Val_SZ_DG > 0.1 & geneSzStats$adj.P.Val_SZ_Hippo < 0.1] = 2
cols[geneSzStats$adj.P.Val_SZ_DG < 0.1 & geneSzStats$adj.P.Val_SZ_Hippo < 0.1] = 4

names(cols)[cols==1] = "None"
names(cols)[cols==3] = "DG"
names(cols)[cols==2] = "HIPPO"
names(cols)[cols==4] = "Both"

pdf("plots/szAnalysis_dgVsHippo_geneLevel.pdf")
par(mar=c(5,6,3,2),cex.axis=2,cex.lab=2)
palette(c("grey", brewer.pal(3, "Set1")))
plot(geneSzStats$t_SZ_DG, geneSzStats$t_SZ_Hippo,
	ylim = c(-6, 6), xlim = c(-6, 6), pch = 21,bg=cols,
	xlab = "DG (SCZD T-stat)", ylab = "HIPPO (SCZD T-stat)")
plot(geneSzStats$logFC_SZ_DG, geneSzStats$logFC_SZ_Hippo,
	ylim = c(-1, 1), xlim = c(-1, 1), pch = 21,
	bg=cols, xlab = "DG (log2FC)", ylab = "HIPPO (log2FC)")
dev.off()

table(abs(geneSzStats$logFC_SZ_DG) > abs(geneSzStats$logFC_SZ_Hippo))

## degradation check
plot(geneSzStats$t_SZ_DG, degradeStats$t)
cor(geneSzStats$logFC_SZ_DG, degradeStats$logFC)
plot(geneSzStats$logFC_SZ_DG, degradeStats$logFC)
plot(geneSzStats$t_SZ_Hippo, degradeStats$t)

## compare to brainseq phase2 results
load("rdas/dxStats_hippo_filtered_qSVA_noHGoldQSV_matchHIPPO.rda",verbose=TRUE)
rm(outExon, outJxn, outTx)
outGene = outGene[names(geneSzStats),]
plot(geneSzStats$t_SZ_Hippo, outGene$t)
abline(0,1,lty=2,col="blue")
plot(geneSzStats$logFC_SZ_Hippo, outGene$logFC)
abline(0,1,lty=2,col="blue")

table(cols, outGene$adj.P.Val < 0.05)


## interaction
table(geneSzStats$adj.P.Val_SZ_Inter < 0.1,names(cols))

# intIndex = which(geneSzStats$adj.P.Val_Sz_Inter < 0.05)
intIndex = which(geneSzStats$adj.P.Val_SZ_Inter < 0.1)
table(names(cols)[intIndex])

pdf("plots/szAnalysis_dgVsHippo_geneLevel_interaction.pdf")
par(mar=c(5,6,3,2),cex.axis=2,cex.lab=2)
palette(c("grey", brewer.pal(3, "Set1")))
plot(geneSzStats$t_SZ_DG[intIndex], geneSzStats$t_SZ_Hippo[intIndex],
	ylim = c(-8, 8), xlim = c(-8, 8), pch = 21,bg=cols[intIndex],
	xlab = "DG", ylab = "HIPPO", main = "T-statistics")
plot(geneSzStats$logFC_SZ_DG[intIndex], geneSzStats$logFC_SZ_Hippo[intIndex],
	ylim = c(-1, 1), xlim = c(-1, 1), pch = 21,
	bg=cols[intIndex], xlab = "DG", ylab = "HIPPO", main = "log2 Fold Changes")
dev.off()

#####################
### gene ontology ###
geneUniverse = geneSzStats$EntrezID
geneUniverse = as.character(geneUniverse[!is.na(geneUniverse)])

## make sets
geneListsDg = split(geneSzStats$EntrezID[cols %in% c(2,4)], 
	sign(geneSzStats$t_SZ_DG[cols %in% c(2,4)]))
names(geneListsDg) = c("DG_Down", "DG_Up")
geneListsHippo = split(geneSzStats$EntrezID[cols %in% c(3,4)], 
	sign(geneSzStats$t_SZ_Hippo[cols %in% c(3,4)]))
names(geneListsHippo) = c("Hippo_Down", "Hippo_Up")
geneList = c(geneListsDg, geneListsHippo)
geneList$Int = geneSzStats$EntrezID[geneSzStats$adj.P.Val_SZ_Inter < 0.05]
geneList = lapply(geneList, function(x) as.character(x[!is.na(x)]))

## GO
goEnr <- compareCluster(geneList, fun = "enrichGO",
                universe = geneUniverse, OrgDb = 'org.Hs.eg.db',
                ont = "ALL", pAdjustMethod = "none",
                pvalueCutoff  = 1, qvalueCutoff  = 1,
				readable= TRUE)
goDf = as.data.frame(goEnr)

## KEGG
keggEnr <- compareCluster(geneList, fun = "enrichKEGG",
                universe = geneUniverse, pAdjustMethod = "none",
                pvalueCutoff  = 1, qvalueCutoff  = 1)
keggDf = as.data.frame(keggEnr)
keggDf$ONTOLOGY="KEGG"

## merge
gsDf = rbind(goDf, keggDf[,colnames(keggDf)])

uniqueIDs = unique(gsDf$ID)
gsDfList = split(gsDf, gsDf$Cluster)
gsDat = sapply(gsDfList, function(x) x[match(uniqueIDs,x$ID), "pvalue"])
rownames(gsDat) = uniqueIDs
gsDat = as.data.frame(gsDat)
gsDat$Description = gsDf$Description[match(uniqueIDs, gsDf$ID)]
save(gsDat, keggEnr, goEnr, file = "rdas/sz_gene_set_results.rda")


##########################
## plots some examples ###
##########################
library(SummarizedExperiment)
library(recount)
library(sva)
library(jaffelab)
library(edgeR)
library(limma)

## load counts
load("count_data/merged_dg_hippo_allSamples_n596.rda", verbose=TRUE)
rm(rse_exon_joint, rse_tx_joint, rse_jxn_joint)
colnames(rse_deg_joint) = colnames(rse_gene_joint)

##### get MDS #####
mds_dg = read.table("genotype_data/Astellas_DG_Genotypes_n263_maf05_geno10_hwe1e6.mds",
	header=TRUE,as.is=TRUE,row.names=1)[,-(1:2)]
mds_dg = mds_dg[rse_gene_joint$BrNum[rse_gene_joint$Dataset == "DG"],1:5]

load("genotype_data/mds_extracted_from_BrainSeq_Phase2_RiboZero_Genotypes_n551.Rdata",ver=TRUE)
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
k = num.sv(log2(covCount+1), mod) # 
qSVs = prcomp(t(log2(covCount+1)))$x[,1:k]
colnames(qSVs) = paste0("qSV", 1:k)

# library(LIBDpheno) # internal db
# pheno = toxicant[[1]]
# pheno$BrNum = paste0("Br", pheno$brnumerical) 
# pheno = pheno[pheno$BrNum%in% rse_gene_joint$BrNum,]
# pheno = pheno[,c("BrNum", "antidepressants_ssris", "antipsychotics")]
# write.csv(pheno, file = "tables/subject_drug_hx.csv", row.names=FALSE)
drug = read.csv("tables/subject_drug_hx.csv",as.is=TRUE)
rse_gene_joint$SSRI = drug$antidepressants_ssris[match(rse_gene_joint$BrNum, drug$BrNum)]
rse_gene_joint$antipsychotics = drug$antipsychotics[match(rse_gene_joint$BrNum, drug$BrNum)]

## split out dg
rse_gene_dg = rse_gene_joint[,rse_gene_joint$Dataset == "DG"]
qSVs_dg = as.matrix(qSVs[colnames(rse_gene_dg),])
mds_dg = as.matrix(mds[rse_gene_dg$BrNum,])

#############################
## modeling for drugs #######
#############################

######### MDD ################
ind = rse_gene_dg$Dx %in% c("Control", "MDD")
rse_gene_mdd = rse_gene_dg[,ind]
rse_gene_mdd$Dx = droplevels(rse_gene_mdd$Dx)
rse_gene_mdd$MDD = NA
rse_gene_mdd$MDD[rse_gene_mdd$Dx == "Control"] = "Control"
rse_gene_mdd$MDD[rse_gene_mdd$Dx == "MDD" &  rse_gene_mdd$SSRI] = "MDD_SSRI"
rse_gene_mdd$MDD[rse_gene_mdd$Dx == "MDD" &  ! rse_gene_mdd$SSRI] = "MDD"

mod_dg_ssri = model.matrix(~ MDD + Age + totalAssignedGene + Sex + mitoRate + 
				mds_dg[ind,] + qSVs_dg[ind,],
				data = colData(rse_gene_mdd))
rse_gene_mdd = rse_gene_mdd[,rownames(mod_dg_ssri)]

dge_ssri <- DGEList(counts = assays(rse_gene_mdd)$counts, genes = rowData(rse_gene_mdd))
dge_ssri <- calcNormFactors(dge_ssri)
vGene_ssri <- voom(dge_ssri, mod_dg_ssri, plot = TRUE)
geneFit_ssri = lmFit(vGene_ssri, mod_dg_ssri)
geneFit_ssri = eBayes(geneFit_ssri)

ssri_effect = topTable(geneFit_ssri, coef = 3, n = nrow(dge_ssri), sort="none")
mdd_effect = topTable(geneFit_ssri, coef = 2, n = nrow(dge_ssri), sort="none")

sum(ssri_effect$adj.P.Val < 0.1)
sum(mdd_effect$adj.P.Val < 0.1)

ssri_sig = topTable(geneFit_ssri, coef = 3, n = nrow(dge_ssri), p.value = 0.1)
cat(ssri_sig$Symbol, sep = ", ")
names(geneDxStats_dg)[geneDxStats_dg$MD_sig] %in% rownames(ssri_sig)

plot(ssri_effect$t, mdd_effect$t)

######### SCZD ################
ind2 = rse_gene_dg$Dx %in% c("Control", "Schizo")
rse_gene_sz = rse_gene_dg[,ind2]
rse_gene_sz$Dx = droplevels(rse_gene_sz$Dx)
rse_gene_sz$SZ = NA
rse_gene_sz$SZ[rse_gene_sz$Dx == "Control"] = "Control"
rse_gene_sz$SZ[rse_gene_sz$Dx == "Schizo" &  rse_gene_sz$antipsychotics] = "SZCD_AP"
rse_gene_sz$SZ[rse_gene_sz$Dx == "Schizo" &  ! rse_gene_sz$antipsychotics] = "SCZD"
table(rse_gene_sz$SZ)

mod_dg_ap = model.matrix(~ SZ + Age + totalAssignedGene + Sex + mitoRate + 
				mds_dg[ind2,] + qSVs_dg[ind2,],
				data = colData(rse_gene_sz))
rse_gene_sz = rse_gene_sz[,rownames(mod_dg_ap)]

dge_ap <- DGEList(counts = assays(rse_gene_sz)$counts, genes = rowData(rse_gene_sz))
dge_ap <- calcNormFactors(dge_ap)
vGene_ap <- voom(dge_ap, mod_dg_ap, plot = TRUE)
geneFit_ap = lmFit(vGene_ap, mod_dg_ap)
geneFit_ap = eBayes(geneFit_ap)

ap_effect = topTable(geneFit_ap, coef = 3, n = nrow(dge_ap), sort="none")

ap_sig = topTable(geneFit_ap, coef = 3, n = nrow(dge_ap), p.value = 0.1)
noap_sig = topTable(geneFit_ap, coef = 2, n = nrow(dge_ap), p.value = 0.1)
cat(ap_sig$Symbol, sep = ", ")

table(names(geneDxStats_dg)[geneDxStats_dg$SZ_sig] %in% rownames(ap_sig))

save(ap_sig, ssri_sig, file = "rdas/de_drug_effects.rda")

######################
### rna fractions ####
######################

## load data
load("rdas/cell_type_fractions.rda")
cellPropEsts = cellPropEsts[colnames(rse_gene_joint),]
cellPropEstsScaled = prop.table(as.matrix(cellPropEsts),1)

## dx
rse_gene_joint$Dx = factor(rse_gene_joint$Dx, levels=c("Control", "Schizo", "MDD", "Bipolar"))
cellTypeListDx = lapply(as.data.frame(cellPropEstsScaled), function(y) {
	summary(lmer(y ~ Dx*Dataset + (1|BrNum), data=as.data.frame(colData(rse_gene_joint))))$coef
})
cellPvalMatDx = t(sapply(cellTypeListDx, function(x) x[,5]))[,-1]
write.csv(cellPvalMatDx, file = "tables/rna_fractions_across_dx_pvalues.csv")


###################################
###### neurogenesis overlap #######

# http://mango.adult-neurogenesis.de/
ng = read.delim("tables/MANGO_result_annotationLevel.txt",
	as.is=TRUE,header=FALSE)
colnames(ng) =c("ID", "Gene", "Process", "CellStage", "Effect","Evidence", "Species", "Type", "Ref")	

## try this data: https://www.genenames.org/tools/hcop/
orth = read.delim("tables/human_all_hcop_seven_column.txt.gz", as.is=TRUE)
orth = orth[grepl("ENSMUS", orth$ortholog_species_ensembl_gene),]

ng$inOrth = ng$ID %in% orth$ortholog_species_entrez_gene
ng$humanEnsembl = orth$human_ensembl_gene[match(ng$ID, orth$ortholog_species_entrez_gene)]


ap_sig$inMango = ap_sig$ensemblID %in% ng$humanEnsembl 
ssri_sig$inMango = ssri_sig$ensemblID %in% ng$humanEnsembl 

table(ap_sig$inMango)
ap_sig[which(ap_sig$inMango),]
table(ssri_sig$inMango)
ssri_sig[which(ssri_sig$inMango),]
#############
## expression
geneExprs = log2(getRPKM(rse_gene_joint, "Length")+1)
geneExprs = geneExprs[names(geneAgeStats),]

