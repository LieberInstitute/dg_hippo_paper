###
library(GenomicRanges)
library(clusterProfiler)
library(limma)

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

#################
## make plots ###
#################
load( "rdas/age_gene_set_results.rda")

gsDat[which(gsDat$DG_Down < 1e-4 & gsDat$Hippo_Down > 0.05),]
gsDat[which(gsDat$DG_Up < 1e-4 & gsDat$Hippo_Up > 0.05),]
gsDat[which(gsDat$Hippo_Down < 1e-4 & gsDat$DG_Down > 0.05),]
gsDat[which(gsDat$Hippo_Up < 1e-4 & gsDat$DG_Up > 0.05),]




table(abs(geneSzStats$logFC_Sz_DG[intIndexMarg]) > 
	abs(geneSzStats$logFC_Sz_Hippo[intIndexMarg]))


hist(geneSzStats$logFC_Sz_DG[intIndexMarg] - geneSzStats$logFC_Sz_Hippo[intIndexMarg])
load("rdas/exonLevel_ageAndSzInteraction.rda")
load("rdas/jxLevel_ageAndSzInteraction.rda")
load("rdas/txLevel_ageAndSzInteraction.rda")
rm(geneAgeStats, exonAgeStats, jxAgeStats, txAgeStats)
