###
library(GenomicRanges)
library(clusterProfiler)
library(jaffelab)
library(VennDiagram)
library(RColorBrewer)
library(limma)
library(lmerTest)

## load results
load("rdas/geneLevel_ageAndSzInteraction.rda")

## load DEqual plot stats
load("rdas/Hippocampus_geneLevel_degradationStats_forDEqual_hg38.rda")
degradeStats = degradeStats[names(geneAgeStats),]

##################
## focus on age ##
head(geneAgeStats)

## write out for CSV
out = mcols(geneAgeStats)
out = out[,c(14,17,1:12, 13,15,16,18,20)]
out = as.data.frame(out)
write.csv(out, file = "tables/ageStats_geneLevel_allTested.csv",row.names=FALSE)

# dg stats
sum(geneAgeStats$adj.P.Val_Age_DG < 0.05)
table(sign(geneAgeStats$t_Age_DG[geneAgeStats$adj.P.Val_Age_DG < 0.05]))
2^quantile(abs(geneAgeStats$logFC_Age_DG[geneAgeStats$adj.P.Val_Age_DG < 0.05])*10, c(.5, 0.25,0.75))

# hippo stats
sum(geneAgeStats$adj.P.Val_Age_Hippo < 0.05)
table(sign(geneAgeStats$t_Age_Hippo[geneAgeStats$adj.P.Val_Age_Hippo < 0.05]))
2^quantile(abs(geneAgeStats$logFC_Age_Hippo[geneAgeStats$adj.P.Val_Age_Hippo < 0.05])*10, c(.5, 0.25,0.75))

## overlap?
table(geneAgeStats$adj.P.Val_Age_DG < 0.05, 
	geneAgeStats$adj.P.Val_Age_Hippo < 0.05,
	dnn = c("DG", "HIPPO"))
getOR(table(geneAgeStats$adj.P.Val_Age_DG < 0.05, 
	geneAgeStats$adj.P.Val_Age_Hippo < 0.05,
	dnn = c("DG", "HIPPO")))

## same direction?
	
eitherIndex = which(geneAgeStats$adj.P.Val_Age_DG < 0.05 | geneAgeStats$adj.P.Val_Age_Hippo < 0.05)
bothIndex = which(geneAgeStats$adj.P.Val_Age_DG < 0.05 & geneAgeStats$adj.P.Val_Age_Hippo < 0.05)
length(eitherIndex)
tt = table(sign(geneAgeStats$t_Age_DG[eitherIndex]),sign(geneAgeStats$t_Age_Hippo[eitherIndex]) )
sum(diag(tt))
sum(diag(tt))/sum(tt)
(sum(tt)-sum(diag(tt)))
(sum(tt)-sum(diag(tt)))/sum(tt)
ttBoth = table(sign(geneAgeStats$t_Age_DG[bothIndex]),sign(geneAgeStats$t_Age_Hippo[bothIndex]) )

## colors
cols = rep(1, length(geneAgeStats))
cols[geneAgeStats$adj.P.Val_Age_DG > 0.05 & geneAgeStats$adj.P.Val_Age_Hippo < 0.05] = 2
cols[geneAgeStats$adj.P.Val_Age_DG < 0.05 & geneAgeStats$adj.P.Val_Age_Hippo > 0.05] = 3
cols[geneAgeStats$adj.P.Val_Age_DG < 0.05 & geneAgeStats$adj.P.Val_Age_Hippo < 0.05] = 4

names(cols)[cols==1] = "None"
names(cols)[cols==2] = "HIPPO"
names(cols)[cols==3] = "DG"
names(cols)[cols==4] = "Both"

pdf("plots/ageAnalysis_dgVsHippo_geneLevel.pdf")
par(mar=c(5,6,3,2),cex.axis=2,cex.lab=2)
palette(c("grey", brewer.pal(3, "Set1")[1:2], "purple"))
plot(geneAgeStats$t_Age_DG, geneAgeStats$t_Age_Hippo,
	ylim = c(-8, 8), xlim = c(-8, 8), pch = 21,bg=cols,
	xlab = "DG-GCL (Age T-stat)", ylab = "HIPPO (Age T-stat)", 
	main = "T-statistics")
legend("bottomright", c("HIPPO", "DG-GCL", "Both"),
	pch=15,col=2:4,cex=1.5)

plot(geneAgeStats$logFC_Age_DG, geneAgeStats$logFC_Age_Hippo,
	ylim = c(-0.04, 0.04), xlim = c(-0.04, 0.04), pch = 21,
	bg=cols, xlab = "DG-GCL", ylab = "HIPPO", main = "log2 Fold Changes")
legend("bottomright", c("HIPPO", "DG-GCL", "Both"),
	pch=15,col=2:4,cex=1.5)
dev.off()

table(abs(geneAgeStats$logFC_Age_DG) > abs(geneAgeStats$logFC_Age_Hippo))
table(abs(geneAgeStats$logFC_Age_DG[eitherIndex]) > abs(geneAgeStats$logFC_Age_Hippo[eitherIndex]))
table(abs(geneAgeStats$logFC_Age_DG[eitherIndex]) > abs(geneAgeStats$logFC_Age_Hippo[eitherIndex]))

## interaction
sum(geneAgeStats$adj.P.Val_Age_Inter < 0.05)
table(geneAgeStats$adj.P.Val_Age_Inter < 0.05,cols)
table(geneAgeStats$P.Value_Age_Inter < 0.05,cols)

# intIndex = which(geneAgeStats$adj.P.Val_Age_Inter < 0.05)
intIndex = which(geneAgeStats$adj.P.Val_Age_Inter < 0.05)
table(names(cols)[intIndex])

pdf("plots/ageAnalysis_dgVsHippo_geneLevel_interaction.pdf")
par(mar=c(5,6,3,2),cex.axis=2,cex.lab=2)
palette(c("grey", brewer.pal(3, "Set1")))
plot(geneAgeStats$t_Age_DG[intIndex], geneAgeStats$t_Age_Hippo[intIndex],
	ylim = c(-8, 8), xlim = c(-8, 8), pch = 21,bg=cols[intIndex],
	xlab = "DG", ylab = "HIPPO", main = "T-statistics")
plot(geneAgeStats$logFC_Age_DG[intIndex], geneAgeStats$logFC_Age_Hippo[intIndex],
	ylim = c(-0.04, 0.04), xlim = c(-0.04, 0.04), pch = 21,
	bg=cols[intIndex], xlab = "DG", ylab = "HIPPO", main = "log2 Fold Changes")
dev.off()

## degradation check
plot(geneAgeStats$t_Age_DG, degradeStats$t)
cor(geneAgeStats$logFC_Age_DG, degradeStats$logFC)
plot(geneAgeStats$logFC_Age_DG, degradeStats$logFC)
plot(geneAgeStats$t_Age_Hippo, degradeStats$t)


## venn diagram

gsPlot = apply(gsDatQ, 2, function(x) rownames(gsDatQ)[which(x < 0.05)])
names(gsPlot) = gsub("_", ":", names(gsPlot))
v = venn.diagram(gsPlot[1:4], 
	fill = c("darkblue", "darkgreen", "lightblue", "lightgreen"), 
	main="", main.pos = c(.5, .2), cat.cex = 1.8, cex=3,
	margin = 0.4, filename = NULL)
pdf("plots/age_geneSet_vennDiagram.pdf")
grid.draw(v)
dev.off()

#####################
### gene ontology ###
geneUniverse = geneAgeStats$EntrezID
geneUniverse = as.character(geneUniverse[!is.na(geneUniverse)])

## make sets
geneListsDg = split(geneAgeStats$EntrezID[cols %in% c(2,4)], 
	sign(geneAgeStats$t_Age_DG[cols %in% c(2,4)]))
names(geneListsDg) = c("DG_Down", "DG_Up")
geneListsHippo = split(geneAgeStats$EntrezID[cols %in% c(3,4)], 
	sign(geneAgeStats$t_Age_Hippo[cols %in% c(3,4)]))
names(geneListsHippo) = c("Hippo_Down", "Hippo_Up")
intList = split(geneAgeStats$EntrezID[geneAgeStats$adj.P.Val_Age_Inter < 0.05],
	paste0(sign(geneAgeStats$t_Age_DG),":",sign(geneAgeStats$t_Age_Hippo))[geneAgeStats$adj.P.Val_Age_Inter < 0.05])
names(intList) = c("Int_DdHd", "Int_DdHu", "Int_DuHd", "Int_DuHu")
geneList = c(geneListsDg, geneListsHippo,intList)
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
gsDatQ = sapply(gsDfList, function(x) x[match(uniqueIDs,x$ID), "qvalue"])
rownames(gsDat) = rownames(gsDatQ) = uniqueIDs
gsDat = as.data.frame(gsDat)
gsDatQ = as.data.frame(gsDatQ)
gsDat$Description = gsDf$Description[match(uniqueIDs, gsDf$ID)]
gsDatQ$Description = gsDf$Description[match(uniqueIDs, gsDf$ID)]
save(gsDat, gsDatQ, keggEnr, goEnr, file = "rdas/age_gene_set_results.rda")

#################
## make plots ###
#################
load( "rdas/age_gene_set_results.rda")
outDat = gsDat
outDat$SetID = rownames(outDat)
outDat = outDat[,c(10,1:9)]
write.table(outDat,  file = "tables/ageStats_Pvalues_geneSet_results.tsv",
	sep = "\t",row.names=FALSE)

colSums(gsDatQ[,-ncol(gsDatQ)] < 0.05,na.rm=TRUE)
table(rowSums(gsDatQ[,-ncol(gsDatQ)] < 0.05,na.rm=TRUE))
table(rowSums(gsDatQ[,1:2] < 0.05,na.rm=TRUE) > 0, 
	rowSums(gsDatQ[,3:4] < 0.05,na.rm=TRUE) > 0)
table(gsDatQ[,"DG_Down"] < 0.05, gsDatQ[,"Hippo_Down"] < 0.05,dnn=c("DG","Hippo"))
table(gsDatQ[,"DG_Up"] < 0.05, gsDatQ[,"Hippo_Up"] < 0.05,dnn=c("DG","Hippo"))

## venn diagram
library(limma)
vennDiagram(vennCounts(gsDatQ[,1:4] < 0.05))

## nicer
gsPlot = apply(gsDatQ, 2, function(x) rownames(gsDatQ)[which(x < 0.05)])
names(gsPlot) = gsub("_", ":", names(gsPlot))
v = venn.diagram(gsPlot[1:4], 
	fill = c("darkblue", "darkgreen", "lightblue", "lightgreen"), 
	main="", main.pos = c(.5, .2), cat.cex = 1.8, cex=3,
	margin = 0.4, filename = NULL)
pdf("plots/age_geneSet_vennDiagram.pdf")
grid.draw(v)
dev.off()

#################	
## examples #####

## down in DG
head(gsDatQ[which(gsDatQ$DG_Down < 0.05 & gsDatQ$Hippo_Down > 0.1),],20)
gsDatQ[which(gsDatQ$DG_Up < 0.05 & gsDatQ$Hippo_Up > 0.2),]

## up in hippo
head(gsDatQ[which(gsDatQ$Hippo_Up < 0.05 & gsDatQ$DG_Up > 0.2),],10)

## up in DG
head(gsDatQ[which(gsDatQ$DG_Up < 0.05 & gsDatQ$Hippo_Up > 0.2),],10)
head(gsDatQ[which(gsDatQ$DG_Up < 0.05 & gsDatQ$Hippo_Up < 0.05),],10)

#### heatmap
gsDat_plot = gsDat[c(which(gsDatQ$DG_Down < 0.05 & gsDatQ$Hippo_Down > 0.1)[1:8],
		which(gsDatQ$DG_Up < 0.05 & gsDatQ$Hippo_Up > 0.2)[1:8]),]
gsDat_plot[is.na(gsDat_plot)] = 1
for(i in 1:8) gsDat_plot[,i] = -log10(gsDat_plot[,i])
rownames(gsDat_plot) = gsDat_plot$Description
gsDat_plot$Description = NULL

## plot
library(lattice)
theSeq = seq(0,10.1,by=0.01) 
my.col <- colorRampPalette(c("white","red"))(length(theSeq))
pdf("plots/age_geneSet_heatmap_logPvalue.pdf",w=14)
print(levelplot(as.matrix(t(gsDat_plot[,1:4])),
	aspect = "fill", at = theSeq,pretty=T,
	panel = panel.levelplot.raster, col.regions = my.col,
	scales=list(x=list(cex=1.6), y=list(cex=1.4)),
	ylab = "", xlab = ""),
	colorkey = list(labels=list(cex=2)))
dev.off()		
		
################
names(gsPlot)[5:8] = c("DG:Down\nHippo:Up", "DG:Down\nHippo:Down", 
			"DG:Up\nHippo:Down", "DG:Up\nHippo:Up")
v2 = venn.diagram(gsPlot[5:8], 
	fill = c("brown", "pink", "purple", "grey"), 
	main="", main.pos = c(.5, .5), cat.cex = 1.7, cex=3,
	margin = 0.6, filename = NULL)
pdf("plots/age_geneSet_vennDiagram_interaction.pdf")
grid.draw(v2)
dev.off()

################
## interaction #
################

## down dg, up hippo
gsDatQ[which(gsDatQ$Int_DdHu < 0.05) ,]

## down dg, down hippo
gsDatQ[which(gsDatQ$Int_DdHd < 0.05)[1:20],]

## up dg, up hippo
gsDatQ[which(gsDatQ$Int_DuHu < 0.05)[1:20],]

#####
table(abs(geneAgeStats$logFC_Age_DG[intIndexMarg]) > 
	abs(geneAgeStats$logFC_Age_Hippo[intIndexMarg]))


hist(geneAgeStats$logFC_Age_DG[intIndexMarg] - geneAgeStats$logFC_Age_Hippo[intIndexMarg])

##########################
## plots some examples ###
##########################
library(SummarizedExperiment)
library(recount)
library(sva)
library(jaffelab)

## load counts
load("count_data/merged_dg_hippo_allSamples_n596.rda", verbose=TRUE)
rm(rse_exon_joint, rse_tx_joint, rse_jxn_joint)

#############
## expression
geneExprs = log2(getRPKM(rse_gene_joint, "Length")+1)
geneExprs = geneExprs[names(geneAgeStats),]

colData(rse_gene_joint)$Dataset = factor(colData(rse_gene_joint)$Dataset)
colData(rse_gene_joint)$Dataset = relevel(colData(rse_gene_joint)$Dataset , "Hippo")

## plots
exampleIndex = which(geneAgeStats$P.Value_Age_DG < 1e-10 & geneAgeStats$P.Value_Age_Hippo > 0.05)
rIndexes = splitit(rse_gene_joint$Dataset)

pdf("plots/ageEffect_example_genes.pdf")
par(mar=c(5,6,2,2), cex.axis=2,cex.lab=2, cex.main=2)
palette(brewer.pal(5,"Set1"))
for(i in exampleIndex) {
	plot(geneExprs[i,] ~ rse_gene_joint$Age,
		pch=21,bg = rse_gene_joint$Dataset,
		ylim = c(0, max(geneExprs[i,])),
		xlab= "Age", ylab = "Expression (log2)")
	for(j in seq(along=rIndexes)) {
		abline(lm(geneExprs[i,] ~ rse_gene_joint$Age, 
			subset=rIndexes[[j]]),col=j,lwd=5)
	}
	legend("topright", rowData(rse_gene_joint)$Symbol[i], bty="n",cex=2.5)
}
dev.off()

## plots
exampleIndexInt = which(geneAgeStats$adj.P.Val_Age_DG > 0.05 & 
	geneAgeStats$adj.P.Val_Age_Hippo > 0.05 & geneAgeStats$adj.P.Val_Age_Inter < 0.05)

pdf("plots/ageEffect_example_genes_interaction.pdf")
par(mar=c(5,6,2,2), cex.axis=2,cex.lab=2, cex.main=2)
palette(brewer.pal(5,"Set1"))
for(i in exampleIndexInt) {
	plot(geneExprs[i,] ~ rse_gene_joint$Age,
		pch=21,bg = rse_gene_joint$Dataset,
		xlab= "Age", ylab = "Expression (log2)")
	for(j in seq(along=rIndexes)) {
		abline(lm(geneExprs[i,] ~ rse_gene_joint$Age, 
			subset=rIndexes[[j]]),col=j,lwd=5)
	}
	legend("topright", rowData(rse_gene_joint)$Symbol[i], bty="n",cex=2.5)
}
dev.off()

###################
## rna fractions ##
###################

## load data
load("rdas/cell_type_fractions.rda")
cellPropEsts = cellPropEsts[colnames(rse_gene_joint),]
cellPropEstsScaled = prop.table(as.matrix(cellPropEsts),1)

## lmerTest
cellTypeList = lapply(as.data.frame(cellPropEstsScaled), function(y) {
	summary(lmer(y ~ Age*Dataset + (1|BrNum), data=as.data.frame(colData(rse_gene_joint))))$coef
})

cellPvalMat = t(sapply(cellTypeList, function(x) x[,5]))[,-1]
write.csv(cellPvalMat, file = "tables/rna_fractions_across_age_pvalues.csv")

###############################
## adult neurogenesis genes ###

# http://mango.adult-neurogenesis.de/
ng = read.delim("tables/MANGO_result_annotationLevel.txt",
	as.is=TRUE,header=FALSE)
colnames(ng) =c("ID", "Gene", "Process", "CellStage", "Effect","Evidence", "Species", "Type", "Ref")	

## match up
ng$inHuman = toupper(ng$Gene) %in% geneAgeStats$Symbol
table(ng$Species[!duplicated(ng$Gene)], ng$inHuman[!duplicated(ng$Gene)])

## try this data: https://www.genenames.org/tools/hcop/
orth = read.delim("tables/human_all_hcop_seven_column.txt.gz", as.is=TRUE)
orth = orth[grepl("ENSMUS", orth$ortholog_species_ensembl_gene),]

ng$inOrth = ng$ID %in% orth$ortholog_species_entrez_gene
ng$humanEnsembl = orth$human_ensembl_gene[match(ng$ID, orth$ortholog_species_entrez_gene)]
ng$exprsHuman = ng$humanEnsembl %in% geneAgeStats$ensemblID

table(ng$inOrth[!duplicated(ng$Gene)], ng$exprsHuman[!duplicated(ng$Gene)])

geneAgeStats$inMango = geneAgeStats$ensemblID %in% ng$humanEnsembl
geneAgeStats$matchMango = match(geneAgeStats$ensemblID, ng$humanEnsembl) 
geneAgeStats$dirMango = ng$Effect[geneAgeStats$matchMango]

table(geneAgeStats$dirMango)

## do gene set tests ##

## chisq
ttDg = table(geneAgeStats$adj.P.Val_Age_DG < 0.05, geneAgeStats$inMango)
getOR(ttDg)
chisq.test(ttDg) # 3.34e-6
ttHippo = table(geneAgeStats$adj.P.Val_Age_Hippo < 0.05, geneAgeStats$inMango)
getOR(ttHippo)
chisq.test(ttHippo) # 0.008
ttInt = table(geneAgeStats$adj.P.Val_Age_Inter < 0.05, geneAgeStats$inMango)
chisq.test(ttInt) # 0.06

## wilcox
geneSetTest(geneAgeStats$inMango, geneAgeStats$t_Age_DG, alternative="either")
geneSetTest(geneAgeStats$inMango, geneAgeStats$t_Age_DG, alternative="mixed")
geneSetTest(geneAgeStats$inMango, geneAgeStats$t_Age_DG, alternative="up")
geneSetTest(geneAgeStats$inMango, geneAgeStats$t_Age_DG, alternative="down")

geneSetTest(geneAgeStats$inMango, geneAgeStats$t_Age_Hippo, alternative="either")
geneSetTest(geneAgeStats$inMango, geneAgeStats$t_Age_Hippo, alternative="mixed")
geneSetTest(geneAgeStats$inMango, geneAgeStats$t_Age_Hippo, alternative="up")
geneSetTest(geneAgeStats$inMango, geneAgeStats$t_Age_Hippo, alternative="down")

pdf("plots/mango_enrichment.pdf")
par(mar=c(5,6,3,2), cex.axis=2,cex.lab=2,cex.main=2)
plot(density(geneAgeStats$t_Age_DG[!geneAgeStats$inMango]), col="black",lwd=3,
	xlab = "T-statistic", main = "DG - Age Effect", ylim = c(0,0.3))
lines(density(geneAgeStats$t_Age_DG[geneAgeStats$inMango]), col="red", lwd=3)
legend("topleft", c("MANGO", "Not"), col = c("red", "black"), pch = 15,cex=2)

plot(density(geneAgeStats$t_Age_Hippo[!geneAgeStats$inMango]), col="black",lwd=3,
	xlab = "T-statistic", main = "HIPPO - Age Effect", ylim = c(0,0.3))
lines(density(geneAgeStats$t_Age_Hippo[geneAgeStats$inMango]), col="red", lwd=3)
legend("topleft", c("MANGO", "Not"), col = c("red", "black"), pch = 15,cex=2)


indNeg = geneAgeStats$inMango & geneAgeStats$dirMango == "Negative"
indPos = geneAgeStats$inMango & geneAgeStats$dirMango == "Positive"
geneSetTest(indNeg, geneAgeStats$t_Age_DG, alternative = "down")
plot(density(geneAgeStats$t_Age_DG[!indNeg | !indNeg]), col="black",lwd=3,
	xlab = "T-statistic", main = "DG - Age Effect", ylim = c(0,0.3))
lines(density(geneAgeStats$t_Age_DG[indNeg]), col="blue", lwd=3)
lines(density(geneAgeStats$t_Age_DG[indPos]), col="green", lwd=3)
legend("topleft", c("MANGO-Pos","MANGO-Neg", "Not"), 
	col = c("green","blue", "black"), pch = 15,cex=1.5)

dev.off()

geneSetTest(indNeg, geneAgeStats$t_Age_Hippo, alternative = "down")
geneSetTest(indPos , geneAgeStats$t_Age_Hippo, alternative = "up")

## examples
xx = geneAgeStats[indNeg,]
xx[order(xx$t_Age_DG),]
#############################
## other features ###########
#############################
load("rdas/exonLevel_ageAndSzInteraction.rda")
load("rdas/jxLevel_ageAndSzInteraction.rda")
load("rdas/txLevel_ageAndSzInteraction.rda")
rm(geneSzStats, exonSzStats, jxSzStats, txSzStats)
