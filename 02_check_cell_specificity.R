###

## load libraries
library(jaffelab)
library(SummarizedExperiment)
library(edgeR)
library(recount)
library(genefilter)
library(clusterProfiler)
library(org.Hs.eg.db)
library(readxl)
library(RColorBrewer)
library(TeachingDemos)
library(lmerTest)

## make tables
dir.create("plots")
dir.create("tables")

## load data
load("count_data/dgPlusHippo_hg38_rseGene_n224.rda")

## make factor
colData(rse_gene_joint)$Dataset = ifelse(colData(rse_gene_joint)$Dataset == "DG", "DG-GCL", "HIPPO")
colData(rse_gene_joint)$Dataset = factor(colData(rse_gene_joint)$Dataset)
colData(rse_gene_joint)$Dataset = relevel(colData(rse_gene_joint)$Dataset , "HIPPO")

## filter on RPKM level
geneIndex = rowMeans(getRPKM(rse_gene_joint[,rse_gene_joint$Dataset == "DG-GCL"], "Length")) > 0.5
rse_gene_joint = rse_gene_joint[geneIndex,]

##### DGE ######
dge = DGEList(counts = assays(rse_gene_joint)$counts, 
	genes = rowData(rse_gene_joint))
dge = calcNormFactors(dge)

####################################
## start w/ gene level

## modeling
mod = model.matrix(~Dataset + mitoRate + rRNA_rate + overallMapRate + totalAssignedGene , 
	data=colData(rse_gene_joint))

## mean-variance
vGene = voom(dge,mod,plot=FALSE)
# gene_dupCorr = duplicateCorrelation(vGene$E, mod, block=rse_gene_joint$BrNum)
# save(gene_dupCorr, file = "rdas/geneLevel_duplicateCorrelation_regionVsCell.rda")
load("rdas/geneLevel_duplicateCorrelation_regionVsCell.rda")

## do analysis
fitGene = lmFit(vGene, block=rse_gene_joint$BrNum, 
	correlation=gene_dupCorr$consensus.correlation)

## top table
eBGene = eBayes(fitGene)
outGene = topTable(eBGene,coef=2,
	p.value = 1,number=nrow(rse_gene_joint), 
	adjust.method = "bonferroni")
outGene$sigColor = as.numeric(outGene$adj.P.Val < 0.01)+1
sigGene = outGene[outGene$adj.P.Val < 0.01,]

# for abstract
dim(sigGene)
max(sigGene$P.Value)

## write out
save(outGene, file = "rdas/DE_output_cellType_lmer.rda")

## subset
theGenes = c("PROX1", "KCNK1", "SOX9","MBP", "MOBP","CAMK1","GABRD")
gStats = outGene[match(theGenes, outGene$Symbol),]

## filter
sigGeneDf = as.data.frame(sigGene)
sigGeneDf$gencodeTx = sapply(sigGeneDf$gencodeTx, paste, collapse=";")

write.csv(sigGeneDf[,c(2,5, 9,11:14,10, 1, 3:4,6)], quote=FALSE,
	row.names=FALSE, file="tables/DE_output_cellType_lmer.csv")

####################
##### plots ########
bIndexes = splitit(rse_gene_joint$BrNum)

## pca
pca = prcomp(t(vGene$E))
pcaVars = getPcaVars(pca)

pdf("plots/pcaPlot_enrichment_n224_paired.pdf")
par(mar=c(5,6,2,2),cex.axis=2,cex.lab=2)
palette(brewer.pal(5,"Set1"))
plot(pca$x, pch = 21, bg = rse_gene_joint$Dataset,cex=1.5,
	xlab=paste0("PC1: ", pcaVars[1], "% Var Expl"),
	ylab=paste0("PC2: ", pcaVars[2], "% Var Expl"))
for(j in seq(along=bIndexes)) {
	ii = bIndexes[[j]]
	lines(pca$x[ii,1], pca$x[ii,2], col ="grey",lwd=0.4)
}
# legend("bottomright", levels(rse_gene_joint$Dataset), 
	# col=1:2,pch=15,cex=2,nc=2)
dev.off()

### target genes #######
pdf("plots/markerGene_enrichment.pdf")
exprs = vGene$E[match(theGenes, vGene$genes$Symbol),]
palette(brewer.pal(5,"Set1"))
par(mar=c(5,6,3,2), cex.axis=2,cex.lab=2,cex.main=3)
for(i in seq(along=theGenes)) {
	boxplot(exprs[i,] ~ rse_gene_joint$Dataset, 
		xlab = "", outline=FALSE,	ylab = "Normalized Expression",
		main = theGenes[i], ylim = c(0,max(exprs[i,])))
	xx = jitter(as.numeric(rse_gene_joint$Dataset),amount=0.1)
	for(j in seq(along=bIndexes)) {
		lines(exprs[i,] ~ xx, data=colData(rse_gene_joint),
			subset=bIndexes[[j]], col ="grey",lwd=0.4)
	}
	points(exprs[i,] ~ xx,	pch = 21, bg = rse_gene_joint$Dataset)
	pv = paste0("p=",signif( gStats$P.Value[i],3))
	legend("bottom", pv, cex=1.6)
}		
dev.off()

up = 2^gStats$logFC  
names(up) = theGenes
up

down = 1/ 2^gStats$logFC  
names(down) = theGenes
down

## volano
theGenes2 = c("PROX1", "SOX9","MBP")
m2 = match(theGenes2, outGene$Symbol)

pdf("plots/volanoPlot_cellType.pdf")
palette(brewer.pal(5, "Dark2"))
par(mar=c(5,6,3,2), cex.axis=2,cex.lab=2,cex.main=2)
plot(-log10(P.Value) ~ logFC, pch = 21, bg=sigColor,
	cex=0.8,data = outGene, xlab = "DG-GCL vs HIPPO log2FC")
shadowtext(outGene$logFC[m2]+0.35, -log10(outGene$P.Value[m2]),
	LETTERS[3:5],font=2,cex=1.25,col="grey")
legend("top", c("FDR>0.05", "FDR<0.05"), 
	pch =15, col = 1:2,cex=1.5)
dev.off()

####################
## gene set 
## split by sign
sigGeneSplit = split(sigGene, sign(sigGene$logFC))
geneListSplit = sapply(sigGeneSplit, function(x) as.character(x$EntrezID[!is.na(x$EntrezID)]))
names(geneListSplit) = c("HIPPO", "DG-GCL")

geneUniverse = as.character(outGene$EntrezID)
geneUniverse = geneUniverse[!is.na(geneUniverse)]

goMF <- compareCluster(geneListSplit, universe = geneUniverse,
				fun = "enrichGO", ont = "MF", 
				OrgDb = org.Hs.eg.db, pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1,
				readable= TRUE)
goBP <- compareCluster(geneListSplit, universe = geneUniverse,
				fun = "enrichGO", ont = "BP", 
				OrgDb = org.Hs.eg.db, pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1,
				readable= TRUE)	
goCC <- compareCluster(geneListSplit, universe = geneUniverse,
				fun = "enrichGO", ont = "CC", 
				OrgDb = org.Hs.eg.db, pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1,
				readable= TRUE)

goMF_sig <- compareCluster(geneListSplit, universe = geneUniverse,
				fun = "enrichGO", ont = "MF", 
				OrgDb = org.Hs.eg.db, pAdjustMethod = "BH",
                pvalueCutoff  = .1, qvalueCutoff  = 0.05,
				readable= TRUE)
goBP_sig <- compareCluster(geneListSplit, universe = geneUniverse,
				fun = "enrichGO", ont = "BP", 
				OrgDb = org.Hs.eg.db, pAdjustMethod = "BH",
                pvalueCutoff  = .1, qvalueCutoff  = 0.05,
				readable= TRUE)
goCC_sig <- compareCluster(geneListSplit, universe = geneUniverse,
				fun = "enrichGO", ont = "CC", 
				OrgDb = org.Hs.eg.db, pAdjustMethod = "BH",
                pvalueCutoff  = .1, qvalueCutoff  = 0.05,
				readable= TRUE)
				
save(goMF,goBP,goCC, goMF_sig, goBP_sig, goCC_sig,
	file = "rdas/geneSetEnrichment_Robjs_lmer.rda")


pdf("plots/geneSetEnrichment_cellType.pdf",
	useDingbats=FALSE, h=8,w=7)
dotplot(goMF_sig, showCategory=10)				
dotplot(goBP_sig, showCategory=10)				
dotplot(goCC_sig, showCategory=10)				
dev.off()

## write out

goOut = rbind(as.data.frame(goBP), as.data.frame(goMF),
	as.data.frame(goCC))
goOut$geneID = NULL
goOut$Ontology = rep(c("BP", "MF","CC"), 
	c(nrow(goBP), nrow(goMF), nrow(goCC)))
goOut$Cluster = relevel(goOut$Cluster, "DG-GCL")
goOut = goOut[order(goOut$Cluster, goOut$pvalue),]	
write.table(goOut, file="tables/geneSetEnrichment_suppTable_lmer.tsv",
	sep = "\t", row.names=FALSE, quote=TRUE)

###############################
## rna fraction enrichments ###
###############################

load("rdas/cell_type_fractions.rda")

cellPropEsts = cellPropEsts[colnames(rse_gene_joint),]
cellPropEstsScaled = prop.table(as.matrix(cellPropEsts),1)

## lmerTest
cellTypeDiffs = t(apply(cellPropEstsScaled, 2, function(y) {
	summary(lmer(y ~ Dataset + (1|BrNum), data=as.data.frame(colData(rse_gene_joint))))$coef[2,c(1,4,5)]
}))
colnames(cellTypeDiffs) = c("PropDiff", "statistic", "P.Value")
cellTypeDiffs = as.data.frame(cellTypeDiffs)

## make plots
cn = c("Excitatory Neurons", "Inhibitory Neurons", "Mixed Neurons",
	"Neural Progenitors", "Oligodendrocytes", "OPCs", "Astrocytes", "Microglia")
	
pdf("plots/cell_fractions_by_dataset.pdf",width=5)
par(mar=c(5,6,2,2), cex.axis=2,cex.lab=2,cex.main=2)
palette(brewer.pal(5,"Set1"))
for(i in 1:ncol(cellPropEstsScaled)) {
	boxplot(cellPropEstsScaled[,i] ~ rse_gene_joint$Dataset, 
			outline=FALSE,	ylab = "RNA Fraction",
			main = cn[i], ylim = c(0,1))
	xx = jitter(as.numeric(rse_gene_joint$Dataset),amount=0.1)
	for(j in seq(along=bIndexes)) {
		lines(cellPropEstsScaled[,i] ~ xx, data=colData(rse_gene_joint),
			subset=bIndexes[[j]], col ="grey",lwd=0.4)
	}
	points(cellPropEstsScaled[,i] ~ xx,	pch = 21, bg = rse_gene_joint$Dataset)
	# ll = ifelse(cellTypeDiffs$PropDiff[i] > 0, "topleft", "topright")
	pv = paste0("p=",signif( cellTypeDiffs$P.Value[i],3))
	legend("topright", pv, cex=1.6)
}
dev.off()

####################
## drug gene sets ##
####################

## drugable genome
drugGenome = read.csv("tables/drug_dbs/IDG_TargetList_20190124.csv",as.is=TRUE)
drugGenome$IDG.Family = factor(drugGenome$IDG.Family, 
	levels = c("GPCR","IC","Kinase"))
drugGeneList = split(drugGenome$HGNC.Symbol, drugGenome$IDG.Family)
names(drugGeneList) = paste0(names(drugGeneList) , "_IDG")

## other sets
minGene = 50
dsigdb = read.delim("tables/drug_dbs/DSigDB_All_detailed.txt",as.is=TRUE)
drugSigList = split(dsigdb$Gene, dsigdb$Drug)
names(drugSigList) = gsub(" ", "-", names(drugSigList))
drugSigList = drugSigList[lengths(drugSigList) > minGene]
names(drugSigList) = paste0(names(drugSigList),"_dsigdb")

dgidb = read.delim("tables/drug_dbs/DGIdb_interactions.tsv", as.is=TRUE)
drugDgiList = split(dgidb$gene_name, dgidb$drug_claim_name)
drugDgiList = drugDgiList[lengths(drugDgiList) > minGene]
names(drugDgiList) = paste0(gsub(" ", "-", names(drugDgiList)), "_DGI")

## overlap to ours
drugInSet = sapply(c(drugGeneList, drugSigList,drugDgiList), function(x) outGene$Symbol %in% x)
rownames(drugInSet) = rownames(outGene)
drugInSet = as.data.frame(drugInSet)

drugTabList = lapply(drugInSet, function(x) tt = table(x, outGene$sigColor))
## stats
chisqList = lapply(drugTabList, chisq.test)

## extract
drugEnr = data.frame(OR = sapply(drugTabList, getOR),
				pval = sapply(chisqList, function(x) x$p.value))
drugEnr$FDR = p.adjust(drugEnr$pval, "fdr")
drugEnr$Bonf = p.adjust(drugEnr$pval, "bonf")
drugEnr$Drug = ss(names(drugTabList), "_")
drugEnr$Database = ss(names(drugTabList), "_",2)
drugEnr$numGenesSet = colSums(drugInSet)
drugEnr$numGenesSig = sapply(drugTabList, function(x) x[2,2])
drugEnr = drugEnr[order(drugEnr$pval),]
write.csv(drugEnr, file = "tables/drug_enrichment_dggclVsHippo.csv",row.names=FALSE)

head(drugEnr[drugEnr$numGenesSet < 1000,],10)


#####################
##### MANGO data ####
#####################

# http://mango.adult-neurogenesis.de/
ng = read.delim("tables/MANGO_result_annotationLevel.txt",
	as.is=TRUE,header=FALSE)
colnames(ng) =c("ID", "Gene", "Process", "CellStage", "Effect","Evidence", "Species", "Type", "Ref")	

## try this data: https://www.genenames.org/tools/hcop/
orth = read.delim("tables/human_all_hcop_seven_column.txt.gz", as.is=TRUE)
orth = orth[grepl("ENSMUS", orth$ortholog_species_ensembl_gene),]

ng$inOrth = ng$ID %in% orth$ortholog_species_entrez_gene
ng$humanEnsembl = orth$human_ensembl_gene[match(ng$ID, orth$ortholog_species_entrez_gene)]
ng$exprsHuman = ng$humanEnsembl %in% outGene$ensemblID

table(ng$inOrth[!duplicated(ng$Gene)], ng$exprsHuman[!duplicated(ng$Gene)])

## add effects
outGene$inMango = outGene$ensemblID %in% ng$humanEnsembl
outGene$matchMango = match(outGene$ensemblID, ng$humanEnsembl) 
outGene$dirMango = ng$Effect[outGene$matchMango]

## chisq
tt = table(outGene$adj.P.Val < 0.05, outGene$inMango)
getOR(tt) # 3.245
chisq.test(tt) # 1.79e-12

with(outGene[outGene$adj.P.Val < 0.05 & outGene$inMango,], 
	table(dirMango, t > 0, dnn = c("Dir", "DG_up")))

## overall
plot(density(outGene$t[!outGene$inMango]), col="black",lwd=3,
	xlab = "T-statistic", main = "DG vs Hippo Effect", ylim = c(0,0.15))
lines(density(outGene$t[outGene$inMango]), col="red", lwd=3)
legend("topleft", c("MANGO", "Not"), col = c("red", "black"), pch = 15,cex=2)

## by direction
indNeg = outGene$inMango & outGene$dirMango == "Negative"
indPos = outGene$inMango & outGene$dirMango == "Positive"
geneSetTest(indNeg, outGene$t, alternative = "down")
plot(density(outGene$t[!indNeg | !indNeg]), col="black",lwd=3,
	xlab = "(HIPPO > DG) T-statistic (HIPPO < DG)", main = "", 
	xlim = c(-20,20), ylim = c(0,0.15))
lines(density(outGene$t[indNeg]), col="blue", lwd=3)
lines(density(outGene$t[indPos]), col="green", lwd=3)
legend("topleft", c("MANGO-Pos","MANGO-Neg", "Not"), 
	col = c("green","blue", "black"), pch = 15,cex=1.5)
	
## wilcox
geneSetTest(outGene$inMango, outGene$t, alternative="either")
geneSetTest(indNeg, outGene$t, alternative="down")
geneSetTest(indPos, outGene$t, alternative="up")

geneSetTest(outGene$inMango, outGene$t, alternative="mixed")
geneSetTest(outGene$inMango, outGene$t, alternative="up")
geneSetTest(outGene$inMango, outGene$t, alternative="down")

