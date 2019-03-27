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

## make tables
dir.create("plots")
dir.create("tables")

## load data
load("count_data/dgPlusHippo_hg38_rseGene_n224.rda")

## make factor
colData(rse_gene_joint)$Dataset = factor(colData(rse_gene_joint)$Dataset)
colData(rse_gene_joint)$Dataset = relevel(colData(rse_gene_joint)$Dataset , "Hippo")

## filter on RPKM level
geneIndex = rowMeans(getRPKM(rse_gene_joint[,rse_gene_joint$Dataset == "DG"], "Length")) > 0.5
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
vGene = voom(dge,mod,plot=TRUE)
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

pdf("plots/pcaPlot_enrichment.pdf")
par(mar=c(5,6,2,2),cex.axis=2,cex.lab=2)
palette(brewer.pal(5,"Set1"))
plot(pca$x, pch = 21, bg = rse_gene_joint$Dataset,cex=1.5,
	xlab=paste0("PC1: ", pcaVars[1], "% Var Expl"),
	ylab=paste0("PC2: ", pcaVars[2], "% Var Expl"))
legend("bottomright", levels(rse_gene_joint$Dataset), 
	col=1:2,pch=15,cex=2,nc=2)
dev.off()

### target genes #######
theGenes = c("PROX1", "KCNK1", "GFAP","MBP", "MOBP","CAMK1","GABRD")
pdf("plots/markerGene_enrichment.pdf")
palette(brewer.pal(5,"Set1"))
par(mar=c(5,6,3,2), cex.axis=2,cex.lab=2,cex.main=2)
for(i in seq(along=theGenes)) {
	matchInd = match(theGenes[i], vGene$genes$Symbol)
	boxplot(vGene$E[matchInd,] ~ rse_gene_joint$Dataset, 
			outline=FALSE,	ylab = "Normalized Expression",
			main = theGenes[i], ylim = range(vGene$E[matchInd,]))
	xx = jitter(as.numeric(rse_gene_joint$Dataset),amount=0.1)
	for(j in seq(along=bIndexes)) {
		lines(vGene$E[matchInd,] ~ xx, data=colData(rse_gene_joint),
			subset=bIndexes[[j]], col ="grey",lwd=0.4)
	}
	points(vGene$E[matchInd,] ~ xx,	pch = 21, bg = rse_gene_joint$Dataset)
	ll = ifelse(fitGene$coef[matchInd,2] > 0, "topleft", "topright")
	pv = paste0("p ", ifelse(outGene$P.Value[matchInd] < 1e-20, "< 1e-20", 
		paste0("= ", signif( outGene$P.Value[matchInd],3))))
	legend(ll, pv, cex=1.6)

}		
dev.off()

2^fitGene$coef[match(theGenes,vGene$genes$Symbol),2] 
1/(2^fitGene$coef[match(theGenes,vGene$genes$Symbol),2] )

## volano
pdf("plots/volanoPlot_cellType.pdf")
palette(brewer.pal(5, "Dark2"))
par(mar=c(5,6,3,2), cex.axis=2,cex.lab=2,cex.main=2)
plot(-log10(P.Value) ~ logFC, pch = 21, bg=sigColor,
	cex=0.8,data = outGene, xlab = "DG vs HIPPO log2FC")
dev.off()

####################
## gene set 
## split by sign
sigGeneSplit = split(sigGene, sign(sigGene$logFC))
geneListSplit = sapply(sigGeneSplit, function(x) as.character(x$EntrezID[!is.na(x$EntrezID)]))
names(geneListSplit) = c("HIPPO", "DG")

geneUniverse = as.character(outGene$EntrezID)
geneUniverse = geneUniverse[!is.na(geneUniverse)]

goEnr <- compareCluster(geneListSplit, universe = geneUniverse,
				fun = "enrichGO", ont = "ALL", 
				OrgDb = org.Hs.eg.db, pAdjustMethod = "BH",
                pvalueCutoff  = 0.1, qvalueCutoff  = 0.05,
				readable= TRUE)
goDf = as.data.frame(goEnr)
save(goEnr,goDf, file = "rdas/geneSetEnrichment_Robjs_lmer.rda")

pdf("plots/geneSetEnrichment_cellType.pdf",
	useDingbats=FALSE, h=8,w=7)
dotplot(goEnr, showCategory=10)				
dotplot(goMF, showCategory=10)				
dotplot(goCC, showCategory=10)				
dev.off()

## write out

goOut = rbind(as.data.frame(goBP), as.data.frame(goMF),
	as.data.frame(goCC))
goOut$geneID = NULL
goOut$Ontology = rep(c("BP", "MF","CC"), 
	c(nrow(goBP), nrow(goMF), nrow(goCC)))
goOut$Cluster = relevel(goOut$Cluster, "DG")
goOut = goOut[order(goOut$Cluster, goOut$pvalue),]	
write.csv(goOut, file="tables/geneSetEnrichment_suppTable_lmer.csv",
	row.names=FALSE, quote=FALSE)

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
