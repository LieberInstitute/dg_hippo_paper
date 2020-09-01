library(limma)
library(jaffelab)
library(VennDiagram)
library(RColorBrewer)

## read in results
load("rdas/twas_gene.Rdata",verbose=TRUE)
load("rdas/twas_exp_gene.Rdata",verbose=TRUE)
load("rdas/tt_objects_gene.Rdata",verbose=TRUE)

tt = as.data.frame(tt) 
write.csv(tt, file = "suppTables/TWAS_allRegions_geneLevel_noH2filter.csv",row.names=FALSE)

## split by region
ttList = split(tt, tt$region)
names(ttList)[1] = "DG-GCL"

## put in same order
g = unique(tt$geneid)
ttList = lapply(ttList, function(x) {
	x = as.data.frame(x[match(g, x$geneid),])
	rownames(x) = g
	x})

## get out matrices
tMat = sapply(ttList, "[[", "TWAS.Z")
pMat = sapply(ttList, "[[", "TWAS.P")
pMat[is.na(pMat)] = 1
fdrMat = sapply(ttList, "[[", "TWAS.FDR")
fdrMat[is.na(fdrMat)] = 1
bonfMat = sapply(ttList, "[[", "TWAS.Bonf")
bonfMat[is.na(bonfMat)] = 1
rownames(tMat) = rownames(pMat) = rownames(fdrMat) = rownames(bonfMat) = g

## tstats
pdf("plots/suppFigure7_hippoVsDGvsDlpfc_TWASz.pdf")
par(mar=c(5,6,2,2), cex.axis=2,cex.lab=2,cex.main=2)
pairs(tMat, pch =21,bg="grey",main = "TWAS Z")
dev.off()

cor(tMat, use="pair")

#########################
## hippo vs dg plot

palette(c("grey", brewer.pal(3, "Set1")))
par(mar=c(5,6,2,2), cex.axis=2,cex.lab=2)
plot(tMat[,"HIPPO"], tMat[,"DG-GCL"], xlab="HIPPO TWAS Z",
	ylab="DG-GCL TWAS Z",pch=21,bg="grey",
	ylim = c(-8,8),xlim=c(-8,8))
plot(tMat[,"DLPFC"], tMat[,"DG-GCL"], xlab="DLPFC TWAS Z",
	ylab="DG-GCL TWAS Z",pch=21,bg="grey",
	ylim = c(-8,8),xlim=c(-8,8))
plot(tMat[,"DLPFC"], tMat[,"HIPPO"], xlab="DLPFC TWAS Z",
	ylab="HIPPO TWAS Z",pch=21,bg="grey",
	ylim = c(-8,8),xlim=c(-8,8))
dev.off()

## fdr
colSums(fdrMat < 0.05)
vennDiagram(vennCounts(fdrMat < 0.05))

fdrList = lapply(ttList, function(x) x[which(x$TWAS.FDR < 0.05),])
v = venn.diagram(lapply(fdrList, "[[", "geneid"), 
	fill = c("darksalmon", "lightgoldenrod", "lightsteelblue2"), main="", main.pos = c(.5, .05), cat.cex = 2.5, cex=3,
	margin = .1, filename = NULL)
pdf("plots/venn_diagram_TWAS_Fdr05.pdf")
grid.draw(v)
dev.off()


## bonf 
colSums(bonfMat < 0.1)
vennDiagram(vennCounts(bonfMat < 0.1))

bonfList = lapply(ttList, function(x) x[which(x$TWAS.Bonf < 0.1),])
v = venn.diagram(lapply(bonfList, "[[", "geneid"), 
	fill = c("darksalmon", "lightgoldenrod", "lightsteelblue2"), main="", main.pos = c(.5, .05), cat.cex = 2.5, cex=3,
	margin = .1, filename = NULL)
pdf("plots/venn_diagram_TWAS_Bonf10.pdf")
grid.draw(v)
dev.off()

#### anchor on dg
dg_fdr_ind = which(fdrMat[,"DG-GCL"] < 0.05)
colMeans(fdrMat[dg_fdr_ind,]  < 0.05)
colMeans(pMat[dg_fdr_ind,]  < 0.05)
colMeans(pMat[dg_fdr_ind,]  < 0.05 & 
	sign(tMat[dg_fdr_ind,]) != sign(tMat[dg_fdr_ind,1]),na.rm=TRUE)

## check 2
tt[tt$genesymbol %in% c("CACNA1C" ,"GRM3"),]

#################
## DE analysis ##
#################
load("rdas/geneLevel_dxEffects_dg.rda")
tt_dg = ttList$`DG-GCL`

tt_dg$SCZD_t = geneDxStats_dg$SZ_t[match(tt_dg$geneid, names(geneDxStats_dg))]
tt_dg$SCZD_pvalue  = geneDxStats_dg$SZ_P.Value[match(tt_dg$geneid, names(geneDxStats_dg))]
tt_dg$SCZD_FDR = geneDxStats_dg$SZ_adj.P.Val[match(tt_dg$geneid, names(geneDxStats_dg))]

sum(!is.na(tt_dg$SCZD_t) & !is.na(tt_dg$TWAS.Z)) 
# 12256

pdf("plots/suppFigure_twasVsDe_moreExprs.pdf")
par(mar=c(5,6,2,2), cex.axis=2,cex.lab=2)
plot(tt_dg$SCZD_t, tt_dg$TWAS.Z,xlab="SCZD DE T-stat",
	ylab="SCZD TWAS Z-score",pch=21,bg="grey",
	ylim = c(-8,8),xlim=c(-8,8))
dev.off()


###########
## drug ###
###########

## drugable genome
drugGenome = read.csv("tables/drug_dbs/IDG_TargetList_20190124.csv",as.is=TRUE)
drugGenome$IDG.Family = factor(drugGenome$IDG.Family, 
	levels = c("GPCR","IC","Kinase"))
drugGeneList = split(drugGenome$HGNC.Symbol, drugGenome$IDG.Family)
names(drugGeneList) = paste0(names(drugGeneList) , "@IDG")

## other sets
dsigdb = read.delim("tables/drug_dbs/DSigDB_All_detailed.txt",as.is=TRUE)
drugSigList = split(dsigdb$Gene, dsigdb$Drug)
names(drugSigList) = gsub(" ", "-", names(drugSigList))
names(drugSigList) = paste0(names(drugSigList),"@dsigdb")

dgidb = read.delim("tables/drug_dbs/DGIdb_interactions.tsv", as.is=TRUE)
drugDgiList = split(dgidb$gene_name, dgidb$drug_claim_name)
names(drugDgiList) = paste0(gsub(" ", "-", names(drugDgiList)), "@DGI")

## overlap to ours
sym = ttList$DLPFC$genesymbol
drugInSet = sapply(c(drugGeneList, drugSigList,drugDgiList), 
	function(x) sym %in% x)
rownames(drugInSet) = rownames(ttList$DLPFC)

## filter to 50 genes per set
drugInSet = drugInSet[,colSums(drugInSet) > 50]
drugInSet = as.data.frame(drugInSet)

## cross tabulate
drugTabList_DG = lapply(drugInSet, function(x) tt = table(x, ttList$`DG-GCL`$TWAS.FDR < 0.05))
drugTabList_DLPFC = lapply(drugInSet, function(x) tt = table(x, ttList$DLPFC$TWAS.FDR < 0.05))
drugTabList_HIPPO= lapply(drugInSet, function(x) tt = table(x, ttList$HIPPO$TWAS.FDR < 0.05))

## stats
chisqList_DG = lapply(drugTabList_DG, chisq.test)
chisqList_DLPFC = lapply(drugTabList_DLPFC, chisq.test)
chisqList_HIPPO = lapply(drugTabList_HIPPO, chisq.test)

## extract
drugEnr = data.frame(OR_DG = sapply(drugTabList_DG, getOR),
		OR_DLPFC = sapply(drugTabList_DLPFC, getOR),
		OR_HIPPO = sapply(drugTabList_HIPPO, getOR),
		pval_DG = sapply(chisqList_DG, "[[", "p.value"),
		pval_DLPFC = sapply(chisqList_DLPFC,  "[[", "p.value"),
		pval_HIPPO = sapply(chisqList_HIPPO,  "[[", "p.value"))
	
drugEnr$FDR_DG = p.adjust(drugEnr$pval_DG, "fdr")
drugEnr$FDR_DLPFC = p.adjust(drugEnr$pval_DLPFC, "fdr")
drugEnr$FDR_HIPPO = p.adjust(drugEnr$pval_HIPPO, "fdr")

drugEnr$Drug = ss(names(drugTabList_DG), "@")
drugEnr$Database = ss(names(drugTabList_DG), "@",2)

drugEnr$numGenesSet = colSums(drugInSet)

drugEnr$numGenesSig_DG = sapply(drugTabList_DG, function(x) x[2,2])
drugEnr$numGenesSig_DLPFC = sapply(drugTabList_DLPFC, function(x) x[2,2])
drugEnr$numGenesSig_HIPPO = sapply(drugTabList_HIPPO, function(x) x[2,2])

drugEnr = drugEnr[order(drugEnr$pval_DG),]

write.csv(drugEnr, file = "tables/drug_enrichment_twas_Fdr05.csv",row.names=FALSE)
