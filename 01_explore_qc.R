######

## load libraries
library(jaffelab)
library(SummarizedExperiment)
library(readxl)
library(RColorBrewer)

## load DG
load("count_data/astellas_dg_hg38_rseGene_n263.rda")
rse_gene_dg = rse_gene
pdDg = colData(rse_gene_dg)
pdDg$Dataset = "DG-GCL"
for(i in grep("integer", sapply(pdDg,class))) pdDg[,i] = as.numeric(pdDg[,i])
pdDg$trimmed = as.character(pdDg$trimmed) # to match
pdDg$Protocol = "RiboZeroGold"

## hippo
load("count_data/hippo_brainseq_phase2_hg38_rseGene_merged_n447.rda")
rse_gene = merge_rse_metrics(rse_gene) # from jaffelab

## add kit info
hipxl <- read_excel('count_data/LIBD_PhaseII_HIPPO_RiboZero_sample_list_01_28_2015.xlsx')
rse_gene$Protocol = hipxl$Protocol[match(rse_gene$RNum, paste0("R", hipxl$RNum))]

pdH = colData(rse_gene)
pdH$Dataset = "Hippo"

## fix character lists
classH = sapply(pdH,class)
for(i in grep("CharacterList", classH)) pdH[,i] = sapply(pdH[,i], paste, collapse=";")
for(i in grep("LogicalList", classH)) pdH[,i] = sapply(pdH[,i], paste, collapse=";")
for(i in grep("NumericList", classH)) pdH[,i] = sapply(pdH[,i], mean)
for(i in grep("IntegerList", classH)) pdH[,i] = as.numeric(sapply(pdH[,i], mean))

## make phenotype columns merge-able 
n = intersect(colnames(pdDg), colnames(pdH))
colData(rse_gene_dg) = pdDg[,n]
colData(rse_gene) = pdH[,n]

### match up
mm = match(pdDg$BrNum, pdH$BrNum)
table(!is.na(mm)) # 129

mcols(rse_gene)$meanExprs = mcols(rse_gene_dg)$meanExprs = NULL
rse_gene_joint = cbind(rse_gene_dg[,!is.na(mm)], rse_gene[,mm[!is.na(mm)]])

rse_gene_joint$Dataset = factor(rse_gene_joint$Dataset, 
	levels = c("Hippo", "DG-GCL"))

############################
## explore, this becomes a figure for QC
bIndexes = splitit(rse_gene_joint$BrNum)

dir.create("qcChecks")
pdf("qcChecks/quality_by_cellType.pdf")
par(mar=c(5,6,3,2), cex.axis=1.5,cex.lab=2,cex.main=2)
palette(brewer.pal(5,"Set2"))
boxplot(mitoRate ~ Dataset, data=colData(rse_gene_joint),
	ylab = "chrM Map Rate", outline = FALSE, 
	ylim = range(rse_gene_joint$mitoRate)) # ribozero regular...hippo > dg
xx = jitter(as.numeric(rse_gene_joint$Dataset),amount=0.15)
for(i in seq(along=bIndexes)) {
	lines(mitoRate ~ xx, data=colData(rse_gene_joint),
		subset=bIndexes[[i]], col ="grey",lwd=0.4)
}
points(	mitoRate ~ xx,	data=colData(rse_gene_joint),
	bg = factor(rse_gene_joint$Protocol), pch =21)
legend("topright", c("Gold", "HMR"), col = 1:2, pch = 15,cex=2,nc=1)
boxplot(rRNA_rate ~ Dataset, data=colData(rse_gene_joint),
	ylab = "rRNA Gene Map Rate",outline=FALSE,
	ylim = range(rse_gene_joint$rRNA_rate)) # ribozero regular...hippo > dg
for(i in seq(along=bIndexes)) {
	lines(rRNA_rate ~ xx, data=colData(rse_gene_joint),
		subset=bIndexes[[i]], col ="grey",lwd=0.4)
}
points(	rRNA_rate ~ xx,	data=colData(rse_gene_joint),
	bg = factor(rse_gene_joint$Protocol), pch =21)

boxplot(RIN ~ Dataset, data=colData(rse_gene_joint), ylab = "RIN",
	outline=FALSE, ylim = range(rse_gene_joint$RIN, na.rm=TRUE)) # hippo > dg
for(i in seq(along=bIndexes)) {
	lines(RIN ~ xx, data=colData(rse_gene_joint),
		subset=bIndexes[[i]], col ="grey",lwd=0.4)
}
points(	RIN ~ xx,	data=colData(rse_gene_joint),
	bg = factor(rse_gene_joint$Protocol), pch =21)


boxplot(overallMapRate ~ Dataset, data=colData(rse_gene_joint),
	ylab = "Overall Map Rate",
	outline=FALSE, ylim = range(rse_gene_joint$overallMapRate)) # hippo > dg
for(i in seq(along=bIndexes)) {
	lines(overallMapRate ~ xx, data=colData(rse_gene_joint),
		subset=bIndexes[[i]], col ="grey",lwd=0.4)
}
points(	overallMapRate ~ xx,	data=colData(rse_gene_joint),
	bg = factor(rse_gene_joint$Protocol), pch =21)

	
boxplot(totalAssignedGene ~ Dataset, data=colData(rse_gene_joint),
	ylab = "Exonic Mapping Rate", outline=FALSE, 
	ylim = range(rse_gene_joint$totalAssignedGene)) # hippo > dg
for(i in seq(along=bIndexes)) {
	lines(totalAssignedGene ~ xx, data=colData(rse_gene_joint),
		subset=bIndexes[[i]], col ="grey",lwd=0.4)
}
points(	totalAssignedGene ~ xx,	data=colData(rse_gene_joint),
	bg = factor(rse_gene_joint$Protocol), pch =21)

dev.off()

