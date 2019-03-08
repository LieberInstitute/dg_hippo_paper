##
library(GenomicRanges)
library(SummarizedExperiment)
library(jaffelab)
library(recount)

## load hippo vs dlpfc
load("rdas/DE_output_cellType_lmer.rda")

## check age effects
load("rdas/geneLevel_ageAndSzInteraction.rda")

## load merged gene counts
load("count_data/dgPlusHippo_hg38_rseGene_n224.rda")

colData(rse_gene_joint)$Dataset = factor(colData(rse_gene_joint)$Dataset)
colData(rse_gene_joint)$Dataset = relevel(colData(rse_gene_joint)$Dataset , "Hippo")
geneIndex = rowMeans(getRPKM(rse_gene_joint[,rse_gene_joint$Dataset == "DG"], "Length")) > 0.5
sum(geneIndex)


## marker genes
g = c("SOX1", "SOX2", "DCX", "MKI67", "NCAM1")
cat(g, sep=", ")

all(g %in% rowData(rse_gene_joint)$Symbol)
g[g %in% outGene$Symbol]
g[ ! g %in% outGene$Symbol]

## check if any reads
geneRpkm = getRPKM(rse_gene_joint, "Length")
geneRpkmSub = geneRpkm[match(g, rowData(rse_gene_joint)$Symbol),]
rIndexes = splitit(rse_gene_joint$Dataset)
gMeans = signif(as.data.frame(sapply(rIndexes, function(ii) rowMeans(geneRpkmSub[,ii]))),3)
gMeans$Symbol = g
gMeans

## filter current results
mm1 = match(rownames(gMeans), rownames(outGene))
mm1 = mm1[!is.na(mm1)]
outGene[mm1,]

## filter current results
mm2 = match(rownames(gMeans), names(geneAgeStats))
mm2 = mm2[!is.na(mm2)]
geneAgeStats[mm2]