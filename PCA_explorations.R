####### libraries #############
library(SummarizedExperiment)
library(jaffelab)
library(MatrixEQTL)
library(sva)
library(limma)
library(edgeR)
library(recount)
library(RColorBrewer)

######################
## load data #########
######################

## filtered counts
load("count_data/merged_dg_hippo_allSamples_n596.rda")
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

colData(rse_gene_joint)$Dataset = ifelse(colData(rse_gene_joint)$Dataset == "DG", "DG-GCL", "HIPPO")
colData(rse_gene_joint)$Dataset = factor(colData(rse_gene_joint)$Dataset)
colData(rse_gene_joint)$Dataset = relevel(colData(rse_gene_joint)$Dataset , "HIPPO")


#################
## get RPKMs ####
#################

geneRpkm = getRPKM(rse_gene_joint, "Length")
geneExprs = log2(geneRpkm+1)

#################
# overall PCA ###
#################
pca = prcomp(t(geneExprs))
pcaVars = getPcaVars(pca)

pdf("plots/pcaPlot_allSamples_n596.pdf")
par(mar=c(5,6,2,2),cex.axis=2,cex.lab=2)
palette(brewer.pal(5,"Set1"))
plot(pca$x, pch = 21, bg = rse_gene_joint$Dataset,cex=1.5,
	xlab=paste0("PC1: ", pcaVars[1], "% Var Expl"),
	ylab=paste0("PC2: ", pcaVars[2], "% Var Expl"))
legend("topright", levels(rse_gene_joint$Dataset), 
	col=1:2,pch=15,cex=2,nc=2)
dev.off()

###################
## by tissues #####
###################

# DG PCA
geneExprs_dg = geneExprs[,rse_gene_joint$Dataset == "DG-GCL"]
pca_dg= prcomp(t(geneExprs_dg))
pcaVars_dg = getPcaVars(pca_dg)

pdf("plots/pcaPlot_dggclSamples_n263.pdf")
par(mar=c(5,6,2,2),cex.axis=2,cex.lab=2)
palette(brewer.pal(5,"Set1"))
plot(pca_dg$x, pch = 21, bg = 2,cex=1.5,
	xlab=paste0("PC1: ", pcaVars_dg[1], "% Var Expl"),
	ylab=paste0("PC2: ", pcaVars_dg[2], "% Var Expl"))
dev.off()

# HIPPO PCA
geneExprs_h = geneExprs[,rse_gene_joint$Dataset == "HIPPO"]
pca_h= prcomp(t(geneExprs_h))
pcaVars_h = getPcaVars(pca_h)

pdf("plots/pcaPlot_hippoSamples_n333.pdf")
par(mar=c(5,6,2,2),cex.axis=2,cex.lab=2)
palette(brewer.pal(5,"Set1"))
plot(pca_h$x, pch = 21, bg = 1,cex=1.5,
	xlab=paste0("PC1: ", pcaVars_h[1], "% Var Expl"),
	ylab=paste0("PC2: ", pcaVars_h[2], "% Var Expl"))
dev.off()

###################
## projections ####
###################

## dg based
geneExprs_h_scaled = scale(t(geneExprs_h), pca_dg$center, pca_dg$scale) 
genePCs_h_projected_on_dg = geneExprs_h_scaled %*% pca_dg$rotation 
xlim_dg = range(c(pca_dg$x[,1], genePCs_h_projected_on_dg[,1]))
ylim_dg = range(c(pca_dg$x[,2], genePCs_h_projected_on_dg[,2]))

pdf("plots/pcaPlot_dggclSamples_n263_withHippoProject.pdf")
par(mar=c(5,6,2,2),cex.axis=2,cex.lab=2)
palette(brewer.pal(5,"Set1"))
plot(genePCs_h_projected_on_dg, pch = 22, bg = 1,cex=1,
	xlim = xlim_dg, ylim = ylim_dg,
	xlab=paste0("PC1: ", pcaVars_dg[1], "% Var Expl"),
	ylab=paste0("PC2: ", pcaVars_dg[2], "% Var Expl"))
points(pca_dg$x, pch = 21, bg=2,cex=1.5)
dev.off()

## hippo based
geneExprs_dg_scaled = scale(t(geneExprs_dg), pca_h$center, pca_h$scale) 
genePCs_dg_projected_on_h = geneExprs_dg_scaled %*% pca_h$rotation 
xlim_h = range(c(pca_h$x[,1], genePCs_dg_projected_on_h[,1]))
ylim_h = range(c(pca_h$x[,2], genePCs_dg_projected_on_h[,2]))

pdf("plots/pcaPlot_hippoSamples_n333_withDggclProject.pdf")
par(mar=c(5,6,2,2),cex.axis=2,cex.lab=2)
palette(brewer.pal(5,"Set1"))
plot(genePCs_dg_projected_on_h, pch = 22, bg = 2,cex=1,
	xlim = xlim_h, ylim = ylim_h,
	xlab=paste0("PC1: ", pcaVars_dg[1], "% Var Expl"),
	ylab=paste0("PC2: ", pcaVars_dg[2], "% Var Expl"))
points(pca_h$x, pch = 21, bg=1,cex=1.5)
dev.off()

