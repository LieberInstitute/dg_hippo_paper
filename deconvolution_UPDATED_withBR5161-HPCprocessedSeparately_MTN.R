### Deconvolute DG-GCL (LCM-seq) / bulk HPC RNA-seq with LIBD snRNA-seq
###   - using `/dcl01/lieber/ajaffe/lab/libd_stem_timecourse/deconvolution/cell_type/check_darmanis_cellType.R`
###     for reference
###   - Generate 'geneCountsSum' object with expression collapsed by cluster/cell type
### MTN 24Aug2019
########################################################################

library(SummarizedExperiment)
library(jaffelab)
library(recount)
library(Seurat)
library(genefilter)
library(RColorBrewer)
library(dplyr)
library(rtracklayer)

## Load 'pilot.br5161.hpc' Seurat object already processed
load("count_data/pilot.br5161.hpc.Rdata")
    
pilot.clusterNamesNums <- data.frame(clusters=c(0:12), clusterNames=
                                       c("Oligo","Oligo","Micro","Astro","OPC","Excit","Astro","Inhib","Inhib",
                                         "mixedNeu","Inhib","NPC","Inhib"),
                                     numberedNames=c("Oligo.1","Oligo.2","Micro","Astro.1","OPC","Excit","Astro.2","Inhib.1","Inhib.2",
                                                     "mixedNeu","Inhib.3","NPC","Inhib.4"))


pilot.br5161.hpc[["cellClusterID"]] <- pilot.clusterNamesNums$clusterNames[match(pilot.br5161.hpc@meta.data$seurat_cluster, pilot.clusterNamesNums$clusters)]
pilot.br5161.hpc[["cellClusterID.numbered"]] <- pilot.clusterNamesNums$numberedNames[match(pilot.br5161.hpc@meta.data$seurat_cluster, pilot.clusterNamesNums$clusters)]

pilotExprs <- as.matrix(pilot.br5161.hpc@assays$RNA@counts)

gCS.pilot.index <- splitit(pilot.br5161.hpc@meta.data$cellClusterID)
# Count sums only
gCS.pilot.byClusterType <- sapply(gCS.pilot.index, function(ii) rowSums(pilotExprs[ ,ii]))


## Now generate RSE/RPKM
    # jump to line 124 - replace 'pilotExprs' with 'gCS.pilot.byClusterType'
load("count_data/astellas_dg_hg38_rseGene_n263.rda")

## update symbols
gencodeGTF = import(con="/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh38_hg38/gencode.v25.annotationGRCh38.gtf", format="gtf")
gencodeGENES = mcols(gencodeGTF)[which(gencodeGTF$type=="gene"),c("gene_id","type","gene_type","gene_name")]
rownames(gencodeGENES) = gencodeGENES$gene_id
rowData(rse_gene)$Symbol = gencodeGENES[rownames(rse_gene),"gene_name"]

## how many genes
table(rownames(gCS.pilot.byClusterType) %in% gr_genes$Symbol)  # 16843

gCS.pilot.byClusterType <- gCS.pilot.byClusterType[which(rownames(gCS.pilot.byClusterType) %in% gr_genes$Symbol), ]
#rownames(gCS.pilot.byClusterType) <- gr_genes$gencodeID[match(rownames(gCS.pilot.byClusterType), gr_genes$Symbol)]
    ## Will keep gene symbol cause will make subsetting on training probes simpler

#gr_genes <- gr_genes[rownames(gCS.pilot.byClusterType), ]
    ## Weirdly this doesn't work anymore... (gives NULL)

length(gr_genes)  # 58037
gr_genes <- gr_genes[match(rownames(gCS.pilot.byClusterType), gr_genes$Symbol), ]
length(gr_genes)  # 16843
# Check order
table(gr_genes$Symbol==rownames(gCS.pilot.byClusterType)) # all TRUE - good.
names(gr_genes) <- gr_genes$Symbol

rse_gene.hpc10x <- SummarizedExperiment(assays = list('counts' = gCS.pilot.byClusterType),
                                        rowRanges = gr_genes)

## Make rpkm object and subset for training probes - do 20 per cluster first
dim(rse_gene.hpc10x)  # 16843 x 8
yExprs.hpc10x = log2(recount::getRPKM(rse_gene.hpc10x, "Length") + 1)
#rownames(yExprs.hpc10x) <- rowRanges(rse_gene.hpc10x)$Symbol

## Now bring in cluster-driving markers from individual-HPC-sample analysis
load("rdas/pilot.markers.br5161.hpc.Rdata")

trainingProbes.20perCluster = pilot.markers.br5161.hpc %>% group_by(cluster) %>%
  top_n(n = 20, wt = avg_logFC) %>% as.data.frame  # of nrow 260 duh

length(unique(trainingProbes.20perCluster$gene))  #[1] 218
trainingProbes.20perCluster <- unique(trainingProbes.20perCluster$gene)

table(trainingProbes.20perCluster %in% rownames(yExprs.hpc10x)) # 216/218
keepProbes <- trainingProbes.20perCluster[trainingProbes.20perCluster %in% rownames(yExprs.hpc10x)]
yExprs.hpc10x <- yExprs.hpc10x[keepProbes, ]

yExprs.hpc10x <- yExprs.hpc10x[ ,c("Excit", "Inhib", "mixedNeu", "NPC",
                                   "Oligo", "OPC", "Astro", "Micro")]

# `minfi::validationCellType()` still needs a formula and identity matrix, even though
#     not modeling at the cell level
form <- as.formula(sprintf("y ~ %s - 1", paste(colnames(yExprs.hpc10x),collapse = "+")))
phenoDF <- as.data.frame(model.matrix(~colnames(yExprs.hpc10x) - 1))
colnames(phenoDF) <- ss(colnames(phenoDF), ".hpc10x)", 2)


## Scale, then do calibration
yExprs.hpc10x.Z <- scale(yExprs.hpc10x)
coefEsts.hpc10x.groupCollapsed <-minfi:::validationCellType(Y = yExprs.hpc10x.Z,
                                                             pheno=phenoDF, modelFix=form)$coefEsts

rownames(coefEsts.hpc10x.groupCollapsed) <- gr_genes$gencodeID[match(rownames(coefEsts.hpc10x.groupCollapsed), gr_genes$Symbol)]


### Deconvolution ====================================
## Load DG & bulk-HPC data
# DG-GCL LCM-seq
## dim 58037 x 263
rse_dg.full <- rse_gene
rm(rse_gene)

# HPC bulk RNA-seq
load("count_data/hippo_brainseq_phase2_hg38_rseGene_merged_n447.rda")
## 'rse_gene' of dim 58037 x 447
rse_hpc.full <- rse_gene
rm(rse_gene)

# Subset for 333 age-matched samples in manuscript
load("count_data/merged_dg_hippo_allSamples_n596.rda", verbose=TRUE)
rm(rse_jxn_joint, rse_tx_joint, rse_exon_joint, rse_deg_joint)

rse_hpc.full <- rse_hpc.full[ ,colnames(rse_gene_joint)[which(rse_gene_joint$Region=="HIPPO")]]

table(rownames(coefEsts.hpc10x.groupCollapsed) %in% rownames(rse_hpc.full))  # all there


### DG-GCL LCM-seq ========================
geneExprs.dg <- log2(recount::getRPKM(rse_dg.full, "Length")+1)
table(rownames(coefEsts.hpc10x.groupCollapsed) %in% rownames(geneExprs.dg))  # all (216) TRUE

geneExprs_scale.dg <- scale(geneExprs.dg[rownames(coefEsts.hpc10x.groupCollapsed)[rownames(coefEsts.hpc10x.groupCollapsed) %in% rownames(geneExprs.dg)], ])

cellPropEsts.dg <- minfi:::projectCellType(geneExprs_scale.dg, coefEsts.hpc10x.groupCollapsed)
  
rse_dg.full$AgeDecade <- paste0(substr(rse_dg.full$Age,1,1),"0's")
  
cellPropEsts.dg_scaled <- prop.table(cellPropEsts.dg,1)
gIndexes = splitit(factor(rse_dg.full$AgeDecade))
cellPropEsts.dg_groupMeans = sapply(gIndexes, function(ii) colMeans(cellPropEsts.dg[ii,]))
cellPropEsts.dg_groupSEs = sapply(gIndexes, function(ii) apply(cellPropEsts.dg[ii,],2,function(x) sd(x)/sqrt(length(x))))
 
# Bar plot
pal = brewer.pal(10, "Set3")

pdf("plots/cellTypeDecon_barplot_DG-GCL-LCM-seq_usingBR5161-HPC-geneCountSums_MTN24Aug2019.pdf",w=20,h=8)
palette(pal)
par(mar = c(14,6,2,2), cex.axis=1.5,cex.lab=2)
oo= order(factor(rse_dg.full$AgeDecade), factor(rse_dg.full$Dx))
bp = barplot(t(cellPropEsts.dg_scaled[oo,]), col = pal,
             ylim = c(0,1.45),xaxt = "n", yaxt= "n",ylab="Class Proportion            ")
g = split(bp, factor(rse_dg.full$AgeDecade)[oo])

axis(1,at=sapply(g,mean), font=2, labels = names(g),las=3, cex.axis=2)
abline(v = sapply(g,min)[-1]-0.6,lwd=3.5,col="black")

legend("top", colnames(cellPropEsts.dg_scaled), pch = 15, 
       col = 1:10,bg="white", nc = 5, cex=2, pt.cex=3)
axis(2, at=seq(0,1,by=0.25))
# Annotate with Dx identifiers
text(x=unname(unlist(g)),y=1.03, labels=substr(rse_dg.full$Dx[oo],1,1), cex=0.6, font=2)
dev.off()



## HPC bulk-RNA-seq ========================
geneExprs.hpc <- log2(recount::getRPKM(rse_hpc.full, "Length")+1)
table(rownames(coefEsts.hpc10x.groupCollapsed) %in% rownames(geneExprs.hpc))  # all TRUE

geneExprs_scale.hpc <- scale(geneExprs.hpc[rownames(coefEsts.hpc10x.groupCollapsed)[rownames(coefEsts.hpc10x.groupCollapsed) %in% rownames(geneExprs.hpc)], ])
cellPropEsts.hpc = minfi:::projectCellType(geneExprs_scale.hpc, coefEsts.hpc10x.groupCollapsed)

rse_hpc.full$AgeDecade <- paste0(substr(rse_hpc.full$Age,1,1),"0's")

cellPropEsts.hpc_scaled <- prop.table(cellPropEsts.hpc,1)
gIndexes = splitit(factor(rse_hpc.full$AgeDecade))
cellPropEsts.hpc_groupMeans = sapply(gIndexes, function(ii) colMeans(cellPropEsts.hpc[ii,]))
cellPropEsts.hpc_groupSEs = sapply(gIndexes, function(ii) apply(cellPropEsts.hpc[ii,],2,function(x) sd(x)/sqrt(length(x))))

### Bar plot
pdf("plots/cellTypeDecon_barplot_HPC-RNA-seq_usingBR5161-HPC-geneCountSums_MTN24Aug2019.pdf",w=20,h=8)
palette(pal)
par(mar = c(14,6,2,2), cex.axis=1.5,cex.lab=2)
oo= order(factor(rse_hpc.full$AgeDecade), factor(rse_hpc.full$Dx))
bp = barplot(t(cellPropEsts.hpc_scaled[oo,]), col = pal,
             ylim = c(0,1.45),xaxt = "n", yaxt= "n",ylab="Class Proportion            ")
g = split(bp, factor(rse_hpc.full$AgeDecade)[oo])

axis(1,at=sapply(g,mean), font=2, labels = names(g),las=3, cex.axis=2)
abline(v = sapply(g,min)[-1]-0.6,lwd=3.5,col="black")

legend("top", colnames(cellPropEsts.hpc_scaled), pch = 15, 
       col = 1:10,bg="white", nc = 5, cex=2, pt.cex=3)
axis(2, at=seq(0,1,by=0.25))
# Annotate with Dx identifiers
text(x=unname(unlist(g)),y=1.03, labels=substr(rse_hpc.full$Dx[oo],1,1), cex=0.6, font=2)
dev.off()
  

## Go ahead and save these stats:
save(coefEsts.hpc10x.groupCollapsed, yExprs.hpc10x.Z, file = "rdas/BR5161-HPC-only_UPDATED_coefEsts_calibration_Zscale_MTN24Aug2019.rda")
write.csv(coefEsts.hpc10x.groupCollapsed, file = "tables/BR5161-HPC-only_UPDATED_coefEsts_calibration_Zscale_MTN24Aug2019.csv")

## save proportions
cellPropEsts = as.data.frame(rbind(cellPropEsts.dg, cellPropEsts.hpc))
save(cellPropEsts, file = "rdas/cell_type_fractions.rda")




