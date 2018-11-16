##
library(jaffelab)
library(IRanges)
library(SummarizedExperiment)


#####################
##### Subset of 881 SNPs from PGC
#####################

################
## load eQTLs

load("eqtl_tables/mergedEqtl_output_dg_raggr_4features.rda", verbose=TRUE)
dg = allEqtl[allEqtl$FDR < 0.05,]
load("eqtl_tables/mergedEqtl_output_hippo_raggr_4features.rda", verbose=TRUE)
hippo = allEqtl[allEqtl$FDR < 0.05,]

# load("eqtl_tables/mergedEqtl_output_interaction_4features.rda", verbose=TRUE)
# inter = allEqtl[allEqtl$FDR < 0.05,]

################
## metrics

## total features
nrow(dg)  ## 68176
nrow(hippo)   ## 54737
# nrow(inter) ## 

## per feature
table(dg$Type)
# Exon   Gene  Jxn   Tx
# 39907  6430 11542 10297
table(hippo$Type)
# Exon 	Gene Jxn Tx
# 31274  4188  9828  9447

# table(inter$Type)
# # Exon  Gene Jxn  Tx
# # 

## unique ensemblIDs
tapply(dg$EnsemblGeneID, dg$Type, function(x) length(unique(x)))
# Exon Gene  Jxn   Tx
#  280  172  175  233
tapply(hippo$EnsemblGeneID, hippo$Type, function(x) length(unique(x)))
# Exon Gene  Jxn   Tx
#   282  135  169  221

# tapply(inter$EnsemblGeneID, inter$Type, function(x) length(unique(x)))
# # Exon Gene  Jxn   Tx
# #   


################
## make csv

# hippo$EnsemblGeneID = ss(hippo$EnsemblGeneID, "\\.")
# dg$EnsemblGeneID = ss(dg$EnsemblGeneID, "\\.")

## snpMap
load("../genotype_data/astellas_dg_genotype_data_n263.rda")
snpMap1 = snpMap
snpMap1$hg19POS = paste0(snpMap1$CHR,":",snpMap1$POS)
snpMap1 = snpMap1[which(rownames(snpMap1) %in% c(hippo$snps,dg$snps) ),c("SNP","chr_hg38","pos_hg38","hg19POS")]

load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/genotype_data/BrainSeq_Phase2_RiboZero_Genotypes_n551.rda")
snpMap2 = snpMap
snpMap2$hg19POS = paste0(snpMap2$CHR,":",snpMap2$POS)
snpMap2 = snpMap2[which(rownames(snpMap2) %in% c(hippo$snps,dg$snps) ),c("SNP","chr_hg38","pos_hg38","hg19POS")]

snpMap = snpMap1[snpMap1$hg19POS %in% snpMap2$hg19POS,]

## featMap
load("/dcl01/ajaffe/data/lab/dg_hippo/count_data/merged_dg_hippo_allSamples_n596.rda", verbose=TRUE)

gMap = as.data.frame(rowRanges(rse_gene_joint))[,c("seqnames","start","end","strand","Class")]
eMap = as.data.frame(rowRanges(rse_exon_joint))[,c("seqnames","start","end","strand","Class")]
jMap = as.data.frame(rowRanges(rse_jxn_joint))[,c("seqnames","start","end","strand","Class")]
txMap = as.data.frame(rowRanges(rse_tx_joint))[,c("seqnames","start","end","strand","source")]
txMap$source = "InGen"

colnames(gMap) = colnames(eMap) = colnames(jMap) = colnames(txMap) = 
	c("feat_chr","feat_start","feat_end","strand","Class")
featMap = rbind(rbind(rbind(gMap, eMap),jMap),txMap)
featMap$Type = c(rep("Gene",nrow(gMap)),rep("Exon",nrow(eMap)),rep("Jxn",nrow(jMap)),rep("Tx",nrow(txMap)))

geneMap = as.data.frame(rowRanges(rse_gene_joint))[,c("gencodeID","Symbol","ensemblID","gene_type")]

## put together

## hippo
snpMap_temp = snpMap[hippo$snps,]
featMap_temp = featMap[hippo$gene,]
geneMap_temp = geneMap[match(hippo$EnsemblGeneID, geneMap$ensemblID),]
hippo2 = cbind(snpMap_temp,featMap_temp,geneMap_temp,hippo)
hippo3 = hippo2[,c(1:4,16,10,5:9,21:22,14,17:20)]
write.csv(hippo3, "raggr_179_snps_hippo_eqtls_fdr05.csv")

## DG
snpMap_temp = snpMap[dg$snps,]
featMap_temp = featMap[dg$gene,]
geneMap_temp = geneMap[match(dg$EnsemblGeneID, geneMap$ensemblID),]
dg2 = cbind(snpMap_temp,featMap_temp,geneMap_temp,dg)
dg3 = dg2[,c(1:4,16,10,5:9,21:22,14,17:20)]
write.csv(dg3, "raggr_179_snps_dg_eqtls_fdr05.csv")


# ## interaction
# snpMap_temp = snpMap[inter$snps,]
# featMap_temp = featMap[inter$gene,]
# geneMap_temp = geneMap[match(inter$EnsemblGeneID, geneMap$ensemblID),]
# inter2 = cbind(snpMap_temp,featMap_temp,geneMap_temp,inter)
# inter3 = inter2[,c(1:4,16,10,5:9,21:22,14,17:20)]
# write.csv(inter3, "raggr_179_snps_inter_eqtls_fdr05.csv")













