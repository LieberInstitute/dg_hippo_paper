###

library(GenomicRanges)

## get DG eQTLs
load("/dcl01/ajaffe/data/lab/dg_hippo/eQTL_GWAS_riskSNPs_n596/eqtl_tables/mergedEqtl_output_dg_raggr_4features.rda",verbose=TRUE)
allEqtl$ID = paste0(allEqtl$snps, ";", allEqtl$gene)

########################
## get hippo eQTLs
load("/dcl01/ajaffe/data/lab/dg_hippo/eQTL_GWAS_riskSNPs_n596/eqtl_tables/matrixEqtl_output_hippo_raggr_4features.rda")

# extract
geneEqtlH = meGene$cis$eqtls
geneEqtlH$gene = as.character(geneEqtlH$gene)
geneEqtlH$snps = as.character(geneEqtlH$snps)
geneEqtlH$ID = paste0(geneEqtlH$snps, ";", geneEqtlH$gene)
geneEqtlH = geneEqtlH[geneEqtlH$ID %in% allEqtl$ID,]

exonEqtlH = meExon$cis$eqtls
exonEqtlH$gene = as.character(exonEqtlH$gene)
exonEqtlH$snps = as.character(exonEqtlH$snps)
exonEqtlH$ID = paste0(exonEqtlH$snps, ";", exonEqtlH$gene)
exonEqtlH = exonEqtlH[exonEqtlH$ID %in% allEqtl$ID,]

jxnEqtlH = meJxn$cis$eqtls
jxnEqtlH$gene = as.character(jxnEqtlH$gene)
jxnEqtlH$snps = as.character(jxnEqtlH$snps)
jxnEqtlH$ID = paste0(jxnEqtlH$snps, ";", jxnEqtlH$gene)
jxnEqtlH = jxnEqtlH[jxnEqtlH$ID %in% allEqtl$ID,]

txEqtlH = meTx$cis$eqtls
txEqtlH$gene = as.character(txEqtlH$gene)
txEqtlH$snps = as.character(txEqtlH$snps)
txEqtlH$ID = paste0(txEqtlH$snps, ";", txEqtlH$gene)
txEqtlH = txEqtlH[txEqtlH$ID %in% allEqtl$ID,]

## join
allEqtlH = rbind(geneEqtlH, exonEqtlH, jxnEqtlH, txEqtlH)
hippoStats = allEqtlH[match(allEqtl$ID, allEqtlH$ID),]
hippoStats = hippoStats[,3:6]
colnames(hippoStats) = paste0("hippo_", colnames(hippoStats))

allEqtl = cbind(allEqtl, hippoStats)

##################################
### interaction ##################
##################################

load("/dcl01/ajaffe/data/lab/dg_hippo/eQTL_GWAS_riskSNPs_n596/eqtl_tables/matrixEqtl_output_interaction_4features.rda")

# extract
geneEqtlInt = meGene$cis$eqtls
geneEqtlInt$gene = as.character(geneEqtlInt$gene)
geneEqtlInt$snps = as.character(geneEqtlInt$snps)
geneEqtlInt$ID = paste0(geneEqtlInt$snps, ";", geneEqtlInt$gene)
geneEqtlInt = geneEqtlInt[geneEqtlInt$ID %in% allEqtl$ID,]

exonEqtlInt = meExon$cis$eqtls
exonEqtlInt$gene = as.character(exonEqtlInt$gene)
exonEqtlInt$snps = as.character(exonEqtlInt$snps)
exonEqtlInt$ID = paste0(exonEqtlInt$snps, ";", exonEqtlInt$gene)
exonEqtlInt = exonEqtlInt[exonEqtlInt$ID %in% allEqtl$ID,]

jxnEqtlInt = meJxn$cis$eqtls
jxnEqtlInt$gene = as.character(jxnEqtlInt$gene)
jxnEqtlInt$snps = as.character(jxnEqtlInt$snps)
jxnEqtlInt$ID = paste0(jxnEqtlInt$snps, ";", jxnEqtlInt$gene)
jxnEqtlInt = jxnEqtlInt[jxnEqtlInt$ID %in% allEqtl$ID,]

txEqtlInt = meTx$cis$eqtls
txEqtlInt$gene = as.character(txEqtlInt$gene)
txEqtlInt$snps = as.character(txEqtlInt$snps)
txEqtlInt$ID = paste0(txEqtlInt$snps, ";", txEqtlInt$gene)
txEqtlInt = txEqtlInt[txEqtlInt$ID %in% allEqtl$ID,]

## join
allEqtlInt = rbind(geneEqtlInt, exonEqtlInt, jxnEqtlInt, txEqtlInt)
intStats = allEqtlInt[match(allEqtl$ID, allEqtlInt$ID),]
intStats = intStats[,3:6]
colnames(intStats) = paste0("inter_", colnames(intStats))

allEqtl = cbind(allEqtl, intStats)
allEqtl$ID = NULL

#####################
## get other annotation stats ###
#####################


## snpMap
load("../genotype_data/astellas_dg_genotype_data_n263.rda")
snpMap1 = snpMap
snpMap1$hg19POS = paste0(snpMap1$CHR,":",snpMap1$POS)
snpMap1 = snpMap1[which(rownames(snpMap1) %in% allEqtl$snps),c("SNP","chr_hg38","pos_hg38","hg19POS")]

load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/genotype_data/BrainSeq_Phase2_RiboZero_Genotypes_n551.rda")
snpMap2 = snpMap
snpMap2$hg19POS = paste0(snpMap2$CHR,":",snpMap2$POS)
snpMap2 = snpMap2[which(rownames(snpMap2) %in% allEqtl$snps),c("SNP","chr_hg38","pos_hg38","hg19POS")]

snpMap = snpMap1[snpMap1$hg19POS %in% snpMap2$hg19POS,]

## risk loci from PGC paper + rAggr proxy markers
riskLoci = read.csv("/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_GWAS_riskSNPs/rAggr_results_179.csv", stringsAsFactors=FALSE)	# 10,981 snps
## only keep African and European results based on races present in data
riskLoci =riskLoci[which(riskLoci$Population %in% c("ACB+ASW+ESN+GWD+LWK+MSL+YRI","CEU+FIN+GBR+IBS+TSI")),]
colnames(riskLoci) = gsub("\\.", "_", colnames(riskLoci))
riskLoci$hg19POS = paste0(riskLoci$SNP2_Chr, ":", riskLoci$SNP2_Pos) 

## merge in stats
snpMapMerge = cbind(snpMap, riskLoci[match(snpMap$hg19POS, riskLoci$hg19POS),c(1,8,13,16:18)])
colnames(snpMapMerge) = gsub("SNP1_", "Index_", colnames(snpMapMerge))
colnames(snpMapMerge) = gsub("SNP2_", "Proxy_", colnames(snpMapMerge))

## featMap
library(SummarizedExperiment)
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

### merge in annotation
allEqtlMerge = cbind(allEqtl, snpMapMerge[match(allEqtl$snps, rownames(snpMapMerge)),])
length(unique(allEqtlMerge$snps))# 6277
length(unique(allEqtlMerge$Index_Name))	# 156 tested
save(allEqtlMerge, file = "/dcl01/ajaffe/data/lab/dg_hippo/eQTL_GWAS_riskSNPs_n596/eqtl_tables/mergedEqtl_output_withHippo.rda")

## FDR in any
fdrMat = as.data.frame(allEqtlMerge[,grep("FDR", colnames(allEqtlMerge))]) < 0.05
sigEqtl = allEqtlMerge[which(rowSums(fdrMat, na.rm=TRUE) > 0),]

### save
save(sigEqtl, file = "/dcl01/ajaffe/data/lab/dg_hippo/eQTL_GWAS_riskSNPs_n596/eqtl_tables/mergedEqtl_output_withHippo_fdr05any.rda")

