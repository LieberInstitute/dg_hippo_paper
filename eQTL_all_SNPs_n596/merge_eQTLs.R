###

library(GenomicRanges)

## get DG eQTLs
load("eqtl_tables/mergedEqtl_output_dg_4features_fdr01.rda")
sigEqtl$ID = paste0(sigEqtl$snps, ";", sigEqtl$gene)

########################
## get hippo eQTLs
load("eqtl_tables/matrixEqtl_output_hippo_4features_dgMatched.rda")

# extract
geneEqtlH = meGene$cis$eqtls
geneEqtlH$gene = as.character(geneEqtlH$gene)
geneEqtlH$snps = as.character(geneEqtlH$snps)
geneEqtlH$ID = paste0(geneEqtlH$snps, ";", geneEqtlH$gene)
geneEqtlH = geneEqtlH[geneEqtlH$ID %in% sigEqtl$ID,]

exonEqtlH = meExon$cis$eqtls
exonEqtlH$gene = as.character(exonEqtlH$gene)
exonEqtlH$snps = as.character(exonEqtlH$snps)
exonEqtlH$ID = paste0(exonEqtlH$snps, ";", exonEqtlH$gene)
exonEqtlH = exonEqtlH[exonEqtlH$ID %in% sigEqtl$ID,]

jxnEqtlH = meJxn$cis$eqtls
jxnEqtlH$gene = as.character(jxnEqtlH$gene)
jxnEqtlH$snps = as.character(jxnEqtlH$snps)
jxnEqtlH$ID = paste0(jxnEqtlH$snps, ";", jxnEqtlH$gene)
jxnEqtlH = jxnEqtlH[jxnEqtlH$ID %in% sigEqtl$ID,]

txEqtlH = meTx$cis$eqtls
txEqtlH$gene = as.character(txEqtlH$gene)
txEqtlH$snps = as.character(txEqtlH$snps)
txEqtlH$ID = paste0(txEqtlH$snps, ";", txEqtlH$gene)
txEqtlH = txEqtlH[txEqtlH$ID %in% sigEqtl$ID,]

## join
allEqtlH = rbind(geneEqtlH, exonEqtlH, jxnEqtlH, txEqtlH)
hippoStats = allEqtlH[match(sigEqtl$ID, allEqtlH$ID),]
hippoStats = hippoStats[,3:6]
colnames(hippoStats) = paste0("hippo_", colnames(hippoStats))

sigEqtl = cbind(sigEqtl, hippoStats)

##################################
### interaction ##################
##################################

load("eqtl_tables/matrixEqtl_output_interaction_4features_dgMatched.rda")

# extract
geneEqtlInt = meGene$cis$eqtls
geneEqtlInt$gene = as.character(geneEqtlInt$gene)
geneEqtlInt$snps = as.character(geneEqtlInt$snps)
geneEqtlInt$ID = paste0(geneEqtlInt$snps, ";", geneEqtlInt$gene)
geneEqtlInt = geneEqtlInt[geneEqtlInt$ID %in% sigEqtl$ID,]

exonEqtlInt = meExon$cis$eqtls
exonEqtlInt$gene = as.character(exonEqtlInt$gene)
exonEqtlInt$snps = as.character(exonEqtlInt$snps)
exonEqtlInt$ID = paste0(exonEqtlInt$snps, ";", exonEqtlInt$gene)
exonEqtlInt = exonEqtlInt[exonEqtlInt$ID %in% sigEqtl$ID,]

jxnEqtlInt = meJxn$cis$eqtls
jxnEqtlInt$gene = as.character(jxnEqtlInt$gene)
jxnEqtlInt$snps = as.character(jxnEqtlInt$snps)
jxnEqtlInt$ID = paste0(jxnEqtlInt$snps, ";", jxnEqtlInt$gene)
jxnEqtlInt = jxnEqtlInt[jxnEqtlInt$ID %in% sigEqtl$ID,]

txEqtlInt = meTx$cis$eqtls
txEqtlInt$gene = as.character(txEqtlInt$gene)
txEqtlInt$snps = as.character(txEqtlInt$snps)
txEqtlInt$ID = paste0(txEqtlInt$snps, ";", txEqtlInt$gene)
txEqtlInt = txEqtlInt[txEqtlInt$ID %in% sigEqtl$ID,]

## join
allEqtlInt = rbind(geneEqtlInt, exonEqtlInt, jxnEqtlInt, txEqtlInt)
intStats = allEqtlInt[match(sigEqtl$ID, allEqtlInt$ID),]
intStats = intStats[,3:6]
colnames(intStats) = paste0("inter_", colnames(intStats))

sigEqtl = cbind(sigEqtl, intStats)
sigEqtl$ID = NULL

### save
save(sigEqtl, file = "eqtl_tables/mergedEqtl_output_dg_4features_fdr01_withHippo.rda")

#####################
# export text file ##
#####################

## for feature coordinates
load("../rdas/exon_names_hg19_hg38.rda", verbose=TRUE)

## load SNP info
load("../genotype_data/BrainSeq_Phase2_RiboZero_Genotypes_n551.rda")
snpMap$pos_hg19 = paste0(snpMap$CHR, ":", snpMap$POS)
snpMap_hippo = snpMap

load("../genotype_data/astellas_dg_genotype_data_n263.rda")
snpMap_dg = snpMap
snpMap_dg$pos_hg19 = paste0("chr", snpMap_dg$CHR, ":", snpMap_dg$POS)
snpMap_dg$pos_hg38 = paste0(snpMap_dg$chr_hg38, ":", snpMap_dg$pos_hg38)

rm(snpMap, snp,mds)

sigEqtlDf = as.data.frame(sigEqtl)
sigEqtlDf$gencodeTx = NULL
sigEqtlDf$snpRsNum = snpMap_dg$name[match(sigEqtlDf$snps, snpMap_dg$SNP)]
sigEqtlDf$snpCoord_hg19 = snpMap_dg$pos_hg19[match(sigEqtlDf$snps, snpMap_dg$SNP)]
sigEqtlDf$snpCoord_hg38 = snpMap_dg$pos_hg38[match(sigEqtlDf$snps, snpMap_dg$SNP)]
sigEqtlDf$alleleCounted = snpMap_dg$newCount[match(sigEqtlDf$snps, snpMap_dg$SNP)]
sigEqtlDf$alleleRef = snpMap_dg$newRef[match(sigEqtlDf$snps, snpMap_dg$SNP)]

## clean up
colnames(sigEqtlDf)[1:2] = c("SNP", "Feature")
sigEqtlDf$exonInternalID = sigEqtlDf$Feature
sigEqtlDf$exonInternalID[sigEqtlDf$Type != "Exon"] = NA

sigEqtlDf$Feature[sigEqtlDf$Type == "Exon"] = hg38_exons$gencode_exonID[
	match(sigEqtlDf$Feature[sigEqtlDf$Type == "Exon"], hg38_exons$hg38_eID)]
sigEqtlDf = sigEqtlDf[order(sigEqtlDf$pvalue),]
rownames(sigEqtlDf) = NULL

write.csv(sigEqtlDf, row.names=FALSE,
	file = gzfile("../genotype_data/DataS1_DGGCL_eQTLs_plusHippo.csv.gz"))
