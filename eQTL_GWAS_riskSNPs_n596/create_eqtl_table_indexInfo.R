##
library(jaffelab)
library(IRanges)
library(SummarizedExperiment)
library(VennDiagram)

################
## load SNP data
## snpMap
load("../genotype_data/astellas_dg_genotype_data_n263.rda")
snpMap1 = snpMap
snpMap1$hg19POS = paste0(snpMap1$CHR,":",snpMap1$POS)
# snpMap1 = snpMap1[which(rownames(snpMap1) %in% c(hippo$snps,dg$snps) ),c("SNP","chr_hg38","pos_hg38","hg19POS")]

load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/genotype_data/BrainSeq_Phase2_RiboZero_Genotypes_n551.rda")
snpMap2 = snpMap
snpMap2$hg19POS = paste0(snpMap2$CHR,":",snpMap2$POS)
# snpMap2 = snpMap2[which(rownames(snpMap2) %in% c(hippo$snps,dg$snps) ),c("SNP","chr_hg38","pos_hg38","hg19POS")]

snpMap = snpMap1[snpMap1$hg19POS %in% snpMap2$hg19POS,]
snpMap$pos_hg19 = snpMap$POS

## risk loci from PGC paper
indexLoci = read.csv("/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_GWAS_riskSNPs/pgc_riskLoci.csv", stringsAsFactors=FALSE) ## 179
indexLoci$hg19POS = paste0(indexLoci$snp_chr,":",indexLoci$snp_pos_hg19)
indexIndex = which(snpMap$hg19POS %in% indexLoci$hg19POS)	# keep 133

## risk loci from PGC paper + rAggr proxy markers
riskLoci = read.csv("/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_GWAS_riskSNPs/rAggr_results_179.csv", stringsAsFactors=FALSE)
riskLoci = riskLoci[which(riskLoci$Population %in% c("ACB+ASW+ESN+GWD+LWK+MSL+YRI","CEU+FIN+GBR+IBS+TSI")),] 
riskLoci_full = riskLoci
colnames(riskLoci) = colnames(riskLoci_full) = gsub("\\.", "_", colnames(riskLoci))
riskLoci$hg19POS1 = paste0(riskLoci$SNP1_Chr, ":", riskLoci$SNP1_Pos) 
riskLoci$hg19POS2 = paste0(riskLoci$SNP2_Chr, ":", riskLoci$SNP2_Pos) 

## keep SNPs from list
keepIndex = which(snpMap$hg19POS %in% riskLoci$hg19POS2)	# keep 10,777 snps from snpMap
snpMap = snpMap[keepIndex,]
keepIndex = which(riskLoci$hg19POS2 %in% snpMap$hg19POS)	# keep 10,754 snps from riskLoci
riskLoci = riskLoci[keepIndex,]

snpMap$Status = ifelse(snpMap$hg19POS %in% indexLoci$hg19POS, "Index","Proxy")
riskLoci$Status1 = ifelse(riskLoci$hg19POS1 %in% indexLoci$hg19POS, "Index","Proxy")
riskLoci$Status2 = ifelse(riskLoci$hg19POS2 %in% indexLoci$hg19POS, "Index","Proxy")

## Also keep track of full list before dropping
riskLoci_full$hg19POS1 = paste0(riskLoci_full$SNP1_Chr, ":", riskLoci_full$SNP1_Pos) 
riskLoci_full$hg19POS2 = paste0(riskLoci_full$SNP2_Chr, ":", riskLoci_full$SNP2_Pos) 
riskLoci_full$SNP2_missing = "missing"
riskLoci_full$SNP2_missing[keepIndex] = "analyzed"
riskLoci_full$Status1 = ifelse(riskLoci_full$hg19POS1 %in% indexLoci$hg19POS, "Index","Proxy")
riskLoci_full$Status2 = ifelse(riskLoci_full$hg19POS2 %in% indexLoci$hg19POS, "Index","Proxy")


### numbers:
# 179 indexSNPs											#  nrow(indexLoci)
# 133 index snps in our SNPs used in eQTLS				#  length(unique(riskLoci$hg19POS2[which(riskLoci$Status2=="Index")]))
# 158 index snps with proxy in our SNPs used in eQTLS	#  length(unique(riskLoci$hg19POS1))
# 7221 rAggr riskLoci (including index SNPs)  			#  length(unique(riskLoci_full$hg19POS2))
# 6330 rAggr riskLoci in our SNPs used in eQTLS 		#  length(unique(riskLoci$hg19POS2)), length(unique(snpMap$hg19POS))



################
## load table
dg = read.csv("raggr_179_snps_dg_eqtls_fdr05.csv", row.names=1)
hippo = read.csv("raggr_179_snps_hippo_eqtls_fdr05.csv", row.names=1)

## unique SNPs
length(unique(dg$hg19POS))
length(unique(hippo$hg19POS))
v = venn.diagram(list(DG = dg$hg19POS, HIPPO = hippo$hg19POS), 
	fill = c("darksalmon", "lightgoldenrod"), main="", main.pos = c(.5, .05), cat.cex = 1.9, cex=3,
	margin = .1, filename = NULL)
pdf("venn_unique_SNP.pdf")
grid.draw(v)
dev.off()

## unique index SNPs
length(unique(dg$hg19POS[dg$hg19POS %in% riskLoci$hg19POS1]))
length(unique(hippo$hg19POS[hippo$hg19POS %in% riskLoci$hg19POS1]))
v = venn.diagram(list(DG = unique(dg$hg19POS[dg$hg19POS %in% riskLoci$hg19POS1]), 
				HIPPO = unique(hippo$hg19POS[hippo$hg19POS %in% riskLoci$hg19POS1])), 
	fill = c("darksalmon", "lightgoldenrod"), main="", main.pos = c(.5, .05), cat.cex = 1.9, cex=3,
	margin = .1, filename = NULL)
pdf("venn_unique_indexsnp.pdf")
grid.draw(v)
dev.off()
	
## unique features
tapply(dg$gene, dg$Type, function(x) length(unique(x)))
tapply(hippo$gene, hippo$Type, function(x) length(unique(x)))
v = venn.diagram(list(DG = dg$gene, HIPPO = hippo$gene), 
	fill = c("darksalmon", "lightgoldenrod"), main="", main.pos = c(.5, .05), cat.cex = 1.9, cex=3,
	margin = .1, filename = NULL)
pdf("venn_unique_features.pdf")
grid.draw(v)
dev.off()

## SNP-feature pairs
nrow(dg)  ## 38609
nrow(hippo)   ## 22572
table(dg$Type)
table(hippo$Type)
v = venn.diagram(list(DG = paste0(dg$SNP, dg$gene), 
				HIPPO = paste0(hippo$SNP, hippo$gene)), 
	fill = c("darksalmon", "lightgoldenrod"), main="", main.pos = c(.5, .05), cat.cex = 1.9, cex=3,
	margin = .1, filename = NULL)
pdf("venn_unique_snpfeatpairs.pdf")
grid.draw(v)
dev.off()


## Unique symbols in SNP-feature pairs
length(unique(dg$Symbol))
length(unique(hippo$Symbol))
tapply(dg$Symbol, dg$Type, function(x) length(unique(x)))
tapply(hippo$Symbol, hippo$Type, function(x) length(unique(x)))
v = venn.diagram(list(DG = dg$Symbol, HIPPO = hippo$Symbol), 
	fill = c("darksalmon", "lightgoldenrod"), main="", main.pos = c(.5, .05), cat.cex = 1.9, cex=3,
	margin = .1, filename = NULL)
pdf("venn_unique_symbol.pdf")
grid.draw(v)
dev.off()
	
	
	
	
################################################################
###### Index SNP info ##########
################################################################

region = hippo

## note which proxy snps have a significant result
riskLoci$proxy_FDRsig = "na"
for (i in 1:nrow(riskLoci)) {
	pos = riskLoci$hg19POS2[i]
	sig = which(region$hg19POS == pos)
	riskLoci$proxy_FDRsig[i] = ifelse(length(sig) > 0, TRUE, FALSE)
}	

## Is SNP index snp or proxy
region$Status = ifelse(region$hg19POS %in% riskLoci$hg19POS1, "Index", "Proxy")

## What is the index snp for each row
proxInd = match(region$hg19POS, riskLoci$hg19POS2)
region$IndexSNP = riskLoci$SNP1_Name[proxInd]
region$IndexSNP_hg19POS = riskLoci$hg19POS1[proxInd]

## was index snp checked in eqtl analysis at all
indexInd = match(region$IndexSNP_hg19POS, riskLoci_full$hg19POS2) ## row of proxy
region$IndexSNP_indata = riskLoci_full$SNP2_missing[indexInd]

## does index snp have any significant eqtl result
region$IndexSNP_fdrSig = "na"
for (i in 1:nrow(region)) {
	if (region$IndexSNP_indata[i] == "analyzed") {
		pos = region$IndexSNP_hg19POS[i]
		sig = which(region$hg19POS == pos)
		region$IndexSNP_fdrSig[i] = ifelse(length(sig) > 0, TRUE, FALSE)
	}
}
## what is the most significant eqtl result for the index snp
region$IndexSNP_mostSigFeat = NA
region$IndexSNP_mostSigFeat_gene = NA
for (i in 1:nrow(region)) {
	if (region$IndexSNP_fdrSig[i] == "TRUE") {
		pos = region$IndexSNP_hg19POS[i]
		tmp = region[which(region$hg19POS == pos),]
		tmp = tmp[order(tmp$pvalue, decreasing=FALSE),]
		region$IndexSNP_mostSigFeat[i] = as.character(tmp$gene[1])
		region$IndexSNP_mostSigFeat_gene[i] = as.character(tmp$Symbol[1])
	}
}

## what is the lead variant for each SNP -- index if sig, else highest LD proxy
region$leadVariant = NA
for (i in 1:nrow(region)) {
	if (region$IndexSNP_fdrSig[i] == "TRUE") {  ## index snp is fdr significant
		region$leadVariant[i] = region$IndexSNP[i]
		
	} else {									## index snp is NOT fdr significant
		## find highest LD proxy snp that is sig
		pos = region$IndexSNP_hg19POS[i]
		tmp = riskLoci[which(riskLoci$hg19POS1 == pos),]
		t_ind = which(tmp$Distance==0 | tmp$proxy_FDRsig==FALSE)
		if (length(t_ind>0)) { tmp = tmp[-t_ind,] }
		tmp = tmp[order(tmp$R_squared, rev(tmp$Distance), decreasing=TRUE),]
		region$leadVariant[i] = as.character(tmp$SNP2_Name[1])
		}
}
## Is the SNP the lead variant? (Check by position)
leadVarInd = match(region$leadVariant, riskLoci$SNP2_Name)
leadVarPos = riskLoci$hg19POS2[leadVarInd]
region$leadVariant_indicator = (region$hg19POS == leadVarPos)


hippo = region


write.csv(dg, file="raggr_179_snps_dg_eqtls_fdr05.csv")
write.csv(hippo, file="raggr_179_snps_hippo_eqtls_fdr05.csv")


dg$pair = paste0(dg$gene, "_snp:", dg$hg19POS)
hippo$pair = paste0(hippo$gene, "_snp:", hippo$hg19POS)
### Numbers
## DG
# unique SNPs:  4181					# length(unique(dg$hg19POS))
# index SNPs with FDR < 0.05:  75		# length(unique(dg$hg19POS[which(dg$Status=="Index")]))
# index SNPs with proxy < 0.05:  116	# length(unique(dg$IndexSNP_hg19POS))
# unique gene symbols:	432				# length(unique(dg$Symbol))-1   # minus 1 to take out NA/""

# unique features:						# tapply(dg$gene, dg$Type, function(x) length(unique(x)))
# Exon Gene  Jxn   Tx
# 1144  172  365  318

# from # gene symbols:					# tapply(dg$Symbol, dg$Type, function(x) length(unique(x)))
# Exon Gene  Jxn   Tx
# 256  146  173  233

# unique snp-feature pairs:				# tapply(dg$pair, dg$Type, function(x) length(unique(x)))
# Exon   Gene  Jxn   Tx
# 39894  6426  11542 10291


## HIPPO
# unique SNPs:  4006					# length(unique(hippo$hg19POS))
# index SNPs with FDR < 0.05:  73		# length(unique(hippo$hg19POS[which(hippo$Status=="Index")]))
# index SNPs with proxy < 0.05:  105	# length(unique(hippo$IndexSNP_hg19POS))
# unique gene symbols:	436				# length(unique(hippo$Symbol))-1   # minus 1 to take out NA/""

# unique features:						# tapply(hippo$gene, hippo$Type, function(x) length(unique(x)))
# Exon Gene  Jxn   Tx
# 1052  135  357  303

# from # gene symbols:					# tapply(hippo$Symbol, hippo$Type, function(x) length(unique(x)))
# Exon Gene  Jxn   Tx
# 259  115  168  221

# unique snp-feature pairs:				# tapply(hippo$pair, hippo$Type, function(x) length(unique(x)))
# Exon   Gene  Jxn   Tx
# 31267  4186  9827  9444







library(VennDiagram)

## unique index SNPs
v = venn.diagram(list(DG = unique(dg$IndexSNP_hg19POS), 
				HIPPO = unique(hippo$IndexSNP_hg19POS)), 
	fill = c("darksalmon", "lightgoldenrod"), main="", main.pos = c(.5, .05), cat.cex = 1.9, cex=3,
	margin = .1, filename = NULL)
pdf("venn_unique_indexsnp_withproxy.pdf")
grid.draw(v)
dev.off()
	
## unique gene symbols
v = venn.diagram(list(DG = unique(dg$Symbol[!is.na(dg$Symbol)]), 
				HIPPO = unique(hippo$Symbol[!is.na(hippo$Symbol)])), 
	fill = c("darksalmon", "lightgoldenrod"), main="", main.pos = c(.5, .05), cat.cex = 1.9, cex=3,
	margin = .1, filename = NULL)
pdf("venn_unique_gene_symbols.pdf")
grid.draw(v)
dev.off()
	









