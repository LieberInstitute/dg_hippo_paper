##
####
library(SummarizedExperiment)
library(jaffelab)
library(MatrixEQTL)
library(sva)
library(recount)

## overlapping snps
load("../genotype_data/BrainSeq_Phase2_RiboZero_Genotypes_n551.rda")
snpMap$pos_hg19 = paste0(snpMap$CHR, ":", snpMap$POS)
snpMap_hippo = snpMap

######################
### load data ####
######################

## load DG
load("../count_data/merged_dg_hippo_allSamples_n596.rda", verbose=TRUE)

## keep adult samples & correct region
keepInd = which(rse_gene_joint$Region == "DentateGyrus")	# 263
rse_gene = rse_gene_joint[,keepInd]
rse_exon = rse_exon_joint[,keepInd]
rse_jxn = rse_jxn_joint[,keepInd]
rse_tx = rse_tx_joint[,keepInd]

pd = colData(rse_gene)

## load SNP data
load("../genotype_data/astellas_dg_genotype_data_n263.rda")
snpMap$pos_hg19 = paste0(snpMap$CHR, ":", snpMap$POS)

## drop rs10708380:150158001:TG:T (missing info in snpMap (and dbSNP))
snpInd = which(rownames(snpMap) == "rs10708380:150158001:TG:T")
snpMap = snpMap[-snpInd,]
snp = snp[-snpInd,]

## drop those with missing hg38 coords
keepSnps = which(!is.na(snpMap$pos_hg38))
snpMap = snpMap[keepSnps,]
snp = snp[keepSnps,]

#####################
# filter brain region
# make mds and snp dimensions equal to N
# (repeat rows or columns for BrNum replicates)
mds = mds[pd$BrNum,]
snp = snp[,pd$BrNum]
rownames(mds) = colnames(snp) = pd$RNum

######################
# statistical model ##
######################
pd$Dx = factor(pd$Dx,
	levels = c("Control", "Schizo", "Bipolar", "MDD"))

mod = model.matrix(~Dx + Sex + as.matrix(mds[,1:5]), data = pd)
colnames(mod)[4:8] = colnames(mds)[1:5]

##########################
### read in HIPPO eQTLs ##
##########################

## from brainseq phase2
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_full/eqtl_tables/mergedEqtl_output_hippo_4features.rda",verbose=TRUE)
sigEqtl = allEqtl[allEqtl$FDR < 0.01,] # just take those significant

######################
# create SNP objects #

table(snpMap_hippo$pos_hg19 %in% snpMap$pos_hg19)
table(snpMap$pos_hg19 %in% snpMap_hippo$pos_hg19)


######################
table(rownames(snp) %in% sigEqtl$snps)
table( sigEqtl$snps %in% rownames(snp))
snp2 = snp[rownames(snp) %in% sigEqtl$snps,]

theSnps = SlicedData$new(as.matrix(snp2))
theSnps$ResliceCombined(sliceSize = 50000)

snpspos = snpMap[rownames(snpMap) %in% sigEqtl$snps,c("SNP","chr_hg38","pos_hg38")]
colnames(snpspos) = c("name","chr","pos")


#######################
####### do PCA ########
#######################

geneRpkm = recount::getRPKM(rse_gene, "Length")
exonRpkm = recount::getRPKM(rse_exon, "Length")
rowRanges(rse_jxn)$Length <- 100
jxnRp10m = recount::getRPKM(rse_jxn, "Length")
txTpm = assays(rse_tx)$tpm

## load PCs 
load("rdas/pcs_4features_dg_263.rda")

## make models
covsGene = SlicedData$new(t(cbind(mod[,-1],genePCs)))
covsExon = SlicedData$new(t(cbind(mod[,-1],exonPCs)))
covsJxn = SlicedData$new(t(cbind(mod[,-1],jxnPCs)))
covsTx = SlicedData$new(t(cbind(mod[,-1],txPCs)))

##########################
### feature annotation ###
##########################

###### gene level
posGene = as.data.frame(rowRanges(rse_gene))[,1:3]
posGene$name = rownames(posGene)
posGene = posGene[,c(4,1:3)]

##### exon level 
posExon = as.data.frame(rowRanges(rse_exon))[,1:3]
posExon$name = rownames(posExon)
posExon = posExon[,c(4,1:3)]

##### junction level 
posJxn = as.data.frame(rowRanges(rse_jxn))[,1:3]
posJxn$name = rownames(posJxn)
posJxn = posJxn[,c(4,1:3)]
names(posJxn)[2:4] = c("Chr", "Start","End")

##### transcript level 
posTx = as.data.frame(rowRanges(rse_tx))[,1:3]
posTx$name = rownames(posTx)
posTx = posTx[,c(4,1:3)]
names(posTx)[2:4] = c("Chr", "Start","End")

#############################
### sliced expression data ##

## filter to eQTLs #####
geneRpkm = geneRpkm[rownames(geneRpkm) %in% sigEqtl$gene,]
exonRpkm = exonRpkm[rownames(exonRpkm) %in% sigEqtl$gene,]
jxnRp10m = jxnRp10m[rownames(jxnRp10m) %in% sigEqtl$gene,]
txTpm = txTpm[rownames(txTpm) %in% sigEqtl$gene,]

## slice up
geneSlice = SlicedData$new(log2(geneRpkm+1))
exonSlice = SlicedData$new(log2(exonRpkm+1))
jxnSlice = SlicedData$new(log2(jxnRp10m+1))
txSlice = SlicedData$new(log2(txTpm+1))

geneSlice$ResliceCombined(sliceSize = 5000)
exonSlice$ResliceCombined(sliceSize = 5000)
jxnSlice$ResliceCombined(sliceSize = 5000)
txSlice$ResliceCombined(sliceSize = 5000)


##########################
### Run EQTLs ############
##########################

print("Starting eQTLs")

meGene = Matrix_eQTL_main(snps=theSnps, gene = geneSlice, 
	cvrt = covsGene, output_file_name.cis =  ".ctxt" ,
	pvOutputThreshold.cis = 1,  pvOutputThreshold=0,
	snpspos = snpspos, genepos = posGene, 
	useModel = modelLINEAR,	cisDist=5e5,
	pvalue.hist = 100,min.pv.by.genesnp = TRUE)	
	
meExon = Matrix_eQTL_main(snps=theSnps, gene = exonSlice, 
	cvrt = covsExon, output_file_name.cis =  ".ctxt" ,
	pvOutputThreshold.cis = 1,  pvOutputThreshold=0,
	snpspos = snpspos, genepos = posExon, 
	useModel = modelLINEAR,	cisDist=5e5,
	pvalue.hist = 100,min.pv.by.genesnp = TRUE)	

meJxn = Matrix_eQTL_main(snps=theSnps, gene = jxnSlice, 
	cvrt = covsJxn, output_file_name.cis =  ".ctxt" ,
	pvOutputThreshold.cis = 1,  pvOutputThreshold=0,
	snpspos = snpspos, genepos = posJxn, 
	useModel = modelLINEAR,	cisDist=5e5,
	pvalue.hist = 100,min.pv.by.genesnp = TRUE)	
	
meTx = Matrix_eQTL_main(snps=theSnps, gene = txSlice, 
	cvrt = covsTx, output_file_name.cis =  ".ctxt" ,
	pvOutputThreshold.cis = 1,  pvOutputThreshold=0,
	snpspos = snpspos, genepos = posTx, 
	useModel = modelLINEAR,	cisDist=5e5,
	pvalue.hist = 100,min.pv.by.genesnp = TRUE)	

save(meGene, meExon, meJxn, meTx,
	file="eqtl_tables/matrixEqtl_output_hippo_4features_dgMatched.rda")

	
######################
###### annotate ######

sigEqtl$snp_feature = paste0(sigEqtl$snps, ";", sigEqtl$gene) # hippo

# extract
geneEqtl = meGene$cis$eqtls
geneEqtl$gene = as.character(geneEqtl$gene)
geneEqtl$snps = as.character(geneEqtl$snps)
geneEqtl$snp_feature = paste0(geneEqtl$snps, ";", geneEqtl$gene)
geneEqtl = geneEqtl[which(geneEqtl$snp_feature %in% sigEqtl$snp_feature),]

exonEqtl = meExon$cis$eqtls
exonEqtl$gene = as.character(exonEqtl$gene)
exonEqtl$snps = as.character(exonEqtl$snps)
exonEqtl$snp_feature = paste0(exonEqtl$snps, ";", exonEqtl$gene)
exonEqtl = exonEqtl[which(exonEqtl$snp_feature %in% sigEqtl$snp_feature),]

jxnEqtl = meJxn$cis$eqtls
jxnEqtl$gene = as.character(jxnEqtl$gene)
jxnEqtl$snps = as.character(jxnEqtl$snps)
jxnEqtl$snp_feature = paste0(jxnEqtl$snps, ";", jxnEqtl$gene)
jxnEqtl = jxnEqtl[which(jxnEqtl$snp_feature %in% sigEqtl$snp_feature),]

txEqtl = meTx$cis$eqtls
txEqtl$gene = as.character(txEqtl$gene)
txEqtl$snps = as.character(txEqtl$snps)
txEqtl$snp_feature = paste0(txEqtl$snps, ";", txEqtl$gene)
txEqtl = txEqtl[which(txEqtl$snp_feature %in% sigEqtl$snp_feature),]

allEqtl = rbind(geneEqtl, exonEqtl, jxnEqtl, txEqtl)

## match up
mm = match(sigEqtl$snp_feature, allEqtl$snp_feature)
sigEqtl$dg_statistic = allEqtl$statistic[mm]
sigEqtl$dg_pvalue = allEqtl$pvalue[mm]
sigEqtl$dg_beta = allEqtl$beta[mm]
sigEqtl$dg_FDR = allEqtl$FDR[mm]

## drop anno columns
sigEqtl$gencodeTx = sigEqtl$snp_feature = NULL
save(sigEqtl, file="eqtl_tables/mergedEqtl_output_hippo_4features_withDg.rda",compress=TRUE)

################
## metrics #####
library(qvalue)

nrow(sigEqtl)
mean(!is.na(sigEqtl$dg_statistic))
table(is.na(sigEqtl$dg_statistic))

length(unique(sigEqtl$gene[!sigEqtl$gene %in% allEqtl$gene]))
length(unique(sigEqtl$snps[!sigEqtl$snps %in% allEqtl$snps]))

length(unique(sigEqtl$gene[is.na(sigEqtl$dg_statistic)]))
length(unique(sigEqtl$snps[is.na(sigEqtl$dg_statistic)]))


mean(sign(sigEqtl$statistic) == sign(sigEqtl$dg_statistic),na.rm=TRUE)
mean(sign(sigEqtl$statistic) == sign(sigEqtl$dg_statistic) & 
	sigEqtl$dg_pvalue < 0.05,na.rm=TRUE)
mean(sign(sigEqtl$statistic) == sign(sigEqtl$dg_statistic) & 
	sigEqtl$dg_FDR < 0.01,na.rm=TRUE)
qv_dg = qvalue(sigEqtl$dg_pvalue)
1 - qv_dg$pi0 # 84.0


