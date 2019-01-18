##### libraries
library(SummarizedExperiment)
library(jaffelab)
library(readxl)

###############
## DG #########
###############

## load DG gene + Tx
load("count_data/astellas_dg_hg38_rseGene_n263.rda")
rse_gene_dg = rse_gene
load("count_data/astellas_dg_hg38_rseTx_n263.rda")
rse_tx_dg = rse_tx
load("count_data/astellas_dg_hg38_rseJxn_n263.rda")
rse_jxn_dg = rse_jxn
load("count_data/astellas_dg_hg38_rseExon_n263.rda")
rse_exon_dg = rse_exon

pdDg = colData(rse_gene_dg)
pdDg$Dataset = "DG"
for(i in grep("integer", sapply(pdDg,class))) pdDg[,i] = as.numeric(pdDg[,i])
pdDg$trimmed = as.character(pdDg$trimmed) # to match


######################
### HIPPO ############
######################

load("count_data/hippo_brainseq_phase2_hg38_rseTx_merged_n447.rda")
rse_tx_hippo = rse_tx
load("count_data/hippo_brainseq_phase2_hg38_rseGene_merged_n447.rda")
rse_gene_hippo = rse_gene
load("count_data/hippo_brainseq_phase2_hg38_rseExon_merged_n447.rda")
rse_exon_hippo = rse_exon
load("count_data/hippo_brainseq_phase2_hg38_rseJxn_merged_n447.rda")
rse_jxn_hippo = rse_jxn

## add kit info
hipxl <- read_excel('count_data/LIBD_PhaseII_HIPPO_RiboZero_sample_list_01_28_2015.xlsx')
rse_gene_hippo$Protocol = hipxl$Protocol[match(rse_gene_hippo$RNum, paste0("R", hipxl$RNum))]

## filter to adults, age 17, no gold
hippoIndex = which(colData(rse_gene_hippo)$Age > 17 & rse_gene_hippo$Protocol == "RiboZeroHMR")
rse_tx_hippo = rse_tx_hippo[,hippoIndex]
rse_gene_hippo = rse_gene_hippo[,hippoIndex]
rse_exon_hippo = rse_exon_hippo[,hippoIndex]
rse_jxn_hippo = rse_jxn_hippo[,hippoIndex]


pdH = colData(rse_gene_hippo)
pdH$Dataset = "Hippo"

## fix numeric list stuff
classH = sapply(pdH,class)
for(i in grep("CharacterList", classH)) pdH[,i] = sapply(pdH[,i], paste, collapse=";")
for(i in grep("LogicalList", classH)) pdH[,i] = sapply(pdH[,i], paste, collapse=";")
for(i in grep("NumericList", classH)) pdH[,i] = sapply(pdH[,i], mean)
for(i in grep("IntegerList", classH)) pdH[,i] = as.numeric(sapply(pdH[,i], mean))

## make phenotype columns merge-able 
n = intersect(colnames(pdDg), colnames(pdH))
colData(rse_gene_dg) = pdDg[,n]
colData(rse_exon_dg) = pdDg[,n]
colData(rse_jxn_dg) = pdDg[,n]
colData(rse_tx_dg) = pdDg[,n]

colData(rse_gene_hippo) = pdH[,n]
colData(rse_exon_hippo) = pdH[,n]
colData(rse_jxn_hippo) = pdH[,n]
colData(rse_tx_hippo) = pdH[,n]

## merge
mcols(rse_gene_hippo)$meanExprs = mcols(rse_gene_dg)$meanExprs = NULL
rse_gene_joint = cbind(rse_gene_dg, rse_gene_hippo)

mcols(rse_exon_hippo)$meanExprs = mcols(rse_exon_dg)$meanExprs = NULL
rse_exon_joint = cbind(rse_exon_dg, rse_exon_hippo)

# merge junction data
mcols(rse_jxn_hippo)$meanExprs = mcols(rse_jxn_dg)$meanExprs = NULL
j = intersect(rownames(rse_jxn_hippo), rownames(rse_jxn_dg))
rse_jxn_joint = cbind(rse_jxn_dg[j,], rse_jxn_hippo[j,])

# merge tx data
rse_tx_joint = cbind(rse_tx_dg, rse_tx_hippo)

###############################
######### filter ##############
###############################
library(recount)
geneIndex = rowMeans(recount::getRPKM(rse_gene_joint[,rse_gene_joint$Dataset == "DG"], "Length")) > 0.5
exonIndex = rowMeans(recount::getRPKM(rse_exon_joint[,rse_exon_joint$Dataset == "DG"], "Length")) > 0.5
rowData(rse_jxn_joint)$Length=100
jxnIndex = rowMeans(recount::getRPKM(rse_jxn_joint[,rse_jxn_joint$Dataset == "DG"], "Length")) > 0.5 &
					rowData(rse_jxn_joint)$Class != "Novel"
txIndex = rowMeans(assays(rse_tx_joint)$tpm) > 0.5 

rse_gene_joint = rse_gene_joint[geneIndex,]
rse_exon_joint = rse_exon_joint[exonIndex,]
rse_jxn_joint = rse_jxn_joint[jxnIndex,]
rse_tx_joint = rse_tx_joint[txIndex,]

###############################
## add qSVs ###################
###############################

load("count_data/degradation_rse_dg_hippo_n263.rda")
rse_deg_dg = cov_rse

load("count_data/degradation_rse_phase2_hippo.rda")
rse_deg_hippo = cov_rse_hippo
colnames(rse_deg_hippo) = ss(colnames(rse_deg_hippo), "_", 1)
identical(granges(rse_deg_dg), granges(rse_deg_hippo)) # TRUE

rse_deg_joint = cbind(rse_deg_dg, rse_deg_hippo)
rse_deg_joint = rse_deg_joint[,ss(colnames(rse_gene_joint), "_B")]

n = ncol(rse_gene_joint)
save(rse_gene_joint, rse_exon_joint, rse_jxn_joint, rse_tx_joint, rse_deg_joint,
	file = paste0("count_data/merged_dg_hippo_allSamples_n", n, ".rda"))
	