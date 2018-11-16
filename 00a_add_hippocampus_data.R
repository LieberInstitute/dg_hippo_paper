##

### libraries
library(SummarizedExperiment)
library(jaffelab)
library(readxl)

## load DG
load("count_data/astellas_dg_hg38_rseTx_n263.rda")
load("count_data/astellas_dg_hg38_rseJxn_n263.rda")
load("count_data/astellas_dg_hg38_rseExon_n263.rda")
load("count_data/astellas_dg_hg38_rseGene_n263.rda")

## rename objects
rse_gene_dg = rse_gene
rse_exon_dg = rse_exon
rse_jxn_dg = rse_jxn
rse_tx_dg = rse_tx
pdDg = colData(rse_gene_dg)
pdDg$Dataset = "DG"
for(i in grep("integer", sapply(pdDg,class))) pdDg[,i] = as.numeric(pdDg[,i])
pdDg$trimmed = as.character(pdDg$trimmed) # to match

## load hippocampus bulk
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/count_data/hippo_brainseq_phase2_hg38_rseTx_merged_n447.rda")
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/count_data/hippo_brainseq_phase2_hg38_rseJxn_merged_n447.rda")
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/count_data/hippo_brainseq_phase2_hg38_rseExon_merged_n447.rda")
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/count_data/hippo_brainseq_phase2_hg38_rseGene_merged_n447.rda")

## add kit info
hipxl <- read_excel('/dcl01/lieber/ajaffe/lab/brainseq_phase2/misc/LIBD_PhaseII_HIPPO_RiboZero_sample_list_01_28_2015.xlsx')
rse_gene$Protocol = hipxl$Protocol[match(rse_gene$RNum, paste0("R", hipxl$RNum))]

## filter to adults, age 17, no gold
hippoIndex = which(rse_gene$Protocol == "RiboZeroHMR")
rse_tx = rse_tx[,hippoIndex]
rse_gene = rse_gene[,hippoIndex]
rse_exon = rse_exon[,hippoIndex]
rse_jxn = rse_jxn[,hippoIndex]

pdH = colData(rse_gene)
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

colData(rse_gene) = pdH[,n]
colData(rse_exon) = pdH[,n]
colData(rse_jxn) = pdH[,n]
colData(rse_tx) = pdH[,n]

## check classes
sapply(colData(rse_gene), class)[
	sapply(colData(rse_gene_dg), class) != sapply(colData(rse_gene), class)]
# looks good

### match up
mm = match(pdDg$BrNum, pdH$BrNum)
table(!is.na(mm)) # 112

mcols(rse_gene)$meanExprs = mcols(rse_gene_dg)$meanExprs = NULL
rse_gene_joint = cbind(rse_gene_dg[,!is.na(mm)], rse_gene[,mm[!is.na(mm)]])

mcols(rse_exon)$meanExprs = mcols(rse_exon_dg)$meanExprs = NULL
rse_exon_joint = cbind(rse_exon_dg[,!is.na(mm)], rse_exon[,mm[!is.na(mm)]])

# merge junction data
mcols(rse_jxn)$meanExprs = mcols(rse_jxn_dg)$meanExprs = NULL
j = intersect(rownames(rse_jxn), rownames(rse_jxn_dg))
rse_jxn_joint = cbind(rse_jxn_dg[j,!is.na(mm)], rse_jxn[j,mm[!is.na(mm)]])

# merge tx data
rse_tx_joint = cbind(rse_tx_dg, rse_tx)

## save
save(rse_gene_joint, file = paste0("count_data/dgPlusHippo_hg38_rseGene_n",
		ncol(rse_gene_joint),".rda"))
save(rse_exon_joint, file = paste0("count_data/dgPlusHippo_hg38_rseExon_n",
		ncol(rse_exon_joint),".rda"))
save(rse_jxn_joint, file = paste0("count_data/dgPlusHippo_hg38_rseJxn_n",
		ncol(rse_jxn_joint),".rda"))
save(rse_tx_joint, file =  paste0("count_data/dgPlusHippo_hg38_rseTx_n",
		ncol(rse_tx_joint),".rda"))