###
####
library(GenomicRanges)
library(VennDiagram)
library(jaffelab)
library(SummarizedExperiment)
library(RColorBrewer)
library(readxl)
library(ggplot2)

## load eQTLs
load("eQTL_all_SNPs_n596/eqtl_tables/mergedEqtl_output_dg_4features_fdr01_withHippo.rda")
sigEqtl$Type = factor(sigEqtl$Type, levels = c("Gene", "Exon",  "Jxn", "Tx"))
sigEqtl$EnsemblGeneID = ss(sigEqtl$EnsemblGeneID, "\\.")

## load data
load("rdas/cleaned_exprs_data_twoData_eQTL.rda")

## which eQTLs are independent vs dependent
sigEqtlFeatureList =split(sigEqtl, sigEqtl$gene)
## maybe correlation of stats by feature
corFeature = sapply(sigEqtlFeatureList[sapply(sigEqtlFeatureList,nrow) > 1],
	function(x) cor(x$statistic, x$hippo_statistic,use="pair"))
	
d = data.frame(FeatureID = names(corFeature), cor = corFeature,
	stringsAsFactors = FALSE)
d$Type = sigEqtl$Type[match(d$FeatureID, sigEqtl$gene)]
d$EnsemblGeneID = sigEqtl$EnsemblGeneID[match(d$FeatureID, sigEqtl$gene)]

g = ggplot(d, aes(x=Type, y=cor)) + geom_violin()
ggsave(g, file="plots/eqtl_corr.pdf")