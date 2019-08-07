##
library(IRanges)

## marginal eQTLs
load("eqtl_tables/mergedEqtl_output_dg_4features_fdr01_withHippo.rda")
geneEqtl = sigEqtl[sigEqtl$Type == "Gene",]
geneEqtl$ID = paste0(geneEqtl$gene, ";", geneEqtl$snps)

## from BootstrapQTL
# load("rdas/DG_genes_Huang-bootstrap_results_BH-BH.Rdata")
# eGenes = eGenes.DG.BHBH 
load("rdas/DG_genes_Huang-bootstrap_results_Bonf-BH.Rdata")
eGenes = eGenes.DG
eGenes$ID = paste0(eGenes$eGene, ";", eGenes$eSNPs)

## add FDR to bootstrap output
mm1 = match(eGenes$ID, geneEqtl$ID)
eGenes$FDR = geneEqtl$FDR[mm1]

length(unique(eGenes$eGene)) # 6342

## how many were missing?
length(unique(eGenes$eGene[is.na(eGenes$FDR)])) # 415
length(unique(eGenes$eGene[is.na(eGenes$FDR)]))/
	length(unique(eGenes$eGene)) # 7.4%
length(unique(eGenes$eSNPs[is.na(eGenes$FDR)])) # 20191
length(unique(eGenes$eSNPs[is.na(eGenes$FDR)]))/
	length(unique(eGenes$eSNPs)) # 5.4%


## other checks
hist(eGenes$winners_curse)
hist(eGenes$corrected_beta/eGenes$nominal_beta)
quantile(eGenes$corrected_beta/eGenes$nominal_beta)

## add bootstrap output to FDR
mm2 =  match(geneEqtl$ID, eGenes$ID)
geneEqtl$corrected_beta = eGenes$corrected_beta[mm2]
geneEqtl$winners_curse = eGenes$winners_curse[mm2]
geneEqtl$eSNP_pval = eGenes$eSNP_pval[mm2]
geneEqtl$eGene_pval = eGenes$eGene_pval[mm2]

g = unique(geneEqtl$gene)
table(g %in% geneEqtl$gene[!is.na(geneEqtl$winners_curse)])
mean(g %in% geneEqtl$gene[!is.na(geneEqtl$winners_curse)])

table(is.na(geneEqtl$winners_curse), geneEqtl$hippo_pvalue < 0.05,
	dnn = c("Boot", "HippoRep"))
table(is.na(geneEqtl$winners_curse), geneEqtl$hippo_pvalue < 0.01,
	dnn = c("Boot", "HippoRep"))
table(is.na(geneEqtl$winners_curse), geneEqtl$hippo_FDR < 0.05,
	dnn = c("Boot", "HippoRep"))