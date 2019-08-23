##
library(IRanges)
library(jaffelab)

## marginal eQTLs
load("eqtl_tables/mergedEqtl_output_dg_4features_fdr01_withHippo.rda")

load("../rdas/num_snps_per_gene_eqtl.rda")

###########################	
## from BootstrapQTL
load("rdas/DG_genes_Huang-bootstrap_results_BH-BH.Rdata",verbose=TRUE)
geneEqtl = sigEqtl[sigEqtl$Type == "Gene",]
geneEqtl$ID = paste0(geneEqtl$gene, ";", geneEqtl$snps)
geneEqtl$numSnps = gOverlap[geneEqtl$gene]

eGenes = eGenes.DG.BHBH
eGenes$ID = paste0(eGenes$eGene, ";", eGenes$eSNPs)

## add FDR to bootstrap output
mm1 = match(eGenes$ID, geneEqtl$ID)
eGenes$FDR = geneEqtl$FDR[mm1]
mean(is.na(eGenes$FDR))
table(unique(eGenes$eGene) %in% geneEqtl$gene)
mean(!unique(eGenes$eGene) %in% geneEqtl$gene)
mean(!unique(eGenes$eGene) %in% geneEqtl$gene)

## add bootstrap output to FDR
mm2 =  match(geneEqtl$ID, eGenes$ID)
geneEqtl$corrected_beta = eGenes$corrected_beta[mm2]
geneEqtl$winners_curse = eGenes$winners_curse[mm2]
geneEqtl$eSNP_pval = eGenes$eSNP_pval[mm2]
geneEqtl$eGene_pval = eGenes$eGene_pval[mm2]

length(table(geneEqtl$gene[is.na(geneEqtl$winners_curse)]))
table(is.na(geneEqtl$winners_curse))
mean(is.na(geneEqtl$winners_curse))

g = unique(geneEqtl$gene)
table(g %in% eGenes$eGene)
mean(! g %in% eGenes$eGene)

table(g %in% geneEqtl$gene[!is.na(geneEqtl$winners_curse)])
mean(g %in% geneEqtl$gene[!is.na(geneEqtl$winners_curse)])

gTab2 = table(geneEqtl$gene, !is.na(geneEqtl$winners_curse), useNA="ifany")
gDf = data.frame(numSig = gTab2[,2], numTotal = rowSums(gTab2), 
				numTest = gOverlap[rownames(gTab2)])
gDf$sigFrac = gDf$numSig/gDf$numTotal
table(gDf$numSig==0)
boxplot(log2(gDf$numTotal) ~ gDf$numSig == 0)
t.test(log2(gDf$numTotal) ~ gDf$numSig == 0)

tapply(gDf$numTotal,  gDf$numSig == 0, quantile)

boxplot(log2(gDf$numTest) ~ gDf$numSig == 0)
plot(log2(gDf$numTotal), gDf$sigFrac)
t.test(log2(rowSums(gTab2)) ~ gTab2[,2] == 0)

## hippo checks
tt_p = table(!is.na(geneEqtl$winners_curse), geneEqtl$hippo_pvalue < 0.05,
	dnn = c("Boot", "HippoRep"))
tt_p
getOR(tt_p)
prop.table(tt_p,2)*100
prop.table(tt_p)*100
mean(geneEqtl$hippo_pvalue < 0.05,na.rm=TRUE)

ii = which(is.na(geneEqtl$winners_curse) & geneEqtl$hippo_pvalue < 0.05)
length(unique(geneEqtl$gene[ii]))

table(!is.na(geneEqtl$winners_curse), geneEqtl$hippo_pvalue < 0.01,
	dnn = c("Boot", "HippoRep"))
tt_fdr = table(!is.na(geneEqtl$winners_curse), geneEqtl$hippo_FDR < 0.05,
	dnn = c("Boot", "HippoRep"))
tt_fdr
getOR(tt_fdr)
prop.table(tt_fdr,2)*100

## reorder
geneEqtl_hOrder = geneEqtl[order(geneEqtl$gene, geneEqtl$hippo_pvalue),]
geneEqtl_hOrder = geneEqtl_hOrder[!duplicated(geneEqtl_hOrder$gene),]
table(geneEqtl_hOrder$hippo_pvalue < 0.05)
table(!is.na(geneEqtl_hOrder$winners_curse), geneEqtl_hOrder$hippo_pvalue < 0.05, 
	dnn = c("Boot", "HippoRep"))
table(!is.na(geneEqtl_hOrder$winners_curse), geneEqtl_hOrder$hippo_FDR < 0.05, 
	dnn = c("Boot", "HippoRep"))

gTab = table(geneEqtl$gene, geneEqtl$hippo_pvalue < 0.05)

geneRepl = gTab[,2]/rowSums(gTab)
geneCorr = gTab2[,2]/rowSums(gTab2)

boxplot(geneRepl ~ geneCorr == 0)
plot( )

## other checks
hist(eGenes$winners_curse)
hist(eGenes$corrected_beta/eGenes$nominal_beta)
quantile(eGenes$corrected_beta/eGenes$nominal_beta)
	
	
	geneEqtl = sigEqtl[sigEqtl$Type == "Gene",]
geneEqtl$ID = paste0(geneEqtl$gene, ";", geneEqtl$snps)

## from BootstrapQTL
load("rdas/DG_genes_Huang-bootstrap_results_Bonf-BH.Rdata",verbose=TRUE)
eGenes = eGenes.DG
eGenes$ID = paste0(eGenes$eGene, ";", eGenes$eSNPs)

## add FDR to bootstrap output
mm1 = match(eGenes$ID, geneEqtl$ID)
eGenes$FDR = geneEqtl$FDR[mm1]
eGenes$hippo_statistic = geneEqtl$hippo_statistic[mm1]
eGenes$hippo_beta = geneEqtl$hippo_beta[mm1]
eGenes$repCorrect_beta = eGenes$nominal_beta - eGenes$hippo_beta
eGenes$hippo_pvalue= geneEqtl$hippo_pvalue[mm1]
length(unique(eGenes$eGene)) # 5564

## how many were missing?
length(unique(eGenes$eGene[is.na(eGenes$FDR)])) # 415
length(unique(eGenes$eGene[is.na(eGenes$FDR)]))/
	length(unique(eGenes$eGene)) # 7.4%
length(unique(eGenes$eSNPs[is.na(eGenes$FDR)])) # 20191
length(unique(eGenes$eSNPs[is.na(eGenes$FDR)]))/
	length(unique(eGenes$eSNPs)) # 5.4%


## other checks
uIndex = !duplicated(eGenes$eGene)
plot(eGenes$repCorrect_beta, eGenes$winners_curse)
plot(eGenes$repCorrect_beta[uIndex], eGenes$winners_curse[uIndex])
plot(eGenes$repCorrect_beta[uIndex], eGenes$winners_curse[uIndex],
	xlim = c(-0.5,0.5), ylim = c(-0.2,0.2))
uSigIndex =  !duplicated(eGenes$eGene) & eGenes$hippo_pvalue < 0.05
plot(eGenes$repCorrect_beta[uSigIndex], eGenes$winners_curse[uSigIndex],
	xlim = c(-0.5,0.5), ylim = c(-0.2,0.2))

boxplot(winners_curse ~ hippo_pvalue < 0.05, data=eGenes)
boxplot(winners_curse ~ hippo_pvalue < 0.05, data=eGenes[uIndex,])

hist(eGenes$winners_curse)
hist(eGenes$corrected_beta- eGenes$nominal_beta)
quantile(eGenes$corrected_beta - eGenes$nominal_beta)

## add bootstrap output to FDR
mm2 =  match(geneEqtl$ID, eGenes$ID)
geneEqtl$corrected_beta = eGenes$corrected_beta[mm2]
geneEqtl$winners_curse = eGenes$winners_curse[mm2]
geneEqtl$eSNP_pval = eGenes$eSNP_pval[mm2]
geneEqtl$eGene_pval = eGenes$eGene_pval[mm2]

g = unique(geneEqtl$gene)
sum(!is.na(geneEqtl$winners_curse))
mean(!is.na(geneEqtl$winners_curse))
table(g %in% geneEqtl$gene[!is.na(geneEqtl$winners_curse)])
mean(g %in% geneEqtl$gene[!is.na(geneEqtl$winners_curse)])

tt_p = table(!is.na(geneEqtl$winners_curse), geneEqtl$hippo_pvalue < 0.05,
	dnn = c("Boot", "HippoRep"))
getOR(tt_p)
prop.table(tt_p,2)*100
prop.table(tt_p)*100
mean(geneEqtl$hippo_pvalue < 0.05,na.rm=TRUE)

table(!is.na(geneEqtl$winners_curse), geneEqtl$hippo_pvalue < 0.01,
	dnn = c("Boot", "HippoRep"))
tt_fdr = table(!is.na(geneEqtl$winners_curse), geneEqtl$hippo_FDR < 0.05,
	dnn = c("Boot", "HippoRep"))
getOR(tt_fdr)
prop.table(tt_fdr,2)*100
	
	
	