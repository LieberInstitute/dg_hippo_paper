####
library(GenomicRanges)
library(VennDiagram)
library(jaffelab)
library(SummarizedExperiment)
library(RColorBrewer)
library(readxl)

## load data
load("count_data/merged_dg_hippo_allSamples_n596.rda")
num_feature = c(Gene=nrow(rse_gene_joint), Exon = nrow(rse_exon_joint),
	Jxn = nrow(rse_jxn_joint), Tx = nrow(rse_tx_joint))

## load eQTLs
load("eQTL_all_SNPs_n596/eqtl_tables/mergedEqtl_output_dg_4features_fdr01_withHippo.rda")
sigEqtl$Type = factor(sigEqtl$Type, levels = c("Gene", "Exon",  "Jxn", "Tx"))
sigEqtl$EnsemblGeneID = ss(sigEqtl$EnsemblGeneID, "\\.")

dim(sigEqtl)
## split by type
sigEqtlList = split(sigEqtl, sigEqtl$Type)

length(unique(sigEqtl$EnsemblGeneID))

## number eQTLs
sapply(sigEqtlList, nrow)

## num unique features
num_eFeature = sapply(sigEqtlList, function(x) length(unique(x$gene)))
num_eFeature
 # Gene  Exon   Jxn    Tx
# 10141 67799 28319 17417
frac_eFeature = num_eFeature/num_feature


## unannoated
sigEqtlUnann = sigEqtl[sigEqtl$Class != "InGen",]
length(unique(sigEqtlUnann$gene))
length(unique(sigEqtlUnann$EnsemblGeneID))
table(sigEqtlUnann$Class[!duplicated(sigEqtlUnann$gene)])

table(unique(sigEqtlUnann$EnsemblGeneID) %in% sigEqtl$EnsemblGeneID[sigEqtl$Class == "InGen"])
## num unique genes
sapply(sigEqtlList, function(x) length(unique(x$EnsemblGeneID)))
sapply(sigEqtlList, function(x) length(unique(x$Symbol)))

## nicer venn diagram
ensPlot = sapply(sigEqtlList, function(x) unique(x$EnsemblGeneID))
v = venn.diagram(ensPlot, fill = brewer.pal(4,"Dark2"), 
	main="", main.pos = c(.5, .2), cat.cex = 1.8, cex=2.5,
	margin = 0.4, filename = NULL)
pdf("plots/eQTL_overall_geneID_vennDiagram.pdf")
grid.draw(v)
dev.off()

## effect size
sapply(sigEqtlList, function(x) quantile(abs(x$beta)))

###############################
## overall hippo replication ##
###############################
mean(sign(sigEqtl$statistic) == sign(sigEqtl$hippo_statistic),na.rm=TRUE)
mean(sign(sigEqtl$statistic) == sign(sigEqtl$hippo_statistic) & 
	sigEqtl$hippo_pvalue < 0.05,na.rm=TRUE)
mean(sign(sigEqtl$statistic) == sign(sigEqtl$hippo_statistic) & 
	sigEqtl$hippo_FDR < 0.01,na.rm=TRUE)

## maybe correlation of stats by feature
sigEqtlFeatureList =split(sigEqtl, sigEqtl$gene)
corFeature = sapply(sigEqtlFeatureList[sapply(sigEqtlFeatureList,nrow) > 1],
	function(x) cor(x$statistic, x$hippo_statistic,use="pair"))
	
d = data.frame(FeatureID = names(corFeature), cor = corFeature,
	stringsAsFactors = FALSE)
d$Type = sigEqtl$Type[match(d$FeatureID, sigEqtl$gene)]
d$EnsemblGeneID = sigEqtl$EnsemblGeneID[match(d$FeatureID, sigEqtl$gene)]

g = ggplot(d, aes(x=Type, y=cor)) + geom_violin()
ggsave(g, file="plots/eqtl_corr.pdf")

#######################
## gene specific
uGene = unique(sigEqtl$EnsemblGeneID)
geneDat = sapply(sigEqtlList, function(x) uGene %in% unlist(x$EnsemblGeneID))
geneDat = as.data.frame(geneDat)
rownames(geneDat) = uGene
geneDat$Symbol =  sigEqtl$Symbol[match(rownames(geneDat), sigEqtl$EnsemblGeneID)]

### dg specific?
sigEqtlHippo = sigEqtl[which(sigEqtl$hippo_pvalue < 0.05),]
geneDat$dgSpecific = ! rownames(geneDat) %in% sigEqtlHippo$EnsemblGeneID

geneList = apply(geneDat[,1:4], 2, which)
v = venn.diagram(geneList, fill = brewer.pal(4,"Set1"), 
	main="", main.pos = c(.5, .2), cat.cex = 1.8, cex=2,
	margin = 0.4, filename = NULL)
pdf("plots/eQTL_overall_geneSpecific_vennDiagram.pdf")
grid.draw(v)
dev.off()

##################
## tx specific?
uTx = unique(unlist(sigEqtl$gencodeTx[sigEqtl$Type != "Gene"]))
txDat = sapply(sigEqtlList[-1], function(x) uTx %in% unlist(x$gencodeTx))
txDat = as.data.frame(txDat)
rownames(txDat) = uTx
load("count_data/astellas_dg_hg38_rseTx_n263.rda")
txDat$EnsemblGeneID =  ss(rowData(rse_tx)$gene_id, "\\.")[match(uTx, rownames(rse_tx))]
txDat$Symbol =  rowData(rse_tx)$gene_name[match(uTx, rownames(rse_tx))]

### dg specific?
txDat$dgSpecific = ! uTx %in% unlist(sigEqtlHippo$gencodeTx[sigEqtlHippo$Type != "Gene"])

txList = apply(txDat[,1:3], 2, which)
v = venn.diagram(txList, fill = brewer.pal(4,"Set1")[2:4], 
	main="", main.pos = c(.5, .2), cat.cex = 1.8, cex=2,
	margin = 0.4, filename = NULL)
pdf("plots/eQTL_overall_txSpecific_vennDiagram.pdf")
grid.draw(v)
dev.off()

##############
## specific?
table(sigEqtl$inter_FDR < 0.05)
length(unique(sigEqtl$EnsemblGeneID[which(sigEqtl$inter_FDR < 0.05)]))
table(sigEqtl$inter_FDR < 0.05 & sigEqtl$hippo_pvalue > 0.05)
length(unique(sigEqtl$EnsemblGeneID[sigEqtl$inter_FDR < 0.05 & sigEqtl$hippo_pvalue > 0.05]))

## opposite dirs
table(sigEqtl$inter_FDR < 0.05 & sign(sigEqtl$hippo_statistic) != sign(sigEqtl$statistic))
sigEqtlDir = sigEqtl[which(sign(sigEqtl$hippo_statistic) != sign(sigEqtl$statistic) & 
	sigEqtl$hippo_FDR < 0.01 & sigEqtl$inter_FDR < 0.05),]
length(unique(sigEqtlDir$gene))
length(unique(sigEqtlDir$EnsemblGeneID))

## effect size?
sigEqtlHippo_FDR = sigEqtlHippo[which(sigEqtlHippo$hippo_FDR < 0.01),]
dim(sigEqtlHippo_FDR)
table(abs(sigEqtlHippo_FDR$beta) > abs(sigEqtlHippo_FDR$hippo_beta) )
mean(abs(sigEqtlHippo_FDR$beta) > abs(sigEqtlHippo_FDR$hippo_beta) )
tt = table(abs(sigEqtlHippo_FDR$beta) > abs(sigEqtlHippo_FDR$hippo_beta),
	sigEqtlHippo_FDR$pvalue <  sigEqtlHippo_FDR$hippo_pvalue,
	dnn = c("DG-bigger", "DG-sigger"))
prop.table(tt,2)

########################################
## examples of different transcripts ###
########################################

gList = split(sigEqtl, sigEqtl$Symbol)
gList = gList[sapply(gList, function(x) sum(
#######################
### plots #############
#######################

# ### libraries
# library(jaffelab)
# library(sva)
# library(recount)

# ## overlapping snps
# load("genotype_data/astellas_dg_genotype_data_n263.rda")
# snpMap$pos_hg19 = paste0(snpMap$CHR, ":", snpMap$POS)
# snpMap_dg = snpMap
# mds_dg = mds

# ######################
# ### load data ####
# ######################

# load("count_data/merged_dg_hippo_allSamples_n596.rda", verbose=TRUE)

# rse_gene = rse_gene_joint
# rse_exon = rse_exon_joint
# rse_jxn = rse_jxn_joint
# rse_tx = rse_tx_joint

# pd = colData(rse_gene)

# ## load SNP data
# load("genotype_data/merged_dg_hippo_allSamples_n596.rda")
# snpMap$pos_hg19 = paste0(snpMap$CHR, ":", snpMap$POS)

# ## drop rs10708380:150158001:TG:T (missing info in snpMap (and dbSNP))
# snpInd = which(rownames(snpMap) == "rs10708380:150158001:TG:T")
# snpMap = snpMap[-snpInd,]
# snp = snp[-snpInd,]
# pd$Dx = factor(pd$Dx,
	# levels = c("Control", "Schizo", "Bipolar", "MDD"))

# mod = model.matrix(~Dx + Sex + as.matrix(mds[,1:5]) + Region, data = pd)
# colnames(mod)[6:10] = colnames(mds)[1:5]

# ## extract pd and rpkms
# geneRpkm = recount::getRPKM(rse_gene, "Length")
# exonRpkm = recount::getRPKM(rse_exon, "Length")
# rowRanges(rse_jxn)$Length <- 100
# jxnRp10m = recount::getRPKM(rse_jxn, "Length")
# txTpm = assays(rse_tx)$tpm

# ## pcs
# load("eQTL_GWAS_riskSNPs_n596/rdas/pcs_4features_interaction_combined_regions_596.rda")

# geneExprsClean = cleaningY(log2(geneRpkm[rownames(geneRpkm) %in% sigEqtl$gene,]+1), cbind(mod, genePCs), P=1)
# exonExprsClean = cleaningY(log2(exonRpkm[rownames(exonRpkm) %in% sigEqtl$gene,]+1), cbind(mod, exonPCs), P=1)
# jxnExprsClean = cleaningY(log2(jxnRp10m[rownames(jxnRp10m) %in% sigEqtl$gene,]+1), cbind(mod, jxnPCs), P=1)
# txExprsClean = cleaningY(log2(txTpm[rownames(txTpm) %in% sigEqtl$gene,]+1), cbind(mod, txPCs), P=1)
# exprsClean = rbind(geneExprsClean, exonExprsClean, jxnExprsClean, txExprsClean)

# ## snp data
# snp2 = snp[rownames(snp) %in% sigEqtl$snps,]
# snpMap2 = snpMap[rownames(snpMap) %in% sigEqtl$snps,]
# save(exprsClean, snp2, snpMap2, pd, file = "rdas/cleaned_exprs_data_twoData_eQTL.rda")

### load data back in
load("rdas/cleaned_exprs_data_twoData_eQTL.rda")

## which to plot??
topIntEqtls = sigEqtl[order(sigEqtl$inter_pvalue),]
topIntEqtls = topIntEqtls[!is.na(topIntEqtls$inter_FDR),]
topIntEqtls = topIntEqtls[!duplicated(topIntEqtls$Symbol),]

## filter
ssnp = snp2[topIntEqtls$snps,]
ssnpMap = snpMap2[topIntEqtls$snps,]
e = exprsClean[topIntEqtls$gene,]
r = pd$Region
r[r=="DentateGyrus"] = "DG-GCL"
r = factor(r, levels = c("HIPPO", "DG-GCL"))

## update exon ID
load("rdas/exon_names_hg19_hg38.rda",verbose=TRUE)
topIntEqtls$FeatLabel = topIntEqtls$gene
mmExon = match(topIntEqtls$FeatLabel, hg38_exons$hg38_eID)
topIntEqtls$FeatLabel[!is.na(mmExon)] = hg38_exons$gencode_exonID[mmExon[!is.na(mmExon)]]

topIntEqtls$SnpLabel = ssnpMap$name[match(topIntEqtls$snps, ssnpMap$SNP)]

pdf("plots/overall_eqtl_boxplots_interaction.pdf",w=8)
par(mar=c(5,6,4,2), cex.axis=2,cex.lab=2, cex.main = 2)
palette(brewer.pal(5,"Set1"))
# for(i in 1:5) { # for testing
for(i in 1:100) {
	s = as.numeric(ssnp[i,])
	s = gsub(0, paste0(ssnpMap$ALT[i],ssnpMap$ALT[i]), s)
	s = gsub(1, paste0(ssnpMap$ALT[i],ssnpMap$COUNTED[i]), s)
	s = gsub(2, paste0(ssnpMap$COUNTED[i],ssnpMap$COUNTED[i]), s)
	
	l = paste0(s, ":", r)
	l = factor(l, levels = paste0(rep(
		c(paste0(ssnpMap$ALT[i],	ssnpMap$ALT[i]),
		paste0(ssnpMap$ALT[i], ssnpMap$COUNTED[i]), 
		paste0(ssnpMap$COUNTED[i], 	ssnpMap$COUNTED[i])), times=2),
			":", rep(c("DG-GCL", "HIPPO"), each=3)))
	boxplot(e[i,] ~ l, outline=FALSE, 
		ylim= range(e[i,])*c(0.95,1.05), 
		ylab= topIntEqtls$FeatLabel[i],
		xlab=topIntEqtls$SnpLabel[i], main = topIntEqtls$Symbol[i],
		names = ss(levels(l),":"))
	points(e[i,] ~ jitter(as.numeric(l),amount=0.1),
		pch = 21, bg=r,cex=1.4)
	legend("topleft", paste0("p=", 
		signif(topIntEqtls$pvalue[i], 3)),cex=1.8, bty="n")
	legend("topright", paste0("p=", 
		signif(topIntEqtls$hippo_pvalue[i], 3)), cex=1.8, bty="n")
	text(x=2, y = min(e[i,])*0.97, "DG-GCL", cex=2)
	text(x=5, y = min(e[i,])*0.97, "HIPPO", cex=2)
	abline(v=3.5,lty=2)
}
dev.off()

pdf("plots/overall_eqtl_boxplots_interaction_sparseLabel.pdf",w=8,useDingbats=FALSE)
par(mar=c(5,6,4,2), cex.axis=2,cex.lab=2, cex.main = 2)
palette(brewer.pal(5,"Set1"))
# for(i in 1:5) { # for testing
for(i in 1:100) {
	s = as.numeric(ssnp[i,])
	s = gsub(0, paste0(ssnpMap$ALT[i],ssnpMap$ALT[i]), s)
	s = gsub(1, paste0(ssnpMap$ALT[i],ssnpMap$COUNTED[i]), s)
	s = gsub(2, paste0(ssnpMap$COUNTED[i],ssnpMap$COUNTED[i]), s)
	
	l = paste0(s, ":", r)
	l = factor(l, levels = paste0(rep(
		c(paste0(ssnpMap$ALT[i],	ssnpMap$ALT[i]),
		paste0(ssnpMap$ALT[i], ssnpMap$COUNTED[i]), 
		paste0(ssnpMap$COUNTED[i], 	ssnpMap$COUNTED[i])), times=2),
			":", rep(c("DG-GCL", "HIPPO"), each=3)))
	boxplot(e[i,] ~ l, outline=FALSE, 
		ylim= range(e[i,])*c(0.95,1.05), 
		ylab= "Expression (log2)",
		xlab=topIntEqtls$SnpLabel[i],
		names = ss(levels(l),":"))
	points(e[i,] ~ jitter(as.numeric(l),amount=0.1),
		pch = 21, bg=r,cex=1.4)
	legend("topleft", paste0("p=", 
		signif(topIntEqtls$pvalue[i], 3)),cex=1.8, bty="n")
	legend("topright", paste0("p=", 
		signif(topIntEqtls$hippo_pvalue[i], 3)), cex=1.8, bty="n")
	legend("top", topIntEqtls$Symbol[i], cex=2,bty="n")
	text(x=2, y = min(e[i,])*0.97, "DG-GCL", cex=2)
	text(x=5, y = min(e[i,])*0.97, "HIPPO", cex=2)
	abline(v=3.5,lty=2)
}
dev.off()


######
## check dopamine results
dop = read_excel("tables/41593_2018_223_MOESM9_ESM.xls",skip=2)
colnames(dop) = gsub(" ", "_", colnames(dop))
dop = as.data.frame(dop)
dop = dop[grep("^rs", dop[,1]),]

## filter down
sigEqtl$snpRsNum = snpMap2$name[match(sigEqtl$snps, snpMap2$SNP)]
sigEqtlDop = sigEqtl[sigEqtl$snpRsNum %in% dop$SNP_ID..1,]

length(unique(sigEqtlDop$snpRsNum))

sigEqtlDop[which(sigEqtlDop$snpRsNum == "rs642803"),]
sigEqtlDop[which(sigEqtlDop$snpRsNum == "rs393152"),]

Indexes = lapply(dop$SNP_ID..1, grep, x=sigEqtl$snps)
