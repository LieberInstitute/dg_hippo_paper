####
library(GenomicRanges)
library(VennDiagram)
library(jaffelab)
library(SummarizedExperiment)
library(RColorBrewer)

## load eQTLs
load("eQTL_GWAS_riskSNPs_n596/eqtl_tables/mergedEqtl_output_withHippo_withDlpfc_fdr05any.rda")

sigEqtl$Type = factor(sigEqtl$Type, levels = c("Gene", "Exon",  "Jxn", "Tx"))
sigEqtl$EnsemblGeneID = ss(sigEqtl$EnsemblGeneID, "\\.")

## write out
sigEqtlOut = as.data.frame(sigEqtl)
sigEqtlOut$gencodeTx = sapply(sigEqtlOut$gencodeTx , paste, collapse=";")
write.csv(sigEqtlOut, file = gzfile("tables/suppTable_sczdEqtl.csv.gz"), row.names=FALSE)

####
## get the region-specific ones
dgEqtl = sigEqtl[sigEqtl$FDR < 0.01,]
hippoEqtl = sigEqtl[sigEqtl$hippo_FDR < 0.01,]
dlpfcEqtl = sigEqtl[which(sigEqtl$dlpfc_FDR < 0.01),]

length(unique(sigEqtl$Index_Name)) # 136, @ FDR< 5%

length(unique(dgEqtl$Index_Name)) # 78
length(unique(hippoEqtl$Index_Name)) # 82
length(unique(dlpfcEqtl$Index_Name)) # 94
sum(unique(dgEqtl$Index_Name) %in% hippoEqtl$Index_Name) # 60
sum(! unique(dgEqtl$Index_Name) %in% c(hippoEqtl$Index_Name,dlpfcEqtl$Index_Name)) # 8

table(unique(dlpfcEqtl$Index_Name) %in% unique(c(dgEqtl$Index_Name, hippoEqtl$Index_Name)))
length(unique(c(dgEqtl$Index_Name, dlpfcEqtl$Index_Name, hippoEqtl$Index_Name)))

## make venn diagram

###########
## dg only

dgEqtlOnly = dgEqtl[! dgEqtl$Index_Name %in% hippoEqtl$Index_Name,]
dgEqtlOnly = dgEqtlOnly[order(dgEqtlOnly$pvalue),]
as.data.frame(dgEqtlOnly[dgEqtlOnly$Distance == 0,])

options(width=150)
colName = c("gene", "Symbol", "pvalue", "hippo_pvalue","dlpfc_pvalue",
	"inter_pvalue", "Index_Name", "R_squared")
as.data.frame(dgEqtlOnly)[dgEqtlOnly$Distance == 0,colName]
dgEqtlOnly[dgEqtlOnly$Index_Name == "rs2007044:2344960:A:G",colName]

dgEqtlOnlyTop = dgEqtlOnly[!duplicated(dgEqtlOnly$Index_Name),]
dgEqtlOnlyTop$lab = paste0("R^2=",signif(dgEqtlOnlyTop$R_squared,2), ", ",dgEqtlOnlyTop$Type, 
	": DGp=", signif(dgEqtlOnlyTop$pvalue,3), " Hp=",
	signif(dgEqtlOnlyTop$hippo_pvalue, 3), " PFCp=",signif(dgEqtlOnlyTop$dlpfc_pvalue, 3))
dgEqtlOnlyTop$ind = which(!duplicated(dgEqtlOnly$Index_Name))
as.data.frame(dgEqtlOnlyTop[,c("Symbol", "lab","ind")])

## check with DLPFC
byLocus = split(sigEqtl, sigEqtl$Index_Name)

#############################
########## plots ############
#############################

# ### libraries
# library(SummarizedExperiment)
# library(jaffelab)
# library(MatrixEQTL)
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
# save(exprsClean, snp2, snpMap2, pd, file = "rdas/cleaned_exprs_data_twoData_eQTL_PGConly.rda")
load("rdas/cleaned_exprs_data_twoData_eQTL_PGConly.rda")

### save grm3 for jooheon
grm3 = data.frame(log2RPKMplus1_cleaned = exprsClean["ENSG00000198822.10",],
		rs6943762 = snp2["rs6943762",], Region = pd$Region,	stringsAsFactors=FALSE)	
grm3$rs6943762_Geno = ifelse(grm3$rs6943762 == 0, "TT", ifelse(grm3$rs6943762 == 1, "CT","CC"))
boxplot(grm3$log2RPKMplus1_cleaned ~ grm3$rs6943762_Geno*grm3$Region)
write.csv(grm3, file = "grm3_data_for_jooheon.csv")

############################
## make a bunch of plots ###
############################

ssnp = snp2[dgEqtlOnly$snps,]
ssnpMap = snpMap2[dgEqtlOnly$snps,]
e = exprsClean[dgEqtlOnly$gene,]
r = pd$Region
r[r=="DentateGyrus"] = "DG-GCL"
r = factor(r, levels = c("HIPPO", "DG-GCL"))

## update exon ID
load("rdas/exon_names_hg19_hg38.rda",verbose=TRUE)
dgEqtlOnly$FeatLabel = dgEqtlOnly$gene
mmExon = match(dgEqtlOnly$FeatLabel, hg38_exons$hg38_eID)
dgEqtlOnly$FeatLabel[!is.na(mmExon)] = hg38_exons$gencode_exonID[mmExon[!is.na(mmExon)]]

dgEqtlOnly$SnpLabel = dgEqtlOnly$Proxy_Name
dgEqtlOnly$SnpLabel[grep("^rs", dgEqtlOnly$SnpLabel)] = ss(dgEqtlOnly$SnpLabel[grep("^rs", dgEqtlOnly$SnpLabel)], ":")
mainTxt = paste0("Index: ", ss(dgEqtlOnly$Index_Name, ":"), 
	" (R^2=",signif(dgEqtlOnly$R_squared,3),")\n",
	dgEqtlOnly$Symbol, " (", dgEqtlOnly$Type, ")")

pdf("plots/gwas_eqtl_boxplots_dgOnly.pdf",w=8)
par(mar=c(5,6,4,2), cex.axis=2,cex.lab=2, cex.main = 2)
palette(brewer.pal(5,"Set1"))
# for(i in 1:5) { # for testing
for(i in 1:nrow(e)) {
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
		ylab= dgEqtlOnly$FeatLabel[i],
		xlab=dgEqtlOnly$SnpLabel[i], main = mainTxt[i],
		names = ss(levels(l),":"))
	points(e[i,] ~ jitter(as.numeric(l),amount=0.1),
		pch = 21, bg=factor(r),cex=1.4)
	legend("topleft", paste0("p=", 
		signif(dgEqtlOnly$pvalue[i], 3)),cex=1.8, bty="n")
	legend("topright", paste0("p=", 
		signif(dgEqtlOnly$hippo_pvalue[i], 3)), cex=1.8, bty="n")
	text(x=2, y = min(e[i,])*0.97, "DG-GCL", cex=2)
	text(x=5, y = min(e[i,])*0.97, "HIPPO", cex=2)
	abline(v=3.5,lty=2)
}
dev.off()

	dim(sigEqtl)
## split by type
sigEqtlList = split(sigEqtl, sigEqtl$Type)

length(unique(sigEqtl$EnsemblGeneID))

## number eQTLs
sapply(sigEqtlList, nrow)

## num unique features
sapply(sigEqtlList, function(x) length(unique(x$gene)))

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