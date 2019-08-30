library(limma)

## load DG-GCL
load("twas/DentateGyrus/gene/heritability/hsq_info.Rdata")
hsq_dg = hsq_info

## load DLPFC/HIPPO
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/rda/hsq_genes.Rdata")
hsq_hippo = hsq$HIPPO
hsq_dlpfc = hsq$DLPFC
identical(hsq_hippo$ID, hsq_dlpfc$ID) # TRUE

## merge
hsq_all = cbind(hsq_hippo[,c("ID", "geneID", "hsq", "hsq.se", "hsq.pv")], 
	hsq_dlpfc[,c("hsq", "hsq.se", "hsq.pv")])
colnames(hsq_all)[3:5] = paste0(colnames(hsq_all)[3:5], "_hippo")
colnames(hsq_all)[6:8] = paste0(colnames(hsq_all)[6:8], "_dlpfc")

## add in dg
hsq_dg = hsq_dg[,c("ID", "geneID", "hsq", "hsq.se", "hsq.pv")]
colnames(hsq_dg)[3:5] = paste0(colnames(hsq_dg)[3:5], "_dg")

## merge
hsq_all = merge(hsq_dg, hsq_all, sort=FALSE, all=TRUE)

## make NAs = 0
hsq_all$hsq_dg[is.na(hsq_all$hsq_dg)] = -0.25
hsq_all$hsq_hippo[is.na(hsq_all$hsq_hippo)] = -0.25
hsq_all$hsq_dlpfc[is.na(hsq_all$hsq_dlpfc)] = -0.25

# number sig
hsq_all$hsq.bh_dg = p.adjust(hsq_all$hsq.pv_dg,"fdr")
hsq_all$hsq.bh_hippo = p.adjust(hsq_all$hsq.pv_hippo,"fdr")
hsq_all$hsq.bh_dlpfc = p.adjust(hsq_all$hsq.pv_dlpfc,"fdr")

fdrs = hsq_all[,grep("bh", colnames(hsq_all))]
colnames(fdrs) = c("DG-GCL", "HIPPO", "DLPFC")
vennDiagram(vennCounts(fdrs < 0.05))
colSums(fdrs < 0.05,na.rm=TRUE)

## filter
hsq_sig = hsq_all[hsq_all$hsq.bh_dg < 0.05,]
table(hsq_sig$hsq.pv_hippo < 0.05)
mean(hsq_sig$hsq.pv_hippo > 0.05,na.rm=TRUE)
table(hsq_sig$hsq.pv_dlpfc < 0.05)

hsq_dlpfc_sig =  hsq_all[hsq_all$hsq.bh_dlpfc < 0.05,]
table(hsq_dlpfc_sig$hsq.pv_hippo < 0.05)
mean(hsq_dlpfc_sig$hsq.pv_hippo > 0.05,na.rm=TRUE)
 #########
## MA ###
#########
pdf("plots/heritability_comparisons_MA.pdf",h=6,w=8)
par(mar = c(5,6,4,2), cex.lab = 2, cex.axis = 2,	cex.main = 2)

## DLPFC vs HIPPO
plot( x = (hsq_all$hsq_dlpfc + hsq_all$hsq_hippo)/2,
		y = hsq_all$hsq_dlpfc - hsq_all$hsq_hippo,
		xlab = '(DLPFC + HIPPO) / 2', bg = "grey",
		ylim = c(-0.7,0.7), xlim = c(-0.1, 1),
		ylab = 'DLPFC - HIPPO',	pch = 21,  
		main = 'Heritability (cis +/- 500kb)')
		
# DG-GCL and HIPPO		
plot( x = (hsq_all$hsq_dg+ hsq_all$hsq_hippo)/2,
		y = hsq_all$hsq_dg - hsq_all$hsq_hippo,
		xlab = '(DG-GCL + HIPPO) / 2',
		ylab = 'DG-GCL - HIPPO',bg = "grey",
		ylim = c(-0.7,0.7), xlim = c(-0.1, 1),
		pch = 21,  main = 'Heritability (cis +/- 500kb)')
dev.off()
		
##############
## scatter ###
##############
pdf("plots/heritability_comparisons_scatter.pdf")
par(mar = c(5,6,4,2), cex.lab = 2, cex.axis = 2,	cex.main = 2)

## DLPFC vs HIPPO
plot( x = hsq_all$hsq_dlpfc, y = hsq_all$hsq_hippo,
		xlab = 'DLPFC' ,  ylab = 'HIPPO', pch = 21, bg="grey",
		xlim = c(-0.1, 1),ylim = c(-0.1, 1),
		main = 'Heritability (cis +/- 500kb)')		
abline(a=0.2, b=1, lty=2,col="blue",lwd=2)
abline(a=0.1, b=1, lty=1,col="blue",lwd=2)
abline(a=-0.1, b=1, lty=1,col="blue",lwd=2)
abline(a=-0.2, b=1, lty=2,col="blue",lwd=2)

## DG-GCL vs HIPPO
plot( x = hsq_all$hsq_dg, y = hsq_all$hsq_hippo,
		xlab = 'DG-GCL' ,  ylab = 'HIPPO',  pch = 21, bg="grey",
		xlim = c(-0.1, 1),ylim = c(-0.1, 1),
		main = 'Heritability (cis +/- 500kb)')
abline(a=0.2, b=1, lty=2,col="blue",lwd=2)
abline(a=0.1, b=1, lty=1,col="blue",lwd=2)
abline(a=-0.1, b=1, lty=1,col="blue",lwd=2)
abline(a=-0.2, b=1, lty=2,col="blue",lwd=2)
dev.off()

## counts
hsq_all$dg_hippo_diff = hsq_all$hsq_dg - hsq_all$hsq_hippo
hsq_all$dlpfc_hippo_diff = hsq_all$hsq_dlpfc - hsq_all$hsq_hippo
table(abs(hsq_all$dg_hippo_diff) > 0.1)
table(abs(hsq_all$dlpfc_hippo_diff) > 0.1)
table(abs(hsq_all$dg_hippo_diff) > 0.2)
table(abs(hsq_all$dlpfc_hippo_diff) > 0.2)
