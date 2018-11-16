###

## check nestin and DCX
ii = match(c("DCX", "NES") ,rowData(rse_gene_joint)$Symbol)
checkTwo  = t(getRPKM(rse_gene_joint, "Length")[ii,])

## ratio
dat = colData(rse_gene_joint)
dat= cbind(dat,checkTwo)
dcx_ratio = tapply(dat$ENSG00000077279.17, dat$BrNum, function(x) log2(x[1]+1) - log2(x[2]+1))
nes_ratio = tapply(dat$ENSG00000132688.10, dat$BrNum, function(x) log2(x[1]+1) - log2(x[2]+1))

############################
## explore, this becomes a figure for QC
bIndexes = splitit(rse_gene_joint$BrNum)

pdf("qcChecks/quality_by_cellType.pdf")
par(mar=c(5,6,3,2), cex.axis=1.5,cex.lab=2,cex.main=2)
palette(brewer.pal(5,"Set2"))
boxplot(mitoRate ~ Dataset, data=colData(rse_gene_joint),
	ylab = "chrM Map Rate", outline = FALSE, 
	ylim = range(rse_gene_joint$mitoRate)) # ribozero regular...hippo > dg
xx = jitter(as.numeric(rse_gene_joint$Dataset),amount=0.15)
for(i in seq(along=bIndexes)) {
	lines(mitoRate ~ xx, data=colData(rse_gene_joint),
		subset=bIndexes[[i]], col ="grey",lwd=0.4)
}
points(	mitoRate ~ xx,	data=colData(rse_gene_joint),
	bg = factor(rse_gene_joint$LibraryType), pch =21)
legend("topright", c("Gold", "HMR"), col = 1:2, pch = 15,cex=2,nc=1)
boxplot(rRNA_rate ~ Dataset, data=colData(rse_gene_joint),
	ylab = "rRNA Gene Map Rate",outline=FALSE,
	ylim = range(rse_gene_joint$rRNA_rate)) # ribozero regular...hippo > dg
for(i in seq(along=bIndexes)) {
	lines(rRNA_rate ~ xx, data=colData(rse_gene_joint),
		subset=bIndexes[[i]], col ="grey",lwd=0.4)
}
points(	rRNA_rate ~ xx,	data=colData(rse_gene_joint),
	bg = factor(rse_gene_joint$LibraryType), pch =21)

boxplot(RIN ~ Dataset, data=colData(rse_gene_joint), ylab = "RIN",
	outline=FALSE, ylim = range(rse_gene_joint$RIN, na.rm=TRUE)) # hippo > dg
for(i in seq(along=bIndexes)) {
	lines(RIN ~ xx, data=colData(rse_gene_joint),
		subset=bIndexes[[i]], col ="grey",lwd=0.4)
}
points(	RIN ~ xx,	data=colData(rse_gene_joint),
	bg = factor(rse_gene_joint$LibraryType), pch =21)


boxplot(overallMapRate ~ Dataset, data=colData(rse_gene_joint),
	ylab = "Overall Map Rate",
	outline=FALSE, ylim = range(rse_gene_joint$overallMapRate)) # hippo > dg
for(i in seq(along=bIndexes)) {
	lines(overallMapRate ~ xx, data=colData(rse_gene_joint),
		subset=bIndexes[[i]], col ="grey",lwd=0.4)
}
points(	overallMapRate ~ xx,	data=colData(rse_gene_joint),
	bg = factor(rse_gene_joint$LibraryType), pch =21)

	
boxplot(totalAssignedGene ~ Dataset, data=colData(rse_gene_joint),
	ylab = "Exonic Mapping Rate", outline=FALSE, 
	ylim = range(rse_gene_joint$totalAssignedGene)) # hippo > dg
for(i in seq(along=bIndexes)) {
	lines(totalAssignedGene ~ xx, data=colData(rse_gene_joint),
		subset=bIndexes[[i]], col ="grey",lwd=0.4)
}
points(	totalAssignedGene ~ xx,	data=colData(rse_gene_joint),
	bg = factor(rse_gene_joint$LibraryType), pch =21)

dev.off()

