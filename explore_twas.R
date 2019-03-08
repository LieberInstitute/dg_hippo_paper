library('tibble')
library('sessioninfo')
library('purrr')
library('SummarizedExperiment')
library('dplyr')
library('limma')
library('readr')
library('jaffelab')

## Load twas data (from read_twas.R)
load('rdas/twas.Rdata', verbose = TRUE)

## filtered counts
load("count_data/merged_dg_hippo_allSamples_n596.rda")
colnames(rse_deg_joint) = colnames(rse_gene_joint)

rse = list(rse_gene_joint, rse_exon_joint,
	rse_jxn_joint, rse_tx_joint)
names(rse) <- unique(twas$all$feature)

## Add the gene id
twas_exp <- map(twas, function(tw) {
    ## For testing the inner part of the function
    # tw <- twas$included
    
    by_feat <- split(tw, tw$feature)
    ## Make sure it's in the right order
    if(!identical(names(by_feat), names(rse))) {
        message(paste(Sys.time(), 'fixing order'))
        by_feat <- by_feat[names(rse)]
    }
    
    ## Now add the gene gencode ID and symbol
    result <- pmap_dfr(
        list(by_feat, rse, names(rse)),
        function(info, rs, feature) {
            
            ## Find the appropriate variables
            gene_var <- case_when(
                feature == 'gene' ~ 'gencodeID',
                feature == 'exon' ~ 'gencodeID',
                feature == 'jxn' ~ 'newGeneID',
                feature == 'tx' ~ 'gene_id'
            )
            symbol_var <- case_when(
                feature == 'gene' ~ 'Symbol',
                feature == 'exon' ~ 'Symbol',
                feature == 'jxn' ~ 'newGeneSymbol',
                feature == 'tx' ~ 'gene_name'
            )
            
            ## Match by id
            m <- match(info$ID, names(rowRanges(rs)))
            stopifnot(!is.na(m))
            
            ## Add the gene id/symbol
            info$geneid <- mcols(rowRanges(rs))[[gene_var]][m]
            info$genesymbol <- mcols(rowRanges(rs))[[symbol_var]][m]
            
            ## Done
            return(info)   
        }
    )
    
    return(result)
})
names(twas_exp) <- names(twas)

## Save the data for later use
save(twas_exp, file = 'rdas/twas_exp.Rdata')

######################
## MORE ANNOTATION ###
######################
tt <- twas_exp$all
## Drop TWAS NA p-values
tt <- tt[!is.na(tt$TWAS.P), ]
## Focus on CLOZUK+PGC2 (psycm) GWAS
tt <- tt[which(tt$type == "psycm"),]

## Add GWAS p-value and OR from the original sumstats file
original <- read_tsv('twas/psycm/clozuk_pgc2.meta.sumstats.txt')
original$CHR[original$CHR == 23] <- 'X'
original$hg19_pos <- with(original, paste0(CHR, ':', BP))

snpmap <- read_tsv('twas/psycm/clozuk_pgc2.meta.reformatted.sumstats_hg38_ourname',
    col_types = cols(
      SNP = col_character(),
      A1 = col_character(),
      A2 = col_character(),
      Z = col_double(),
      N = col_double(),
      chr = col_character(),
      basepair = col_double(),
      basepairhg19 = col_double(),
      originalSNP = col_character()
    )
)
snpmap$chrpos_hg19 = paste0(snpmap$chr, ":", snpmap$basepairhg19)


## Load big snpMap table
bim_hg38 = read_delim("/dcl01/ajaffe/data/lab/dg_hippo_paper/twas/filter_data/LIBD_Brain_Illumina_h650_1M_Omni5M_Omni2pt5_Macrogen_imputed_run2_LDfiltered_DentateGyrus.bim", 
	delim = "\t", col_names=  FALSE)
bim_hg19 = read_delim("/dcl01/ajaffe/data/lab/dg_hippo_paper/twas/filter_data/LIBD_Brain_Illumina_h650_1M_Omni5M_Omni2pt5_Macrogen_imputed_run2_LDfiltered_DentateGyrus.bim.original", 
	delim = "\t", col_names=  FALSE)
colnames(bim_hg19) = colnames(bim_hg38) = c("CHR", "SNP", "CM", "POS", "A1", "A2")

## add chr and pos
bim_hg19$chrpos_hg19 = paste0(bim_hg19$CHR, ":", bim_hg19$POS)
bim_hg38$chrpos_hg38 = paste0(bim_hg38$CHR, ":", bim_hg38$POS)
bim_hg38$chrpos_hg19 = bim_hg19$chrpos_hg19[match(bim_hg38$SNP, bim_hg19$SNP)]

## Match to the big snpMap table
m_to_fullMap <- match(tt$BEST.GWAS.ID, bim_hg38$SNP)
stopifnot(!any(is.na(m_to_fullMap)))
m_to_fullMap_qtl <- match(tt$EQTL.ID, bim_hg38$SNP)
stopifnot(!any(is.na(m_to_fullMap_qtl)))

## Add pos hg19 and pos hg38
tt$BEST.GWAS.chrpos_hg19 <- bim_hg38$chrpos_hg19[m_to_fullMap]
tt$BEST.GWAS.chrpos_hg38 <- bim_hg38$chrpos_hg38[m_to_fullMap]

tt$EQTL.chrpos_hg19 <- bim_hg38$chrpos_hg19[m_to_fullMap_qtl]
tt$EQTL.chrpos_hg38 <- bim_hg38$chrpos_hg38[m_to_fullMap_qtl]

m_to_map <- match(tt$BEST.GWAS.ID, snpmap$SNP)
table(is.na(m_to_map))
## Hm... I'm not sure why some are NAs
#  FALSE   TRUE
#  41212   467
m_to_map <- match(tt$BEST.GWAS.chrpos_hg19, snpmap$chrpos_hg19)
table(is.na(m_to_map))

print(tt[head(which(is.na(m_to_map))), ], width = 200)

## Hm....
m_to_map_qtl <- match(tt$EQTL.ID, snpmap$SNP)
table(is.na(m_to_map_qtl))
# FALSE  TRUE
# 41392   287

addmargins(table(
    'By BEST.GWAS.ID' = is.na(m_to_map),
    'By EQTL.ID' = is.na(m_to_map_qtl)
))
               # By EQTL.ID
# By BEST.GWAS.ID FALSE  TRUE   Sum
          # FALSE 40964   248 41212
          # TRUE    428    39   467
          # Sum   41392   287 41679
		  
## Well, after that it all looks ok
m_to_ori <- match(snpmap$originalSNP[m_to_map], original$SNP)
table(is.na(m_to_ori))
# FALSE  TRUE
# 41212   467

## Hm... it's odd that the same number don't match by either name or chr position
m_to_ori2 <- match(tt$BEST.GWAS.chrpos_hg19, original$hg19_pos)
table(is.na(m_to_ori), is.na(m_to_ori2))
        # FALSE  TRUE
  # FALSE 41212     0
  # TRUE      0   467
## Supplement m_to_ori (by name) with the matching by chr and position in hg19
m_to_ori[is.na(m_to_ori)] <- m_to_ori2[is.na(m_to_ori)]
table(is.na(m_to_ori))
#  FALSE   TRUE
# 127238     13

m_to_map_qtl2 <- match(tt$EQTL.chrpos_hg19, original$hg19_pos)
table(is.na(m_to_map_qtl), is.na(m_to_map_qtl2))
        # FALSE  TRUE
  # FALSE 41392     0
  # TRUE      0   287



## Lets get the originally reported summarized p-values
BEST.GWAS.P <- original$P[m_to_ori]
## This is how the calculate the p-values displayed by FUSION-TWAS
# https://github.com/gusevlab/fusion_twas/blob/master/FUSION.post_process.R#L641
BEST.GWAS.P.computed <- 2*(pnorm( abs(tt$BEST.GWAS.Z ) , lower.tail=F ))

## Ok, assign them to our table
tt$BEST.GWAS.P <- original$P[m_to_ori]
tt$BEST.GWAS.OR <- original$OR[m_to_ori]
tt$BEST.GWAS.SE <- original$SE[m_to_ori]
tt$BEST.GWAS.P.computed <- 2*(pnorm( abs(tt$BEST.GWAS.Z ) , lower.tail=F ))
tt$EQTL.P.computed <- 2*(pnorm( abs(tt$EQTL.GWAS.Z ) , lower.tail=F ))

## Compute FDR by region for each feature 4 features
tt <- map_dfr(split(tt, tt$region), function(reg) {
    res <- map_dfr(split(reg, reg$feature), function(reg_feat) {
        reg_feat$TWAS.FDR <- p.adjust(reg_feat$TWAS.P, 'fdr')
        reg_feat$TWAS.Bonf <- p.adjust(reg_feat$TWAS.P, 'bonf')
        return(reg_feat)
    })
    return(res[order(res$TWAS.P), ])
})

## Add raggr data
## Code based on https://github.com/LieberInstitute/brainseq_phase2/blob/master/eQTL_GWAS_riskSNPs/create_eqtl_table_indexInfo.R

## risk loci from PGC paper
indexLoci <- read.csv("eQTL_GWAS_riskSNPs_n596/pgc_riskLoci.csv", stringsAsFactors=FALSE)
indexLoci$hg19POS = paste0(indexLoci$Chromosome, ":", indexLoci$snp_pos_hg19)

## risk loci from PGC paper + rAggr proxy markers
riskLoci <- read.csv("eQTL_GWAS_riskSNPs_n596/rAggr_results_179.csv", stringsAsFactors=FALSE)
colnames(riskLoci) = gsub("\\.", "_", colnames(riskLoci))
riskLoci$SNP1_rs = ss(riskLoci$SNP1_Name, ":")
riskLoci$SNP2_rs = ss(riskLoci$SNP2_Name, ":")

## match index
riskLoci$matchLocus = match(riskLoci$SNP1_rs, indexLoci$Index_SNP_dbSNP_b141_ )

riskLoci$SNP1_chrpos_hg19 = paste0(riskLoci$SNP1_Chr, ":", riskLoci$SNP1_Pos) 
riskLoci$SNP2_chrpos_hg19 = paste0(riskLoci$SNP2_Chr, ":", riskLoci$SNP2_Pos) 

mm = match(riskLoci$SNP2_chrpos_hg19, snpmap$chrpos_hg19)

#### annotate
tt$Status = "Other"
tt$Status[tt$BEST.GWAS.chrpos_hg19 %in% riskLoci$SNP2_chrpos_hg19] = "Proxy"
tt$Status[tt$BEST.GWAS.chrpos_hg19 %in% riskLoci$SNP1_chrpos_hg19] = "Index"

#### match back to risk locus
tt$riskLocus = riskLoci$matchLocus[match(tt$BEST.GWAS.chrpos_hg19, riskLoci$SNP2_chrpos_hg19)]
length(table(tt$riskLocus))
tt$distProxyToIndex = riskLoci$Distance[match(tt$BEST.GWAS.chrpos_hg19, riskLoci$SNP2_chrpos_hg19)]
table(!is.na(tt$riskLocus), tt$Status)
        # Index Other Proxy
  # FALSE     0 39377     0
  # TRUE    756     0  1546
quantile(abs(tt$distProxyToIndex),na.rm=TRUE)
    # 0%    25%    50%    75%   100%
     # 0      0   4799  22462 251788

## rename for hippo compare
tt_DG = tt
save(tt_DG, riskLoci, file = 'rdas/twas_exp_DG_GWASannotated.Rdata') 

###################
## read in hippo ##
###################
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/rda/tt_objects.Rdata", verbose=TRUE)
tt_Hippo = tt[tt$region == "HIPPO",]
tt_Hippo$riskLocus = riskLoci$matchLocus[match(tt_Hippo$BEST.GWAS.pos_hg19, riskLoci$SNP2_chrpos_hg19)]

table(tt_Hippo$BEST.GWAS.status)

### merge 
colnames(tt_DG)[c(38,40)] = c("BEST.GWAS.status","BEST.GWAS.indexSNP_distance")
tt_DG$Region = "DG"
tt_Hippo$Region = "HIPPO"
nam = colnames(tt_DG)[colnames(tt_DG) %in% colnames(tt_Hippo)]

tt_both = rbind(tt_DG[,nam], tt_Hippo[,nam])

## if GWAS significant, make into proxy
tt_both$BEST.GWAS.status[tt_both$BEST.GWAS.P < 5e-8] = "Proxy"
save(tt_both, riskLoci, file = 'rdas/twas_exp_DGandHippoMerged_GWASannotated.Rdata') 

####################
## explore #########
####################
tt_both = as.data.frame(tt_both)
tt_DG = tt_both[tt_both$Region == "DG",]
tt_Hippo = tt_both[tt_both$Region == "HIPPO",]

tt_DG_fdr = tt_DG[tt_DG$TWAS.FDR < 0.05,]
tt_DG_bonf = tt_DG[tt_DG$TWAS.Bonf < 0.05,]

## for big table
nrow(tt_DG)
table(tt_DG$feature)
sapply(split(tt_DG, tt_DG$feature), function(x) length(unique(x$geneid)))
length(unique(tt_DG$geneid))

## overall numbers
nrow(tt_DG_fdr)
table(tt_DG_fdr$feature)
length(unique(tt_DG_fdr$geneid))

nrow(tt_DG_bonf)
table(tt_DG_bonf$feature)
length(unique(tt_DG_bonf$geneid))

#####
## cross tab

prop.table(table(tt_DG_bonf$BEST.GWAS.status != "Other"))
table(tt_DG_bonf$BEST.GWAS.status != "Other")
## examples
tt_DG_bonf[tt_DG_bonf$genesymbol %in% c("CACNA1C", "GRM3"),]

table(tt_DG_fdr$BEST.GWAS.status != "Other")
length(unique(tt_DG_fdr$geneid[tt_DG_fdr$BEST.GWAS.status == "Other"]))
prop.table(table(tt_DG_fdr$BEST.GWAS.status != "Other"))

bonfOut = tt_DG_bonf[which(tt_DG_bonf$BEST.GWAS.status == "Other"),]
bonfOut = bonfOut[order(bonfOut$BEST.GWAS.P, decreasing=TRUE),]
bonfOut$genesymbol[!duplicated(bonfOut$geneid)]
fdrOut = tt_DG_fdr[which(tt_DG_fdr$BEST.GWAS.status == "Other"),]
fdrOut = fdrOut[order(fdrOut$BEST.GWAS.P, decreasing=TRUE),]
fdrOut$genesymbol[!duplicated(fdrOut$geneid)]

table(tt_DG_fdr$BEST.GWAS.status != "Other", tt_DG_fdr$feature)
rowSums(table(tt_DG_fdr$BEST.GWAS.status != "Other", tt_DG_fdr$feature))
table(tt_DG_bonf$BEST.GWAS.status != "Other", tt_DG_bonf$feature)
rowSums(table(tt_DG_bonf$BEST.GWAS.status != "Other", tt_DG_bonf$feature))
###################
# hippo compare ###
###################
tt_DG$hippo_TWASp = tt_Hippo$TWAS.P[match(tt_DG$ID, tt_Hippo$ID)]
tt_DG$hippo_TWAS_FDR = tt_Hippo$TWAS.FDR[match(tt_DG$ID, tt_Hippo$ID)] < 0.05
tt_DG$hippo_TWAS_Bonf = tt_Hippo$TWAS.Bonf[match(tt_DG$ID, tt_Hippo$ID)] < 0.05

tt_DG_fdr$hippo_TWASp = tt_Hippo$TWAS.P[match(tt_DG_fdr$ID, tt_Hippo$ID)]
tt_DG_fdr$hippo_TWAS_FDR = tt_Hippo$TWAS.FDR[match(tt_DG_fdr$ID, tt_Hippo$ID)] < 0.05
tt_DG_fdr$hippo_TWAS_Bonf = tt_Hippo$TWAS.Bonf[match(tt_DG_fdr$ID, tt_Hippo$ID)] < 0.05

tt_DG_bonf$hippo_TWASp = tt_Hippo$TWAS.P[match(tt_DG_bonf$ID, tt_Hippo$ID)]
tt_DG_bonf$hippo_TWAS_FDR = tt_Hippo$TWAS.FDR[match(tt_DG_bonf$ID, tt_Hippo$ID)] < 0.05
tt_DG_bonf$hippo_TWAS_Bonf = tt_Hippo$TWAS.Bonf[match(tt_DG_bonf$ID, tt_Hippo$ID)] < 0.05

table(!is.na(tt_DG$hippo_TWASp))
prop.table(table(!is.na(tt_DG$hippo_TWASp)))
table(!is.na(tt_DG$hippo_TWASp), tt_DG$feature)
rowSums(table(!is.na(tt_DG$hippo_TWASp), tt_DG$feature))
table(tt_DG_fdr$hippo_TWASp < 0.05, tt_DG_fdr$feature)
table(tt_DG_fdr$hippo_TWAS_FDR, tt_DG_fdr$feature)

table(is.na(tt_DG_fdr$hippo_TWASp))
table(tt_DG_fdr$hippo_TWASp < 0.05, useNA="ifany")
table(tt_DG_fdr$hippo_TWAS_FDR, useNA="ifany")

table(tt_DG_bonf$hippo_TWASp < 0.05, tt_DG_bonf$feature, useNA="ifany")
table(tt_DG_bonf$hippo_TWAS_FDR, tt_DG_bonf$feature, useNA="ifany")
table(tt_DG_bonf$hippo_TWAS_FDR, useNA="ifany")


### extract
tt = as.data.frame(twas_exp$all)
tt = tt[which(tt$type == "psycm"),]
# > table(tt$feature)
 # exon  gene   jxn    tx
# 37115  5242 16522  9320

tt$FDR = p.adjust(tt$TWAS.P,"fdr")
sum(tt$FDR < 0.05, na.rm=TRUE) # 3009
# [1] 3099
ttSig = tt[which(tt$FDR < 0.05),]
ttSig = ttSig[order(ttSig$TWAS.P),]

## number of unique genes
length(unique(ttSig$geneid)) #  1072

## quick plots
vennDiagram(vennCounts(table(ttSig$geneid, ttSig$feature) > 0))

## raggr output?
load('eQTL_GWAS_riskSNPs_n596/eqtl_tables/mergedEqtl_output_dg_raggr_4features.rda')
sigEqtl = allEqtl[allEqtl$FDR < 0.01,]

ttSig$in_eqtl = ttSig$ID %in% sigEqtl$gene
mean(ttSig$in_eqtl) # 0.13
sigEqtl$in_twas_sig = sigEqtl$gene %in% ttSig$ID
sigEqtl$in_twas_tested = sigEqtl$gene %in% tt$ID
mean(sigEqtl$in_twas) # 0.54
mean(sigEqtl$in_twas_tested) # 0.935
sapply(10^seq(-4,-12,by=-2), function(x) 
	mean(sigEqtl$in_twas_sig[sigEqtl$pvalue < x])/
		mean(sigEqtl$in_twas_tested[sigEqtl$pvalue < x]))

## included?
ttInc = as.data.frame(twas_exp$included)
ttInc$FDR = tt$FDR[match(ttInc$ID, tt$ID)]
vennDiagram(vennCounts(table(ttInc$geneid, ttInc$feature) > 0))

library('data.table')
psycm_ori <- fread('/dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/psycm/clozuk_pgc2.meta.sumstats.txt')
psycm <- fread('/dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/psycm/clozuk_pgc2.meta.reformatted.sumstats')

summary(twas_exp$included$TWAS.P)
summary(twas_exp$all$TWAS.P)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
#    0.00    0.09    0.33    0.39    0.66    1.00  147856

## load eQTL raggr results
raggr_files <- c(
    'dg_raggr' = 'eQTL_GWAS_riskSNPs_n596/eqtl_tables/mergedEqtl_output_dg_raggr_4features.rda',
    'hippo_raggr' = 'eQTL_GWAS_riskSNPs_n596/eqtl_tables/mergedEqtl_output_hippo_raggr_4features.rda')
raggr <- map(raggr_files, function(x) {
    load(x, verbose = TRUE)
    return(allEqtl)
})
names(raggr) <- c('DG', 'HIPPO')


map_dbl(raggr, ~ sum(.x$FDR < 0.01 ))
   # DG HIPPO
# 39714 34325


# riskLoci <- read.csv("/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_GWAS_riskSNPs/pgc_riskLoci.csv", stringsAsFactors=FALSE)


## Read in the files that Emily cleaned up at
## https://github.com/LieberInstitute/brainseq_phase2/blob/master/eQTL_GWAS_riskSNPs/create_eqtl_table_indexInfo.R
raggr_clean_files <- c(
    'DG' = '/dcl01/ajaffe/data/lab/dg_hippo_paper/eQTL_GWAS_riskSNPs_n596/raggr_179_snps_dg_eqtls_fdr05.csv',
    'HIPPO' = '/dcl01/ajaffe/data/lab/dg_hippo_paper/eQTL_GWAS_riskSNPs_n596/raggr_179_snps_hippo_eqtls_fdr05.csv'
)
raggr_clean <- map(raggr_clean_files, read.csv, stringsAsFactors = FALSE)
names(raggr_clean) <- names(raggr_clean_files)
raggr_clean = map(raggr_clean, function(x) x[x$FDR < 0.01,])

## Start looking into how to match tables by SNP id..
hmm <- twas_exp$all$BEST.GWAS.ID
length(hmm)
# [1] 68199
length(unique(hmm))
# [1] 3395
hmm <- hmm[!is.na(hmm)]
length(hmm)
# [1] 260187
length(unique(hmm))
# [1] 7975

h <- unique(hmm)
table(grepl(':', h))
# FALSE  TRUE
 # 2593   801
head(h[grepl(':', h)])

h_eqtl <- unique(twas_exp$all$EQTL.ID)
h_eqtl <- h_eqtl[!is.na(h_eqtl)]


h_psycm <- unique(twas_exp$all$EQTL.ID[twas_exp$all$type == 'psycm'])
h_psycm <- h_psycm[!is.na(h_psycm)]

table(raggr_clean$HIPPO$SNP %in% h)
# FALSE  TRUE
# 33683   642

table(raggr_clean$DG$SNP %in% h)
# FALSE  TRUE
# 38961   753

table(raggr_clean$HIPPO$SNP %in% h_eqtl)
# FALSE  TRUE
# 33022  1303
table(raggr_clean$DG$SNP %in% h_eqtl)
# FALSE  TRUE
# 37840  1874

table(raggr_clean$HIPPO$SNP %in% h_psycm)
# FALSE  TRUE
# 33022  1303
table(raggr_clean$DG$SNP %in% h_psycm)
# FALSE  TRUE
# 37840  1874



table(unique(raggr_clean$HIPPO$SNP) %in% h)
# FALSE  TRUE
#  5353   157
table(unique(raggr_clean$DLPFC$SNP) %in% h)
# FALSE  TRUE
#  6600   180

table(unique(raggr_clean$HIPPO$SNP) %in% h_eqtl)
# FALSE  TRUE
#  5302   208
table(unique(raggr_clean$DLPFC$SNP) %in% h_eqtl)
# FALSE  TRUE
#  6555   225

table(unique(raggr_clean$HIPPO$SNP) %in% h_psycm)
# FALSE  TRUE
#  5305   205
table(unique(raggr_clean$DLPFC$SNP) %in% h_psycm)
# FALSE  TRUE
#  6558   222


table(unique(raggr_clean$HIPPO$IndexSNP) %in% h)
# FALSE  TRUE
#   100     3
table(unique(raggr_clean$DLPFC$IndexSNP) %in% h)
# FALSE  TRUE
#   112     4

table(unique(raggr_clean$HIPPO$IndexSNP) %in% h_eqtl)
# FALSE
  # 103
table(unique(raggr_clean$DLPFC$IndexSNP) %in% h_eqtl)
# FALSE
#   116

table(unique(raggr_clean$HIPPO$IndexSNP) %in% h_psycm)
# FALSE
#   103
table(unique(raggr_clean$DLPFC$IndexSNP) %in% h_psycm)
# FALSE
#   116

## Read in the BIM files used for computing the weights (one per region)
## https://github.com/LieberInstitute/brainseq_phase2/blob/master/twas/filter_data/filter_snps.R#L128-L174
library('data.table')
bims <- map(c('DLPFC', 'HIPPO'), function(region) {
    bimfile <- paste0(
       '/dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/filter_data/LIBD_Brain_Illumina_h650_1M_Omni5M_Omni2pt5_Macrogen_imputed_run2_LDfiltered_',
       region,
       '.bim'
   )
   fread(bimfile,
       col.names = c('chr', 'snp', 'position', 'basepair', 'allele1', 'allele2'),
       colClasses = c('character', 'character', 'numeric', 'integer', 'character', 'character')
   )
})
names(bims) <- c('DLPFC', 'HIPPO')


## Check matches by SNP name
table(unique(raggr_clean$HIPPO$SNP) %in% bims$HIPPO$snp)
# FALSE  TRUE
#  4483  1027
table(unique(raggr_clean$DLPFC$SNP) %in% bims$DLPFC$snp)
# FALSE  TRUE
#  5586  1194


table(unique(raggr_clean$HIPPO$IndexSNP) %in% bims$HIPPO$snp)
# FALSE  TRUE
#   100     3
table(unique(raggr_clean$DLPFC$IndexSNP) %in% bims$DLPFC$snp)
# FALSE  TRUE
#   112     4

## Check with pairs of chr:position
## instead of SNP names
make_pairs <- function(rag, bim) {
    pairs_rag <- paste(
        gsub('chr', '', rag$chr_hg38),
        ':',
        rag$pos_hg38,
        sep = ''
    )
    pairs_bim <- paste(
        bim$chr,
        ':',
        bim$basepair,
        sep = ''
    )
    return(list(rag = pairs_rag, bim = pairs_bim))
}
pairs <- map2(raggr_clean, bims, make_pairs)

## Values of "trues" match the searches by SNP name...
with(pairs$HIPPO, table(unique(rag) %in% unique(bim)))
# FALSE  TRUE
#  4470  1027
with(pairs$DLPFC, table(unique(rag) %in% unique(bim)))
# FALSE  TRUE
#  5564  1194


## Read in the original BSP2 bim file with hg19 coordinates
bfile <- '/dcl01/lieber/ajaffe/Brain/Imputation/Merged/LIBD_Brain_Illumina_h650_1M_Omni5M_Omni2pt5_Macrogen_imputed_run2.bim'
bsp2_bim <- fread(
    bfile,
    col.names = c('chr', 'snp', 'position', 'basepair', 'allele1', 'allele2')
)

## Check by SNP name...
table(unique(raggr_clean$HIPPO$SNP) %in% bsp2_bim$snp)
# TRUE
# 5510
table(unique(raggr_clean$DLPFC$SNP) %in% bsp2_bim$snp)
# TRUE
# 6780

## Ok, not all index snps were in our data to begin with
## which is why we did the raggr analysis in the
## first place
table(unique(raggr_clean$HIPPO$IndexSNP) %in% bsp2_bim$snp)
# FALSE  TRUE
#    30    73
table(unique(raggr_clean$DLPFC$IndexSNP) %in% bsp2_bim$snp)
# FALSE  TRUE
#    31    85


table(bims$HIPPO$snp %in% bsp2_bim$snp)
#    TRUE
# 1022527
table(bims$DLPFC$snp %in% bsp2_bim$snp)
#    TRUE
# 1022527



check_by_locus <- function(rag, ref) {
    by_loc <- split(rag$SNP, rag$IndexSNP)
    map_dbl(by_loc, ~ sum(.x %in% ref))
}

by_locus <- map(raggr_clean, check_by_locus, ref = h)
map(by_locus, ~ table(.x > 0))
# $HIPPO
#
# FALSE  TRUE
#    40    63
#
# $DLPFC
#
# FALSE  TRUE
#    48    68

by_locus_considered <- map2(
    raggr_clean,
    list(HIPPO = bims$HIPPO$snp, DLPFC = bims$DLPFC$snp),
    check_by_locus
)
map(by_locus_considered, ~ table(.x > 0))
# $HIPPO
#
# FALSE  TRUE
#    16    87
#
# $DLPFC
#
# FALSE  TRUE
#    17    99

map2(
    by_locus,
    by_locus_considered,
    ~ addmargins(table(
        'among TWAS-all best' = .x > 0,
        'among SNPs for TWAS weights' = .y > 0
    ))
)
# $HIPPO
#                    among SNPs for TWAS weights
# among TWAS-all best FALSE TRUE Sum
#               FALSE    16   24  40
#               TRUE      0   63  63
#               Sum      16   87 103
#
# $DLPFC
#                    among SNPs for TWAS weights
# among TWAS-all best FALSE TRUE Sum
#               FALSE    17   31  48
#               TRUE      0   68  68
#               Sum      17   99 116


by_locus_eqtl <- map(raggr_clean, check_by_locus, ref = h_eqtl)
map(by_locus_eqtl, ~ table(.x > 0))
# $HIPPO
#
# FALSE  TRUE
#    57    46
#
# $DLPFC
#
# FALSE  TRUE
#    67    49

map2(
    by_locus_eqtl,
    by_locus_considered,
    ~ addmargins(table(
        'among TWAS-all best' = .x > 0,
        'among SNPs for TWAS weights' = .y > 0
    ))
)
# $HIPPO
#                    among SNPs for TWAS weights
# among TWAS-all best FALSE TRUE Sum
#               FALSE    16   41  57
#               TRUE      0   46  46
#               Sum      16   87 103
#
# $DLPFC
#                    among SNPs for TWAS weights
# among TWAS-all best FALSE TRUE Sum
#               FALSE    17   50  67
#               TRUE      0   49  49
#               Sum      17   99 116


by_locus_psycm <- map(raggr_clean, check_by_locus, ref = h_psycm)
map(by_locus_psycm, ~ table(.x > 0))
# $HIPPO
#
# FALSE  TRUE
#    58    45
#
# $DLPFC
#
# FALSE  TRUE
#    68    48
map2(
    by_locus_psycm,
    by_locus_considered,
    ~ addmargins(table(
        'among TWAS-all best' = .x > 0,
        'among SNPs for TWAS weights' = .y > 0
    ))
)
# $HIPPO
#                    among SNPs for TWAS weights
# among TWAS-all best FALSE TRUE Sum
#               FALSE    16   42  58
#               TRUE      0   45  45
#               Sum      16   87 103
#
# $DLPFC
#                    among SNPs for TWAS weights
# among TWAS-all best FALSE TRUE Sum
#               FALSE    17   51  68
#               TRUE      0   48  48
#               Sum      17   99 116


check_with_p <- function(threshold = 0.05) {
    h_psycm_p5 <- unique(twas_exp$all$EQTL.ID[twas_exp$all$type == 'psycm' & twas_exp$all$TWAS.P < threshold])
    h_psycm_p5 <- h_psycm_p5[!is.na(h_psycm_p5)]
    by_locus_psycm_p5 <- map(raggr_clean, check_by_locus, ref = h_psycm_p5)
    map2(
        by_locus_psycm_p5,
        by_locus_considered,
        ~ addmargins(table(
            'among TWAS-all best' = .x > 0,
            'among SNPs for TWAS weights' = .y > 0
        ))
    )
}
check_with_p()
# $HIPPO
#                    among SNPs for TWAS weights
# among TWAS-all best FALSE TRUE Sum
#               FALSE    16   43  59
#               TRUE      0   44  44
#               Sum      16   87 103
#
# $DLPFC
#                    among SNPs for TWAS weights
# among TWAS-all best FALSE TRUE Sum
#               FALSE    17   53  70
#               TRUE      0   46  46
#               Sum      17   99 116


check_with_p(1e-4)
# $HIPPO
#                    among SNPs for TWAS weights
# among TWAS-all best FALSE TRUE Sum
#               FALSE    16   45  61
#               TRUE      0   42  42
#               Sum      16   87 103
#
# $DLPFC
#                    among SNPs for TWAS weights
# among TWAS-all best FALSE TRUE Sum
#               FALSE    17   55  72
#               TRUE      0   44  44
#               Sum      17   99 116

check_with_p(1e-6)
# $HIPPO
#                    among SNPs for TWAS weights
# among TWAS-all best FALSE TRUE Sum
#               FALSE    16   50  66
#               TRUE      0   37  37
#               Sum      16   87 103
#
# $DLPFC
#                    among SNPs for TWAS weights
# among TWAS-all best FALSE TRUE Sum
#               FALSE    17   61  78
#               TRUE      0   38  38
#               Sum      17   99 116

check_with_p(1e-8)
# $HIPPO
#                    among SNPs for TWAS weights
# among TWAS-all best FALSE TRUE Sum
#               FALSE    16   62  78
#               TRUE      0   25  25
#               Sum      16   87 103
#
# $DLPFC
#                    among SNPs for TWAS weights
# among TWAS-all best FALSE TRUE Sum
#               FALSE    17   73  90
#               TRUE      0   26  26
#               Sum      17   99 116


gene_by_locus <- function(rag, ref) {
    by_loc <- split(rag$gene, rag$IndexSNP)
    map_dbl(by_loc, ~ sum(.x %in% ref))
}


g_by_locus <- map(names(rse), function(feature) {
    map(
        map(raggr_clean, ~ subset(.x, tolower(Type) == feature)),
        gene_by_locus,
        ref = twas_exp$all$ID[twas_exp$all$feature == feature]
    )
})
names(g_by_locus) <- names(rse)
map(g_by_locus, function(x) {
    r <- map_dfr(x, ~ table(.x > 0))
    r$state <- c(FALSE, TRUE)
    return(r)
})
# $gene
# # A tibble: 2 x 3
#   HIPPO DLPFC state
#   <int> <int> <lgl>
# 1     5    11 FALSE
# 2    45    56 TRUE
#
# $exon
# # A tibble: 2 x 3
#   HIPPO DLPFC state
#   <int> <int> <lgl>
# 1     9    12 FALSE
# 2    68    75 TRUE
#
# $jxn
# # A tibble: 2 x 3
#   HIPPO DLPFC state
#   <int> <int> <lgl>
# 1    12    19 FALSE
# 2    76    79 TRUE
#
# $tx
# # A tibble: 2 x 3
#   HIPPO DLPFC state
#   <int> <int> <lgl>
# 1     9     9 FALSE
# 2    56    67 TRUE


### Ehem.... matching gene ids vs snp ids... ehem... T_T
# g_by_locus_considered <- map(names(rse), function(feature) {
#     map2(
#         map(raggr_clean, ~ subset(.x, tolower(Type) == feature)),
#         list(HIPPO = bims$HIPPO$snp, DLPFC = bims$DLPFC$snp),
#         gene_by_locus
#     )
# })
# names(g_by_locus_considered) <- names(rse)
# map(g_by_locus_considered, function(x) {
#     r <- map_dfr(x, ~ table(factor(.x > 0, levels = c('FALSE', 'TRUE'))))
#     r$state <- c(FALSE, TRUE)
#     return(r)
# })
#

clean_tabs <- function(l) {
    map2_dfr(l, names(l), function(x, y) {
        x$feature <- y
        return(x)
    })
}

g_by_locus_reg <- map(names(rse), function(feature) {
    map2(
        map(raggr_clean, ~ subset(.x, tolower(Type) == feature)),
        map(names(raggr_clean), ~ twas_exp$all$ID[twas_exp$all$feature == feature & twas_exp$all$region == .x]),
        gene_by_locus
    )
})
names(g_by_locus_reg) <- names(rse)
clean_tabs(map(g_by_locus_reg, function(x) {
    r <- map_dfr(x, ~ table(.x > 0))
    r$state <- c(FALSE, TRUE)
    return(r)
}))
# # A tibble: 8 x 4
#   HIPPO DLPFC state feature
#   <int> <int> <lgl> <chr>
# 1     6    11 FALSE gene
# 2    44    56 TRUE  gene
# 3    11    13 FALSE exon
# 4    66    74 TRUE  exon
# 5    13    19 FALSE jxn
# 6    75    79 TRUE  jxn
# 7    12     9 FALSE tx
# 8    53    67 TRUE  tx

clean_tabs_type <- function(l) {
    map2_dfr(l, names(l), function(x, y) {
        x$type <- y
        return(x)
    })
}


g_by_type <- map(c('psycm', 'pgc2'), function(type) {
    g_reg <- map(names(rse), function(feature) {
        map2(
            map(raggr_clean, ~ subset(.x, tolower(Type) == feature)),
            map(names(raggr_clean), ~ twas_exp$all$ID[twas_exp$all$feature == feature & twas_exp$all$region == .x & twas_exp$all$type == type]),
            gene_by_locus
        )
    })
    names(g_reg) <- names(rse)
    map(g_reg, function(x) {
        r <- map_dfr(x, ~ table(.x > 0))
        r$state <- c(FALSE, TRUE)
        return(r)
    })
})
names(g_by_type) <- c('psycm', 'pgc2')
clean_tabs_type(map(g_by_type, clean_tabs))
# # A tibble: 16 x 5
#    HIPPO DLPFC state feature type
#    <int> <int> <lgl> <chr>   <chr>
#  1     6    11 FALSE gene    psycm
#  2    44    56 TRUE  gene    psycm
#  3    17    13 FALSE exon    psycm
#  4    60    74 TRUE  exon    psycm
#  5    13    19 FALSE jxn     psycm
#  6    75    79 TRUE  jxn     psycm
#  7    12     9 FALSE tx      psycm
#  8    53    67 TRUE  tx      psycm
#  9     6    11 FALSE gene    pgc2
# 10    44    56 TRUE  gene    pgc2
# 11    11    13 FALSE exon    pgc2
# 12    66    74 TRUE  exon    pgc2
# 13    13    19 FALSE jxn     pgc2
# 14    75    79 TRUE  jxn     pgc2
# 15    12     9 FALSE tx      pgc2
# 16    53    67 TRUE  tx      pgc2


g_by_type_inc <- map(c('psycm', 'pgc2'), function(type) {
    g_reg <- map(names(rse), function(feature) {
        map2(
            map(raggr_clean, ~ subset(.x, tolower(Type) == feature)),
            map(names(raggr_clean), ~ twas_exp$included$ID[twas_exp$included$feature == feature & twas_exp$included$region == .x & twas_exp$included$type == type]),
            gene_by_locus
        )
    })
    names(g_reg) <- names(rse)
    map(g_reg, function(x) {
        r <- map_dfr(x, ~ table(.x > 0))
        r$state <- c(FALSE, TRUE)
        return(r)
    })
})
names(g_by_type_inc) <- c('psycm', 'pgc2')
clean_tabs_type(map(g_by_type_inc, clean_tabs))
# # A tibble: 16 x 5
#    HIPPO DLPFC state feature type
#    <int> <int> <lgl> <chr>   <chr>
#  1    30    38 FALSE gene    psycm
#  2    20    29 TRUE  gene    psycm
#  3    46    45 FALSE exon    psycm
#  4    31    42 TRUE  exon    psycm
#  5    48    48 FALSE jxn     psycm
#  6    40    50 TRUE  jxn     psycm
#  7    39    41 FALSE tx      psycm
#  8    26    35 TRUE  tx      psycm
#  9    32    42 FALSE gene    pgc2
# 10    18    25 TRUE  gene    pgc2
# 11    51    58 FALSE exon    pgc2
# 12    26    29 TRUE  exon    pgc2
# 13    61    57 FALSE jxn     pgc2
# 14    27    41 TRUE  jxn     pgc2
# 15    43    44 FALSE tx      pgc2
# 16    22    32 TRUE  tx      pgc2




check_by_feature <- function(tw, cut) {
    res <- map(c('psycm', 'pgc2'), function(type) {
        g_reg <- map(names(rse), function(feature) {
            map2(
                map(raggr_clean, ~ subset(.x, tolower(Type) == feature)),
                map(names(raggr_clean), ~ tw$ID[tw$feature == feature & tw$region == .x & tw$type == type & tw$TWAS.P < cut]),
                ~ match(unique(.x$gene), .y)
            )
        })
        names(g_reg) <- names(rse)
        map(g_reg, function(x) {
            # r <- map_dfr(x, ~ table(factor(!is.na(.x), levels = c('FALSE', 'TRUE'))))
            r <- map_dfr(x, ~ table(!is.na(.x)))
            r$state <- c(FALSE, TRUE)
            return(r)
        })
    })
    names(res) <- c('psycm', 'pgc2')
    
    clean_tabs_type(map(res, clean_tabs))
}

check_by_feature(twas_exp$all, cut = 1.01)
# # A tibble: 16 x 5
#    HIPPO DLPFC state feature type
#    <int> <int> <lgl> <chr>   <chr>
#  1    53    66 FALSE gene    psycm
#  2    70   105 TRUE  gene    psycm
#  3   408   669 FALSE exon    psycm
#  4   449   694 TRUE  exon    psycm
#  5   232   295 FALSE jxn     psycm
#  6   275   364 TRUE  jxn     psycm
#  7   104   161 FALSE tx      psycm
#  8   140   171 TRUE  tx      psycm
#  9    53    66 FALSE gene    pgc2
# 10    70   105 TRUE  gene    pgc2
# 11   395   669 FALSE exon    pgc2
# 12   462   694 TRUE  exon    pgc2
# 13   232   295 FALSE jxn     pgc2
# 14   275   364 TRUE  jxn     pgc2
# 15   104   161 FALSE tx      pgc2
# 16   140   171 TRUE  tx      pgc2

check_by_feature(twas_exp$all, cut = 0.01)
# # A tibble: 16 x 5
#    HIPPO DLPFC state feature type
#    <int> <int> <lgl> <chr>   <chr>
#  1    75    98 FALSE gene    psycm
#  2    48    73 TRUE  gene    psycm
#  3   511   841 FALSE exon    psycm
#  4   346   522 TRUE  exon    psycm
#  5   307   386 FALSE jxn     psycm
#  6   200   273 TRUE  jxn     psycm
#  7   140   210 FALSE tx      psycm
#  8   104   122 TRUE  tx      psycm
#  9    83   106 FALSE gene    pgc2
# 10    40    65 TRUE  gene    pgc2
# 11   520   907 FALSE exon    pgc2
# 12   337   456 TRUE  exon    pgc2
# 13   331   419 FALSE jxn     pgc2
# 14   176   240 TRUE  jxn     pgc2
# 15   148   225 FALSE tx      pgc2
# 16    96   107 TRUE  tx      pgc2

check_by_feature(twas_exp$included, cut = 0.01)
# # A tibble: 16 x 5
#    HIPPO DLPFC state feature type
#    <int> <int> <lgl> <chr>   <chr>
#  1   105   141 FALSE gene    psycm
#  2    18    30 TRUE  gene    psycm
#  3   827  1321 FALSE exon    psycm
#  4    30    42 TRUE  exon    psycm
#  5   467   610 FALSE jxn     psycm
#  6    40    49 TRUE  jxn     psycm
#  7   220   300 FALSE tx      psycm
#  8    24    32 TRUE  tx      psycm
#  9   107   145 FALSE gene    pgc2
# 10    16    26 TRUE  gene    pgc2
# 11   833  1335 FALSE exon    pgc2
# 12    24    28 TRUE  exon    pgc2
# 13   480   619 FALSE jxn     pgc2
# 14    27    40 TRUE  jxn     pgc2
# 15   224   302 FALSE tx      pgc2
# 16    20    30 TRUE  tx      pgc2



check_by_feature_rev <- function(tw, cut) {
    res <- map(c('psycm', 'pgc2'), function(type) {
        g_reg <- map(names(rse), function(feature) {
            map2(
                map(raggr_clean, ~ subset(.x, tolower(Type) == feature)),
                map(names(raggr_clean), ~ tw$ID[tw$feature == feature & tw$region == .x & tw$type == type & tw$TWAS.P < cut]),
                ~ match(unique(.y[!is.na(.y)]), .x$gene)
            )
        })
        names(g_reg) <- names(rse)
        map(g_reg, function(x) {
            # r <- map_dfr(x, ~ table(factor(!is.na(.x), levels = c('FALSE', 'TRUE'))))
            r <- map_dfr(x, ~ table(!is.na(.x)))
            r$state <- c(FALSE, TRUE)
            return(r)
        })
    })
    names(res) <- c('psycm', 'pgc2')
    
    clean_tabs_type(map(res, clean_tabs))
}

check_by_feature_rev(twas_exp$all, cut = 1.01)
# # A tibble: 16 x 5
#    HIPPO DLPFC state feature type
#    <int> <int> <lgl> <chr>   <chr>
#  1  3907  5377 FALSE gene    psycm
#  2    70   105 TRUE  gene    psycm
#  3 25237 38437 FALSE exon    psycm
#  4   449   694 TRUE  exon    psycm
#  5 15832 20289 FALSE jxn     psycm
#  6   275   364 TRUE  jxn     psycm
#  7  7014  8890 FALSE tx      psycm
#  8   140   171 TRUE  tx      psycm
#  9  3960  5450 FALSE gene    pgc2
# 10    70   105 TRUE  gene    pgc2
# 11 29404 39031 FALSE exon    pgc2
# 12   462   694 TRUE  exon    pgc2
# 13 16087 20556 FALSE jxn     pgc2
# 14   275   364 TRUE  jxn     pgc2
# 15  7132  9035 FALSE tx      pgc2
# 16   140   171 TRUE  tx      pgc2

check_by_feature_rev(twas_exp$all, cut = 0.01)
# # A tibble: 16 x 5
#    HIPPO DLPFC state feature type
#    <int> <int> <lgl> <chr>   <chr>
#  1   354   514 FALSE gene    psycm
#  2    48    73 TRUE  gene    psycm
#  3  2389  3689 FALSE exon    psycm
#  4   346   522 TRUE  exon    psycm
#  5  1533  1913 FALSE jxn     psycm
#  6   200   273 TRUE  jxn     psycm
#  7   720   897 FALSE tx      psycm
#  8   104   122 TRUE  tx      psycm
#  9   264   391 FALSE gene    pgc2
# 10    40    65 TRUE  gene    pgc2
# 11  2142  2904 FALSE exon    pgc2
# 12   337   456 TRUE  exon    pgc2
# 13  1235  1519 FALSE jxn     pgc2
# 14   176   240 TRUE  jxn     pgc2
# 15   608   691 FALSE tx      pgc2
# 16    96   107 TRUE  tx      pgc2

check_by_feature_rev(twas_exp$included, cut = 0.01)
# # A tibble: 16 x 5
#    HIPPO DLPFC state feature type
#    <int> <int> <lgl> <chr>   <chr>
#  1    60    96 FALSE gene    psycm
#  2    18    30 TRUE  gene    psycm
#  3    62    74 FALSE exon    psycm
#  4    30    42 TRUE  exon    psycm
#  5    95    93 FALSE jxn     psycm
#  6    40    49 TRUE  jxn     psycm
#  7    89    98 FALSE tx      psycm
#  8    24    32 TRUE  tx      psycm
#  9    40    55 FALSE gene    pgc2
# 10    16    26 TRUE  gene    pgc2
# 11    45    57 FALSE exon    pgc2
# 12    24    28 TRUE  exon    pgc2
# 13    67    53 FALSE jxn     pgc2
# 14    27    40 TRUE  jxn     pgc2
# 15    58    54 FALSE tx      pgc2
# 16    20    30 TRUE  tx      pgc2








## put this on pause...


library('UpSetR')

tw <- twas_exp$included

library('VennDiagram')
library('RColorBrewer')
venn_cols <- brewer.pal('Set1', n = 4)

genes <- with(subset(tw, type == 'psycm'), split(geneid, feature))
genes <- genes[names(rse)]
genes <- map(genes, ~ .x[!is.na(.x)])

names(venn_cols) <- names(genes)

make_venn <- function(genes, title = 'DE features grouped by gene id') {
    v <- venn.diagram(genes, filename = NULL,
        main = title,
        col = 'transparent', fill = venn_cols[names(genes)],
        alpha = 0.5, margin = 0,
        main.cex = 2, cex = 2, cat.fontcase = 'bold', cat.cex = 2,
        cat.col = venn_cols[names(genes)])
    grid.newpage()
    grid.draw(v)
}

dir.create('pdf', showWarnings = FALSE)

pdf('pdf/test.pdf', useDingbats = FALSE)
make_venn(genes)
dev.off()
system('rm VennDiagram*')

library('gplots')
genes2 <- map(genes, unique)
names(genes2) <- names(genes)

x <- venn(genes2, show.plot = FALSE)

## Get the matrix out
y <- matrix(x, ncol = ncol(x), dimnames = attr(x, 'dimnames'))
z <- map_dbl(names(genes), ~ sum(y[ y[, .x] > 0, 'num']))
names(z) <- names(genes)
z



tw_up <- map_dfr(split(tw, tw$geneid), function(i) {
    
    data.frame(
        geneid = unique(i$geneid),
        genesymbol = unique(i$genesymbol),
        gene = ifelse('gene' %in% i$feature, 1, 0),
        exon = ifelse('exon' %in% i$feature, 1, 0),
        jxn = ifelse('jxn' %in% i$feature, 1, 0),
        tx = ifelse('tx' %in% i$feature, 1, 0),
        DLPFC = ifelse('DLPFC' %in% i$region, 1, 0),
        HIPPO = ifelse('HIPPO' %in% i$region, 1, 0),
        psycm = ifelse('psycm' %in% i$type, 1, 0),
        pgc2 = ifelse('pgc2' %in% i$type, 1, 0),
        stringsAsFactors = FALSE
    )
    
})

set_meta <- data.frame(
    sets = c(names(rse), 'DLPFC', 'HIPPO', 'psycm', 'pgc2'),
    colors = c(brewer.pal('Set1', n = 4), 'dark orange', 'skyblue3', 'sienna1', 'springgreen4'),
    stringsAsFactors = FALSE
)
cols <- set_meta$colors
names(cols) <- set_meta$colors


pdf('pdf/test.pdf', useDingbats = FALSE, width = 14)
upset(tw_up,
    nset = 8,
    sets = c(names(rse), 'DLPFC', 'HIPPO', 'psycm', 'pgc2'),
    # empty.intersections = 'on',
    group.by = 'sets',
    # nintersects = NA,
    nintersects = 70,
    keep.order = TRUE,
    set.metadata = list(
        data = set_meta,
        plots = list(
            list(
                type = 'matrix_rows',
                column = 'colors',
                colors = cols
            )
        )
    )
)
dev.off()

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()