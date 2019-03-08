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

# #######################
# ### add  sz control? ##
# load("rdas/jxLevel_ageAndSzInteraction.rda")
# jxSzStats = jxSzStats[,c("t_SZ_DG", "P.Value_SZ_DG")]
# load("rdas/exonLevel_ageAndSzInteraction.rda")
# exonSzStats = exonSzStats[,c("t_SZ_DG", "P.Value_SZ_DG")]
# load("rdas/geneLevel_ageAndSzInteraction.rda")
# geneSzStats = geneSzStats[,c("t_SZ_DG", "P.Value_SZ_DG")]
# load("rdas/txLevel_ageAndSzInteraction.rda")
# txSzStats=txSzStats[,c("t_SZ_DG", "P.Value_SZ_DG")]
# szStats = c(geneSzStats, exonSzStats, jxSzStats, txSzStats)

# tt_DG$SZ_t = szStats$t_SZ_DG[match(tt_DG$ID, names(szStats))]

# ## anything?
# plot(tt_DG$SZ_t, tt_DG$TWAS.Z)
# plot(tt_DG$SZ_t[tt_DG$TWAS.FDR<0.05], tt_DG$TWAS.Z[tt_DG$TWAS.FDR<0.05])
# cor(tt_DG$SZ_t[tt_DG$TWAS.FDR<0.05], tt_DG$TWAS.Z[tt_DG$TWAS.FDR<0.05],use="comp")
# plot(tt_DG$SZ_t[tt_DG$TWAS.Bonf<0.05], tt_DG$TWAS.Z[tt_DG$TWAS.Bonf<0.05])
# cor(tt_DG$SZ_t[tt_DG$TWAS.Bonf<0.05], tt_DG$TWAS.Z[tt_DG$TWAS.Bonf<0.05],use="comp")
# # not really

# ## check bonf again 
# bonf = tt_DG[tt_DG$TWAS.Bonf  <0.05,]
# bonf[bonf$genesymbol %in% c("CACNA1C","GRM3"),]

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()