#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=100G,h_vmem=100G,h_fsize=200G
#$ -N hippo_eQTL
#$ -o logs/eQTL_hippo_4_features.txt
#$ -e logs/eQTL_hippo_4_features.txt
#$ -m a
echo "**** Job starts ****"
date

Rscript /dcl01/ajaffe/data/lab/dg_hippo/eQTL_GWAS_riskSNPs_n596/raggr_run_eqtls_hippo.R

echo "**** Job ends ****"
date
