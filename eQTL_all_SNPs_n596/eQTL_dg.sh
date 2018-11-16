#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=100G,h_vmem=100G,h_fsize=200G
#$ -N dg_eQTL
#$ -o logs/eQTL_dg_4_features.txt
#$ -e logs/eQTL_dg_4_features.txt
#$ -m a
echo "**** Job starts ****"
date

Rscript /dcl01/ajaffe/data/lab/dg_hippo/eQTL_all_SNPs_n596/run_eqtls_dg.R

echo "**** Job ends ****"
date
