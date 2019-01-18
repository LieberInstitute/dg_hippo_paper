#!/bin/bash
#$ -cwd
#$ -l mem_free=100G,h_vmem=100G,h_fsize=200G
#$ -N inter_eQTL
#$ -o logs/eQTL_inter_4_features.txt
#$ -e logs/eQTL_inter_4_features.txt
#$ -m a
echo "**** Job starts ****"
date

Rscript /dcl01/ajaffe/data/lab/dg_hippo_paper/eQTL_GWAS_riskSNPs_n596/run_eqtls_interaction.R

echo "**** Job ends ****"
date
