---
layout: page
title: Data
group: navigation
permalink: "data.html"
---
{% include JB/setup %}


Data Availability
------------

### Processed Data

The easiest way to access the underlying data from this paper is to use [this R object](http://LieberInstitute.github.io/as3mt-paper/rdas/rawAndRpkmCounts_plusGenotype_10q24_DLPFC_n738.rda) which contains:

1. Phenotype data: `pd`,  rows are subjects and columns are covariates and identifiers. Identifiers for the expression data are RNA Numbers (`Rxxxx`) and identifiers for the genotype data are Brain Numbers (`Brxxx`). Corresponding text file: [pd](http://LieberInstitute.github.io/as3mt-paper/data/phenotype_n738_LIBD.csv)
2. Gene-level data: `geneRpkm` and `geneCounts` correspond to gene RPKM and count values respectively and have genes/features down the rows and subjects/samples across the columns (in the same order as the rows of the phenotype data). `geneMap` is the corresponding annotation information. Corresponding text files: [geneRpkm](http://LieberInstitute.github.io/as3mt-paper/data/geneRpkm_n738_LIBD.csv), [geneMap](http://LieberInstitute.github.io/as3mt-paper/data/geneAnnotation_Ensembl75.csv)
3. Exon-level data: `exonRpkm` and `exonCounts` correspond to exon RPKM and count values respectively and have exons/features down the rows and subjects/samples across the columns (in the same order as the rows of the phenotype data). `exonMap` is the corresponding annotation information. Corresponding text files: [exonRpkm](http://LieberInstitute.github.io/as3mt-paper/data/exonRpkm_n738_LIBD.csv), [exonMap](http://LieberInstitute.github.io/as3mt-paper/data/exonAnnotation_Ensembl75.csv)
4. Junction-level data: `jRpkm` and `jCounts` correspond to junction RPMs and count values respectively and have exons/features down the rows and subjects/samples across the columns (in the same order as the rows of the phenotype data). `jMap` is the corresponding annotation information, which is a `GRanges` object. Corresponding text files: [jRpm](http://LieberInstitute.github.io/as3mt-paper/data/junctionRpm_n738_LIBD.csv), [jMap](http://LieberInstitute.github.io/as3mt-paper/data/junctionAnnotation_LIBD.csv)
5. Genotype data: `snpMat` contains the number of copies of the minor allele, and have variants down the rows and subjects/samples across the columns (in the same order as the rows of the phenotype data). `snpInfo` contains the corresponding annotation and GWAS statistics information. Corresponding text files: [snpMap](http://LieberInstitute.github.io/as3mt-paper/data/snpMinorCounts_LIBD.csv), [snpInfo](http://LieberInstitute.github.io/as3mt-paper/data/snp_annotation.csv)

### Raw Data

Raw data is available through the Sequence Read Archive (SRA) at accession [SRPxxxx]. 

------------------
<a href="http://libd.org">
<img src="images/LIBD_logo.jpg" alt="Drawing" style="width: 250px;"/>
