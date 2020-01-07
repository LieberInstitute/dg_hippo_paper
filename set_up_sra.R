###
library(SummarizedExperiment)
library(parallel)

## load data
load("count_data/astellas_dg_hg38_rseGene_n263.rda")
pd = colData(rse_gene)

####################
## add read info ###
####################
new_fq_path = "/dcl02/lieber/ajaffe/DGGCL_RNAseq_Reads"
dir.create(new_fq_path)

pd$leftReads = paste0("/dcl01/ajaffe/data/lab/dg_hippo/preprocessed_data/paired_end_n292/merged_fastq/",
	pd$SAMPLE_ID, ".fastq.gz")
pd$rightReads = paste0("/dcl01/ajaffe/data/lab/dg_hippo/preprocessed_data/paired_end_n292/merged_fastq/",
	pd$SAMPLE_ID, "_read2.fastq.gz")
all(file.exists(c(pd$leftReads,pd$rightReads)))

newLeft = paste0(new_fq_path, "/", pd$RNum, "_R1.fastq.gz")
newRight= paste0(new_fq_path, "/", pd$RNum, "_R2.fastq.gz")

mclapply(1:nrow(pd), function(i) {
	file.copy(pd$leftReads[i], newLeft[i])
	file.copy(pd$rightReads[i], newRight[i])
}, mc.cores=12)

## make biosample attributes
biosample_data = data.frame(
	sample_name = pd$RNum,
	sample_title = paste("DG-GCL RNA Sample", pd$RNum),
	bioproject_accession = "",
	organism = "Homo Sapiens",
	isolate = pd$BrNum,
	age = pd$Age,
	biomaterial_provider = "Lieber Institute for Brain Development; 855 N Wolfe St, Ste 300; Baltimore, MD 21205",
	sex = pd$Sex,
	tissue = "Human Hippocampus",
	cell_line = "",
	cell_subtype = "",
	cell_type = "Granule Cell of Human Dentate Gyrus",
	culture_collection = "",
	dev_stage = "Adult",
	disease = pd$Dx,
	disease_stage = "",
	ethnicity = pd$Race,
	health_state = "Deceased",
	karyotype = "Normal",
	phenotype = "",
	population = "",
	race = pd$Race,
	sample_type = "Laser Capture Microdissection followed by RNA-seq",
	treatment = "",
	description = "",
		stringsAsFactors=FALSE)

write.table(biosample_data, file = "tables/jaffe_Human_BioSample_DG-GCL.tsv",
		sep = "\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

#### make read metadata
sra_metadata = data.frame(
	sample_name = pd$RNum,
	library_ID = pd$RNum, 
	title = paste("DG-GCL RNA Sample", pd$RNum),
	library_strategy = "RNA-Seq",
	library_source = "TRANSCRIPTOMIC",
	library_selection = "Inverse rRNA",
	library_layout = "paired",
	platform = "ILLUMINA",
	instrument_model = "Illumina HiSeq 2000",
	design_description = "The granule cell layer of the human hippocampus was obtained with laser capture microdissection, RNA was extracted and subjected to RNA sequencing after creating libraries with the Illumina Ribo-Zero protocol, and sequenced on an Illumina HiSeq 2000 with 2x100bp reads." ,
	filetype = "fastq",
	filename = paste0(pd$RNum, "_R1.fastq.gz"),
	filename2 = paste0(pd$RNum, "_R2.fastq.gz"),
		stringsAsFactors = FALSE)
write.table(sra_metadata, file = "tables/jaffe_Human_SraMetadata_DG-GCL.tsv",
		sep = "\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

## run ascp
thecall = "ascp -i /users/ajaffe/.aspera/aspera.openssh -QT -l100m -k1 -d /dcl02/lieber/ajaffe/DGGCL_RNAseq_Reads/ subasp@upload.ncbi.nlm.nih.gov:uploads/andrew.jaffe@libd.org_tqi4SAFC/DG_GCL"
system(thecall)

## and aws s3 for count objects
thecall = paste("aws s3 cp", list.files("count_data",full=TRUE), "s3://jaffe-nat-neuro-dggcl/")
sapply(thecall, system)
