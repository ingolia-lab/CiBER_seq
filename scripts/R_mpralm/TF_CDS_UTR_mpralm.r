#!/usr/bin/env Rscript

#You will need to load in three dataframes before you start. 
#First: grna.assign.barcode.grna.good (This dataframe links barcodes with gRNAs)
#second: guide.good.targets (This one links gRNA with target genes)
#third: sgd (this one links gene with description/function and other useful info)

#download sgd, bc_gRNA_assignment table, and gRNA_re-assignment table

if (!file.exists("SGD_features.tab")) {
  sgd <- download.file('https://downloads.yeastgenome.org/curation/chromosomal_feature/SGD_features.tab', destfile="SGD_features.tab")
}
sgd <- read.delim("SGD_features.tab", header=FALSE, quote="",
                  col.names=c("sgdid", "type", "qual", "name", "gene", "alias",
                              "parent", "sgdid2", "chrom", "start", "end",
                              "strand", "genpos", "cver", "sver", "desc"))

options(stringsAsFactors=FALSE)
if (!exists("grna.assign.barcode.grna.good")) {
  grna.assign.barcode.grna.good <- read.delim("~/CiBER_seq_package/all_raw_fasta_gz/PE_bc_gRNA_assignment/grna-assign-barcode-grna-good.txt",
                                stringsAsFactors=FALSE)
}

grna.assign.barcode.grna.good$guide <- as.character(grna.assign.barcode.grna.good$guide)
grna.assign.barcode.grna.good$guide[grna.assign.barcode.grna.good$guide == "No_gRNA"] <- paste("No_gRNA", seq(1:788), sep=" ")

options(stringsAsFactors=FALSE)
if (!exists("guide.good.targets")) {
  guide.good.targets <- read.delim("~/CiBER_seq_package/scripts/guide.good.targets_plus_empty.txt",
                                stringsAsFactors=FALSE)
}

# Load in raw count data and relevant gRNA barcode assignment and gRNA assignment of target files
options(stringsAsFactors=FALSE)
if (!exists("TF_CDS")) {
  TF_CDS <- read.delim("~/CiBER_seq_package/all_raw_fasta_gz/GCN4_CDS_UTR/all_TF_CDS_counts.txt",
                                stringsAsFactors=FALSE)
}

options(stringsAsFactors=FALSE)
if (!exists("TF_UTR")) {
  TF_UTR <- read.delim("~/CiBER_seq_package/all_raw_fasta_gz/GCN4_CDS_UTR/all_TF_UTR_counts.txt",
                                stringsAsFactors=FALSE)
}

#start by loading in all the analysis packages you'll use with the code below
#loading packages and libraries
if (!requireNamespace("mpra", quietly = TRUE))
  source("https://bioconductor.org/biocLite.R")
BiocInstaller::biocLite(c("mpra"))

if (!requireNamespace("dplyr", quietly = TRUE))
  install.packages("dplyr")

if (!requireNamespace("data.table", quietly = TRUE))
  install.packages("data.table")

library(mpra)
library(dplyr)
library(data.table)

#Incorporate a lookup table strategy to assign long column names to easier variale names
old_utr_vars <- c("barcode", "IVT2_CUpostL_S19_L008_R1_001utr.count.txt", "IVT2_CUpreL_S18_L008_R1_001utr.count.txt",
                  "RNA_CU_PostL_S58_L007_R1_001utr.count.txt", "RNA_CU_PreL_S57_L007_R1_001utr.count.txt",
                  "IVT2_CUpostR_S106_L003_R1_001utr.count.txt", "IVT2_CUpreR_S105_L003_R1_001utr.count.txt",
                  "RNA_CU_PostR_S66_L008_R1_001utr.count.txt", "RNA_CU_PreR_S65_L008_R1_001utr.count.txt")

old_cds_vars <- c("barcode", "IVT2_CUpostL_S19_L008_R1_001cds.count.txt", "IVT2_CUpreL_S18_L008_R1_001cds.count.txt",
                  "RNA_CU_PostL_S58_L007_R1_001cds.count.txt", "RNA_CU_PreL_S57_L007_R1_001cds.count.txt",
                  "IVT2_CUpostR_S106_L003_R1_001cds.count.txt", "IVT2_CUpreR_S105_L003_R1_001cds.count.txt",
                  "RNA_CU_PostR_S66_L008_R1_001cds.count.txt", "RNA_CU_PreR_S65_L008_R1_001cds.count.txt")

newvars <- c("barcode", "IVT_postL", "IVT_preL",
             "RNA_postL", "RNA_preL",
             "IVT_postR", "IVT_preR",
             "RNA_postR", "RNA_preR")

lookup = data.frame(old_utr_vars, old_cds_vars, newvars)

names(TF_UTR) <- lookup[match(names(TF_UTR), lookup$old_utr_vars),"newvars"]
names(TF_CDS) <- lookup[match(names(TF_CDS), lookup$old_cds_vars),"newvars"]

head(TF_UTR)
head(TF_CDS)

#create dataframes that mpralm can use, filter with a 32 count cut off for DNA pre samples
#filter for either left or right > 32 read counts and merge with barcode assignment dataframe
cds_32 <- filter(TF_CDS,
                  TF_CDS$IVT_preL > 32 | TF_CDS$IVT_preR > 32)

utr_32 <- filter(TF_UTR,
                  TF_UTR$IVT_preL > 32 | TF_UTR$IVT_preR > 32)

cds_32_guides <- merge(cds_32, grna.assign.barcode.grna.good, by="barcode")

utr_32_guides <- merge(utr_32, grna.assign.barcode.grna.good, by="barcode")

#Preparing the DNA, RNA, and Element ID (eid) dataframes for His4 mpralm analysis
dna <- cds_32_guides
dna <-data.frame(row.names=dna$barcode, dna$IVT_preL,
                 dna$IVT_preR,
                 dna$IVT_postL,
                 dna$IVT_postR)
names(dna) <- c("preL", "preR", "postL", "postR")

rna <- cds_32_guides
rna <-data.frame(row.names=rna$barcode, rna$RNA_preL,
                 rna$RNA_preR,
                 rna$RNA_postL,
                 rna$RNA_postR)
names(rna) <- c("preL", "preR", "postL", "postR")

eid <- as.character(cds_32_guides$guide)

#mpra dataset construction, this will take the three above dataframes and create an mpra matrix
mpraset <- MPRASet(DNA = dna, RNA = rna, eid = eid, eseq = NULL, barcode = NULL)
mpraset

#mpra data set analysis, this block of code may be expected to take ~10ish mins to run
design <- data.frame(intcpt = 1,
                     pre = grepl("pre", colnames(mpraset)),
                     turb = grepl("L", colnames(mpraset)))

mpralm_fit <- mpralm(object = mpraset, design = design,
                     aggregate = "sum", normalize = TRUE,
                     model_type = "indep_groups", plot = TRUE)
toptab <- topTable(mpralm_fit, coef = "pre", number = Inf)
head(toptab)

#This example code will make modifications to the output analysis dataframe. 
#You'll need to replace "gene_name_mpralm" variable with the name you've chosen for your dataframe

cds_sum_mpralm <- toptab
setDT(cds_sum_mpralm, keep.rownames = TRUE)[]
names(cds_sum_mpralm)[names(cds_sum_mpralm) == "rn"] <- "Guide"
cds_sum_mpralm <- merge(cds_sum_mpralm, guide.good.targets, by="Guide")
cds_sum_mpralm$gene <- sgd[match(cds_sum_mpralm$Yorf1, sgd$name), "gene"]
cds_sum_mpralm$desc <- sgd[match(cds_sum_mpralm$Yorf1, sgd$name), "desc"]
head(cds_sum_mpralm)

#Same mpralm analysis for pgk1 dataset
dna <- utr_32_guides
dna <-data.frame(row.names=dna$barcode, dna$IVT_preL,
                 dna$IVT_preR,
                 dna$IVT_postL,
                 dna$IVT_postR)
names(dna) <- c("preL", "preR", "postL", "postR")

rna <- utr_32_guides
rna <-data.frame(row.names=rna$barcode, rna$RNA_preL,
                 rna$RNA_preR,
                 rna$RNA_postL,
                 rna$RNA_postR)
names(rna) <- c("preL", "preR", "postL", "postR")

eid <- as.character(utr_32_guides$guide)

#mpra dataset construction, this will take the three above dataframes and create an mpra matrix
mpraset <- MPRASet(DNA = dna, RNA = rna, eid = eid, eseq = NULL, barcode = NULL)
mpraset

#mpra data set analysis, this block of code may be expected to take ~10ish mins to run
design <- data.frame(intcpt = 1,
                     pre = grepl("pre", colnames(mpraset)),
                     turb = grepl("L", colnames(mpraset)))

mpralm_fit <- mpralm(object = mpraset, design = design,
                     aggregate = "sum", normalize = TRUE,
                     model_type = "indep_groups", plot = TRUE)
toptab <- topTable(mpralm_fit, coef = "pre", number = Inf)
head(toptab)

#This example code will make modifications to the output analysis dataframe. 
#You'll need to replace "gene_name_mpralm" variable with the name you've chosen for your dataframe

utr_sum_mpralm <- toptab
setDT(utr_sum_mpralm, keep.rownames = TRUE)[]
names(utr_sum_mpralm)[names(utr_sum_mpralm) == "rn"] <- "Guide"
utr_sum_mpralm <- merge(utr_sum_mpralm, guide.good.targets, by="Guide")
utr_sum_mpralm$gene <- sgd[match(utr_sum_mpralm$Yorf1, sgd$name), "gene"]
utr_sum_mpralm$desc <- sgd[match(utr_sum_mpralm$Yorf1, sgd$name), "desc"]
head(utr_sum_mpralm)

#saving mpralm analysis
write.table(cds_sum_mpralm, "~/CiBER_seq_package/all_raw_fasta_gz/GCN4_CDS_UTR/cds_sum_mpralm.txt", sep="\t")
write.table(utr_sum_mpralm, "~/CiBER_seq_package/all_raw_fasta_gz/GCN4_CDS_UTR/utr_sum_mpralm.txt", sep="\t")

