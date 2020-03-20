#!/usr/bin/env Rscript

#download sgd, bc_gRNA_assignment table, and gRNA_re-assignment table, and dataframes

if (!file.exists("SGD_features.tab")) {
  sgd <- download.file('https://downloads.yeastgenome.org/curation/chromosomal_feature/SGD_features.tab', destfile="SGD_features.tab")
}
sgd <- read.delim("SGD_features.tab", header=FALSE, quote="",
                  col.names=c("sgdid", "type", "qual", "name", "gene", "alias",
                              "parent", "sgdid2", "chrom", "start", "end",
                              "strand", "genpos", "cver", "sver", "desc"))

options(stringsAsFactors=FALSE)
if (!exists("grna.assign.barcode.grna.good")) {
  grna.assign.barcode.grna.good <- read.delim("~/CiBER_seq_package/all_raw_fastq/PE_bc_gRNA_assignment/grna-assign-barcode-grna-good.txt",
                                stringsAsFactors=FALSE)
}

grna.assign.barcode.grna.good$guide <- as.character(grna.assign.barcode.grna.good$guide)
grna.assign.barcode.grna.good$guide[grna.assign.barcode.grna.good$guide == "No_gRNA"] <- paste("No_gRNA", seq(1:788), sep=" ")

options(stringsAsFactors=FALSE)
if (!exists("guide.good.targets")) {
  guide.good.targets <- read.delim("./guide.good.targets_plus_empty.txt",
                                stringsAsFactors=FALSE)
}

# Load in raw count data and relevant gRNA barcode assignment and gRNA assignment of target files
options(stringsAsFactors=FALSE)
if (!exists("his4_seq1")) {
  his4_seq1 <- read.delim("~/CiBER_seq_package/all_raw_fastq/HIS4_PGK1_3AT/all_his4_seq1_counts.txt",
                                stringsAsFactors=FALSE)
}

options(stringsAsFactors=FALSE)
if (!exists("pgk1_seq1")) {
  pgk1_seq1 <- read.delim("~/CiBER_seq_package/all_raw_fastq/HIS4_PGK1_3AT/all_pgk1_seq1_counts.txt",
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
old_pgk1_vars <- c("barcode", "IVT_3AT_L_S21_L008_R1_001pgk1.count.txt", "IVT_3AT_R_S22_L008_R1_001pgk1.count.txt", "IVT_PostL_S19_L008_R1_001pgk1.count.txt",
                   "IVT_PostR_S20_L008_R1_001pgk1.count.txt", "IVT_PreL_S17_L008_R1_001pgk1.count.txt", "IVT_PreR_S18_L008_R1_001pgk1.count.txt",
                   "RNA_3AT_L_S27_L008_R1_001pgk1.count.txt", "RNA_3AT_R_S28_L008_R1_001pgk1.count.txt", "RNA_PostL_S25_L008_R1_001pgk1.count.txt",
                   "RNA_PostR_S26_L008_R1_001pgk1.count.txt", "RNA_PreL_S23_L008_R1_001pgk1.count.txt", "RNA_PreR_S24_L008_R1_001pgk1.count.txt")

old_his4_vars <- c("barcode", "IVT_3AT_L_S21_L008_R1_001his4.count.txt", "IVT_3AT_R_S22_L008_R1_001his4.count.txt", "IVT_PostL_S19_L008_R1_001his4.count.txt",
                   "IVT_PostR_S20_L008_R1_001his4.count.txt", "IVT_PreL_S17_L008_R1_001his4.count.txt", "IVT_PreR_S18_L008_R1_001his4.count.txt",
                   "RNA_3AT_L_S27_L008_R1_001his4.count.txt", "RNA_3AT_R_S28_L008_R1_001his4.count.txt", "RNA_PostL_S25_L008_R1_001his4.count.txt",
                   "RNA_PostR_S26_L008_R1_001his4.count.txt", "RNA_PreL_S23_L008_R1_001his4.count.txt", "RNA_PreR_S24_L008_R1_001his4.count.txt")

newvars <- c("barcode", "IVT_3AT_L", "IVT_3AT_R", "IVT_postL", "IVT_postR", "IVT_preL", "IVT_preR", "RNA_3AT_L", "RNA_3AT_R", "RNA_postL", "RNA_postR", "RNA_preL", "RNA_preR")

lookup = data.frame(old_pgk1_vars, old_his4_vars, newvars)

names(his4_seq1) <- lookup[match(names(his4_seq1), lookup$old_his4_vars),"newvars"]
names(pgk1_seq1) <- lookup[match(names(pgk1_seq1), lookup$old_pgk1_vars),"newvars"]

head(his4_seq1)
head(pgk1_seq1)

#create dataframes that mpralm can use, filter with a 32 count cut off for DNA pre samples
#filter for either left or right > 32 read counts and merge with barcode assignment dataframe
pgk1_32 <- filter(pgk1_seq1, 
                  pgk1_seq1$IVT_preL > 32 | pgk1_seq1$IVT_preR > 32)

his4_32 <- filter(his4_seq1,
                  his4_seq1$IVT_preL > 32 | his4_seq1$IVT_preR > 32)

pgk1_32_guides <- merge(pgk1_32, grna.assign.barcode.grna.good, by="barcode")

his4_32_guides <- merge(his4_32, grna.assign.barcode.grna.good, by="barcode")

#Preparing the DNA, RNA, and Element ID (eid) dataframes for His4 mpralm analysis
dna <- his4_32_guides
dna <-data.frame(row.names=dna$barcode, dna$IVT_3AT_L, 
                 dna$IVT_3AT_R, 
                 dna$IVT_preL, 
                 dna$IVT_preR)
names(dna) <- c("3AT_L", "3AT_R", "preL", "preR")

rna <- his4_32_guides
rna <-data.frame(row.names=rna$barcode, rna$RNA_3AT_L, 
                 rna$RNA_3AT_R, 
                 rna$RNA_preL, 
                 rna$RNA_preR)
names(rna) <- c("3AT_L", "3AT_R", "preL", "preR")

eid <- as.character(his4_32_guides$guide)

#mpra dataset construction, this will take the three above dataframes and create an mpra matrix
mpraset <- MPRASet(DNA = dna, RNA = rna, eid = eid, eseq = NULL, barcode = NULL)
mpraset

#mpra data set analysis, this block of code may be expected to take ~10ish mins to run
design <- data.frame(intcpt = 1,
                     pre = grepl("pre", colnames(mpraset)),
                     right = grepl("R", colnames(mpraset)))

mpralm_fit <- mpralm(object = mpraset, design = design,
                     aggregate = "sum", normalize = TRUE,
                     model_type = "indep_groups", plot = TRUE)
toptab <- topTable(mpralm_fit, coef = 2, number = Inf)
head(toptab)

#This example code will make modifications to the output analysis dataframe. 
#You'll need to replace "gene_name_mpralm" variable with the name you've chosen for your dataframe

his4_sum_mpralm <- toptab 
setDT(his4_sum_mpralm, keep.rownames = TRUE)[]
names(his4_sum_mpralm)[names(his4_sum_mpralm) == "rn"] <- "Guide"
his4_sum_mpralm <- merge(his4_sum_mpralm, guide.good.targets, by="Guide")
his4_sum_mpralm$gene <- sgd[match(his4_sum_mpralm$Yorf1, sgd$name), "gene"]
his4_sum_mpralm$desc <- sgd[match(his4_sum_mpralm$Yorf1, sgd$name), "desc"]
head(his4_sum_mpralm)

#Same mpralm analysis for pgk1 dataset
dna <- pgk1_32_guides
dna <-data.frame(row.names=dna$barcode, dna$IVT_3AT_L, 
                 dna$IVT_3AT_R, 
                 dna$IVT_preL, 
                 dna$IVT_preR)
names(dna) <- c("3AT_L", "3AT_R", "preL", "preR")

rna <- pgk1_32_guides
rna <-data.frame(row.names=rna$barcode, rna$RNA_3AT_L, 
                 rna$RNA_3AT_R, 
                 rna$RNA_preL, 
                 rna$RNA_preR)
names(rna) <- c("3AT_L", "3AT_R", "preL", "preR")

eid <- as.character(pgk1_32_guides$guide)

#mpra dataset construction, this will take the three above dataframes and create an mpra matrix
mpraset <- MPRASet(DNA = dna, RNA = rna, eid = eid, eseq = NULL, barcode = NULL)
mpraset

#mpra data set analysis, this block of code may be expected to take ~10ish mins to run
design <- data.frame(intcpt = 1,
                     pre = grepl("pre", colnames(mpraset)),
                     right = grepl("R", colnames(mpraset)))

mpralm_fit <- mpralm(object = mpraset, design = design,
                     aggregate = "sum", normalize = TRUE,
                     model_type = "indep_groups", plot = TRUE)
toptab <- topTable(mpralm_fit, coef = 2, number = Inf)
head(toptab)

#This example code will make modifications to the output analysis dataframe. 
#You'll need to replace "gene_name_mpralm" variable with the name you've chosen for your dataframe

pgk1_sum_mpralm <- toptab 
setDT(pgk1_sum_mpralm, keep.rownames = TRUE)[]
names(pgk1_sum_mpralm)[names(pgk1_sum_mpralm) == "rn"] <- "Guide"
pgk1_sum_mpralm <- merge(pgk1_sum_mpralm, guide.good.targets, by="Guide")
pgk1_sum_mpralm$gene <- sgd[match(pgk1_sum_mpralm$Yorf1, sgd$name), "gene"]
pgk1_sum_mpralm$desc <- sgd[match(pgk1_sum_mpralm$Yorf1, sgd$name), "desc"]
head(pgk1_sum_mpralm)

#saving mpralm analysis
write.table(pgk1_sum_mpralm, "~/CiBER_seq_package/all_raw_fastq/HIS4_PGK1_3AT/pgk1_prev3AT_sum_mpralm.txt", sep="\t")
write.table(his4_sum_mpralm, "~/CiBER_seq_package/all_raw_fastq/HIS4_PGK1_3AT/his4_prev3AT_sum_mpralm.txt", sep="\t")
