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
if (!exists("HTS1_RPC31_rep1")) {
  HTS1_RPC31 <- read.delim("~/CiBER_seq_package/all_raw_fasta_gz/HTS1_RPC31/all_HTS1_RPC31_counts.txt",
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

#merge dataframes
old <- c("barcode", "IVT_HTS1_Post_S2_R1_001.count.txt", "IVT_RPC31_Post_S4_R1_001.count.txt",
         "IVT_HTS1_Pre_S1_R1_001.count.txt", "IVT_RPC31_Pre_S3_R1_001.count.txt",
         "RNA_HTS1_Post_S6_R1_001.count.txt",
         "RNA_RPC31_Post_S8_R1_001.count.txt", "RNA_HTS1_Pre_S5_R1_001.count.txt",
         "RNA_RPC31_Pre_S7_R1_001.count.txt", "HTS1_DNA_post.count.txt", "RPC31_DNA_post.count.txt", "HTS1_DNA_pre.count.txt",
         "RPC31_DNA_pre.count.txt", "HTS1_RNA_post.count.txt", "RPC31_RNA_post.count.txt", "HTS1_RNA_pre.count.txt",
         "RPC31_RNA_pre.count.txt")

newvars <- c("barcode", "IVT_HTS1_post1", "IVT_RPC31_post1",
         "IVT_HTS1_pre1", "IVT_RPC31_pre1",
         "RNA_HTS1_post1",
         "RNA_RPC31_post1", "RNA_HTS1_pre1",
         "RNA_RPC31_pre1", "IVT_HTS1_post3", "IVT_RPC31_post3", "IVT_HTS1_pre3",
         "IVT_RPC31_pre3", "RNA_HTS1_post3", "RNA_RPC31_post3", "RNA_HTS1_pre3",
         "RNA_RPC31_pre3")

lookup = data.frame(old, newvars)
head(lookup)

names(HTS1_RPC31) <- lookup[match(names(HTS1_RPC31), lookup$old),"newvars"]

head(HTS1_RPC31)
all <- HTS1_RPC31

all_HTS1 <- data.frame(all$barcode, all$IVT_HTS1_post1, all$IVT_HTS1_pre1, all$RNA_HTS1_post1, 
                       all$RNA_HTS1_pre1, all$IVT_HTS1_post3, 
                       all$IVT_HTS1_pre3, all$RNA_HTS1_post3, 
                       all$RNA_HTS1_pre3)

all_RPC31 <- data.frame(all$barcode, all$IVT_RPC31_post1, all$IVT_RPC31_pre1, all$RNA_RPC31_post1, 
                        all$RNA_RPC31_pre1, all$IVT_RPC31_post3, 
                        all$IVT_RPC31_pre3, all$RNA_RPC31_post3, all$RNA_RPC31_pre3)

names(all_HTS1) <- gsub(x = names(all_HTS1), pattern = "all.", replacement = "")
names(all_RPC31) <- gsub(x = names(all_RPC31), pattern = "all.", replacement = "")

all_HTS1[is.na(all_HTS1)] <- 0
all_RPC31[is.na(all_RPC31)] <- 0

head(all_HTS1)
head(all_RPC31)

#create dataframes that mpralm can use, pool seq runs, filter with a 32 count cut off for DNA pre samples
#filter for either left or right > 32 read counts and merge with barcode assignment dataframe
all_HTS1_32 <- filter(all_HTS1, 
                  all_HTS1$IVT_HTS1_pre1 > 32 | all_HTS1$IVT_HTS1_pre3 > 32)

all_RPC31_32 <- filter(all_RPC31,
                  all_RPC31$IVT_RPC31_pre1 > 32 | all_RPC31$IVT_RPC31_pre3 > 32)

HTS1_32_guides <- merge(all_HTS1_32, grna.assign.barcode.grna.good, by="barcode")

RPC31_32_guides <- merge(all_RPC31_32, grna.assign.barcode.grna.good, by="barcode")

#Preparing the DNA, RNA, and Element ID (eid) dataframes for His4 mpralm analysis

dna <- HTS1_32_guides
dna <-data.frame(row.names=dna$barcode, dna$IVT_HTS1_pre1,
                 dna$IVT_HTS1_pre3,
                 dna$IVT_HTS1_post1,
                 dna$IVT_HTS1_post3)
names(dna) <- c("pre1", "pre3", "post1", "post3")

rna <- HTS1_32_guides
rna <-data.frame(row.names=rna$barcode, rna$RNA_HTS1_pre1,
                 rna$RNA_HTS1_pre3,
                 rna$RNA_HTS1_post1,
                 rna$RNA_HTS1_post3)
names(rna) <- c("pre1", "pre3", "post1", "post3")

eid <- as.character(HTS1_32_guides$guide)

#mpra dataset construction, this will take the three above dataframes and create an mpra matrix
mpraset <- MPRASet(DNA = dna, RNA = rna, eid = eid, eseq = NULL, barcode = NULL)
mpraset

#mpra data set analysis, this block of code may be expected to take ~10ish mins to run
design <- data.frame(intcpt = 1,
                     pre = grepl("pre", colnames(mpraset)),
                     rep3 = grepl("3", colnames(mpraset)))

mpralm_fit <- mpralm(object = mpraset, design = design,
                     aggregate = "sum", normalize = TRUE,
                     model_type = "indep_groups", plot = TRUE)
toptab <- topTable(mpralm_fit, coef = 2, number = Inf)
head(toptab)

#Modifications to the output analysis dataframe. 

HTS1_sum_mpralm <- toptab
setDT(HTS1_sum_mpralm, keep.rownames = TRUE)[]
names(HTS1_sum_mpralm)[names(HTS1_sum_mpralm) == "rn"] <- "Guide"
HTS1_sum_mpralm <- merge(HTS1_sum_mpralm, guide.good.targets, by="Guide")
HTS1_sum_mpralm$gene <- sgd[match(HTS1_sum_mpralm$Yorf1, sgd$name), "gene"]
HTS1_sum_mpralm$desc <- sgd[match(HTS1_sum_mpralm$Yorf1, sgd$name), "desc"]
head(HTS1_sum_mpralm)

#Same mpralm analysis for RPC31 dataset
dna <- RPC31_32_guides
dna <-data.frame(row.names=dna$barcode, dna$IVT_RPC31_pre1,
                 dna$IVT_RPC31_pre3,
                 dna$IVT_RPC31_post1,
                 dna$IVT_RPC31_post3)
names(dna) <- c("pre1", "pre3", "post1", "post3")

rna <- RPC31_32_guides
rna <-data.frame(row.names=rna$barcode, rna$RNA_RPC31_pre1,
                 rna$RNA_RPC31_pre3,
                 rna$RNA_RPC31_post1,
                 rna$RNA_RPC31_post3)
names(rna) <- c("pre1", "pre3", "post1", "post3")

eid <- as.character(RPC31_32_guides$guide)

#mpra dataset construction, this will take the three above dataframes and create an mpra matrix
mpraset <- MPRASet(DNA = dna, RNA = rna, eid = eid, eseq = NULL, barcode = NULL)
mpraset

#mpra data set analysis, this block of code may be expected to take ~10ish mins to run
design <- data.frame(intcpt = 1,
                     pre = grepl("pre", colnames(mpraset)),
                     rep3 = grepl("3", colnames(mpraset)))

mpralm_fit <- mpralm(object = mpraset, design = design,
                     aggregate = "sum", normalize = TRUE,
                     model_type = "indep_groups", plot = TRUE)
toptab <- topTable(mpralm_fit, coef = 2, number = Inf)
head(toptab)

#This example code will make modifications to the output analysis dataframe. 
#You'll need to replace "gene_name_mpralm" variable with the name you've chosen for your dataframe

RPC31_sum_mpralm <- toptab
setDT(RPC31_sum_mpralm, keep.rownames = TRUE)[]
names(RPC31_sum_mpralm)[names(RPC31_sum_mpralm) == "rn"] <- "Guide"
RPC31_sum_mpralm <- merge(RPC31_sum_mpralm, guide.good.targets, by="Guide")
RPC31_sum_mpralm$gene <- sgd[match(RPC31_sum_mpralm$Yorf1, sgd$name), "gene"]
RPC31_sum_mpralm$desc <- sgd[match(RPC31_sum_mpralm$Yorf1, sgd$name), "desc"]
head(RPC31_sum_mpralm)

#saving mpralm analysis
write.table(HTS1_sum_mpralm, "~/CiBER_seq_package/all_raw_fasta_gz/HTS1_RPC31/HTS1_sum_mpralm.txt", sep="\t")
write.table(RPC31_sum_mpralm, "~/CiBER_seq_package/all_raw_fasta_gz/HTS1_RPC31/RPC31_sum_mpralm.txt", sep="\t")
