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
  grna.assign.barcode.grna.good <- read.delim("~/CiBER_seq_package/all_raw_fastq/PE_bc_gRNA_assignment/grna-assign-barcode-grna-good.txt",
                                              stringsAsFactors=FALSE)
}

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

# Load in raw count data and relevant gRNA barcode assignment and gRNA assignment of target files
#options(stringsAsFactors=FALSE)
#if (!exists("his4_moreseq")) {
#  his4_moreseq <- read.delim("~/CiBER_seq_package/all_raw_fastq/HIS4_PGK1_pooled/more_seq/all_his4_moreseq_counts.txt",
#                             stringsAsFactors=FALSE)
#}

# Load in raw count data and relevant gRNA barcode assignment and gRNA assignment of target files
#options(stringsAsFactors=FALSE)
#if (!exists("HTS1_RPC31_rep1")) {
#  HTS1_RPC31 <- read.delim("~/CiBER_seq_package/all_raw_fastq/HTS1_RPC31/all_HTS1_RPC31_counts.txt",
#                           stringsAsFactors=FALSE)
#}

# Load in raw count data and relevant gRNA barcode assignment and gRNA assignment of target files
options(stringsAsFactors=FALSE)
if (!exists("HTS1_counts")) {
  all_HTS1 <- read.delim("~/CiBER_seq_package/all_raw_fastq/HTS1_RPC31/all_HTS1counts.txt",
                         stringsAsFactors=FALSE)
}

if (!exists("RPC31_counts")) {
  all_RPC31 <- read.delim("~/CiBER_seq_package/all_raw_fastq/HTS1_RPC31/all_RPC31counts.txt",
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
grna.assign.barcode.grna.good$guide <- as.character(grna.assign.barcode.grna.good$guide)
grna.assign.barcode.grna.good$guide[grna.assign.barcode.grna.good$guide == "No_gRNA"] <- paste("No_gRNA", seq(1:788), sep=" ")

#old_his4_vars <- c("barcode", "IVT_3AT_L_S21_L008_R1_001his4.count.txt", "IVT_3AT_R_S22_L008_R1_001his4.count.txt", "IVT_PostL_S19_L008_R1_001his4.count.txt",
#                   "IVT_PostR_S20_L008_R1_001his4.count.txt", "IVT_PreL_S17_L008_R1_001his4.count.txt", "IVT_PreR_S18_L008_R1_001his4.count.txt",
#                   "RNA_3AT_L_S27_L008_R1_001his4.count.txt", "RNA_3AT_R_S28_L008_R1_001his4.count.txt", "RNA_PostL_S25_L008_R1_001his4.count.txt",
#                   "RNA_PostR_S26_L008_R1_001his4.count.txt", "RNA_PreL_S23_L008_R1_001his4.count.txt", "RNA_PreR_S24_L008_R1_001his4.count.txt")
#
#newvars <- c("barcode", "IVT_3AT_L", "IVT_3AT_R", "IVT_postL", "IVT_postR", "IVT_preL", "IVT_preR", "RNA_3AT_L", "RNA_3AT_R", "RNA_postL", "RNA_postR", "RNA_preL", "RNA_preR")
#
#old_his4_moreseq <- c("barcode", "IVT_PostL_S52_L008_R1_001his4.count.txt", "IVT_PostR_S53_L008_R1_001his4.count.txt",
#                      "IVT_PreL_S50_L008_R1_001his4.count.txt", "IVT_PreR_S51_L008_R1_001his4.count.txt", "RNA_PostL_S56_L008_R1_001his4.count.txt",
#                      "RNA_PostR_S57_L008_R1_001his4.count.txt", "RNA_PreL_S54_L008_R1_001his4.count.txt",
#                      "RNA_PreR_S55_L008_R1_001his4.count.txt")
#
#newvars2 <- c("barcode", "IVT_postL2", "IVT_postR2", "IVT_preL2", "IVT_preR2", "RNA_postL2", "RNA_postR2", "RNA_preL2", "RNA_preR2")
#
#lookup1 = data.frame(old_his4_vars, newvars)
#lookup2 = data.frame(old_his4_moreseq, newvars2)
#
#names(his4_seq1) <- lookup1[match(names(his4_seq1), lookup1$old_his4_vars),"newvars"]
#names(his4_moreseq) <- lookup2[match(names(his4_moreseq), lookup2$old_his4_moreseq),"newvars2"]

names(his4_seq1) <- gsub(x = names(his4_seq1), pattern = "his4.count.txt", replacement = "")
names(his4_seq1) <- gsub(x = names(his4_seq1), pattern = "PH_", replacement = "")
names(his4_seq1) <- gsub(x = names(his4_seq1), pattern = "postL", replacement = "post1")
names(his4_seq1) <- gsub(x = names(his4_seq1), pattern = "postR", replacement = "post3")
head(his4_seq1)

#HTS1 and RPC31 column rename and split by dataset
#old <- c("barcode", "IVT_HTS1_Post_S2_R1_001.count.txt", "IVT_RPC31_Post_S4_R1_001.count.txt",
#         "IVT_HTS1_Pre_S1_R1_001.count.txt", "IVT_RPC31_Pre_S3_R1_001.count.txt",
#         "RNA_HTS1_Post_S6_R1_001.count.txt",
#         "RNA_RPC31_Post_S8_R1_001.count.txt", "RNA_HTS1_Pre_S5_R1_001.count.txt",
#         "RNA_RPC31_Pre_S7_R1_001.count.txt", "HTS1_DNA_post.count.txt", "RPC31_DNA_post.count.txt", "HTS1_DNA_pre.count.txt",
#         "RPC31_DNA_pre.count.txt", "HTS1_RNA_post.count.txt", "RPC31_RNA_post.count.txt", "HTS1_RNA_pre.count.txt",
#         "RPC31_RNA_pre.count.txt")
#
#newvars <- c("barcode", "IVT_HTS1_post1", "IVT_RPC31_post1",
#             "IVT_HTS1_pre1", "IVT_RPC31_pre1",
#             "RNA_HTS1_post1",
#             "RNA_RPC31_post1", "RNA_HTS1_pre1",
#             "RNA_RPC31_pre1", "IVT_HTS1_post3", "IVT_RPC31_post3", "IVT_HTS1_pre3",
#             "IVT_RPC31_pre3", "RNA_HTS1_post3", "RNA_RPC31_post3", "RNA_HTS1_pre3",
#             "RNA_RPC31_pre3")
#
#lookup = data.frame(old, newvars)
#head(lookup)
#
#names(HTS1_RPC31) <- lookup[match(names(HTS1_RPC31), lookup$old),"newvars"]
#
#head(HTS1_RPC31)
#all <- HTS1_RPC31
#
#all_HTS1 <- data.frame(all$barcode, all$IVT_HTS1_post1, all$IVT_HTS1_pre1, all$RNA_HTS1_post1, 
#                       all$RNA_HTS1_pre1, all$IVT_HTS1_post3, 
#                       all$IVT_HTS1_pre3, all$RNA_HTS1_post3, 
#                       all$RNA_HTS1_pre3)
#
#all_RPC31 <- data.frame(all$barcode, all$IVT_RPC31_post1, all$IVT_RPC31_pre1, all$RNA_RPC31_post1, 
#                        all$RNA_RPC31_pre1, all$IVT_RPC31_post3, 
#                        all$IVT_RPC31_pre3, all$RNA_RPC31_post3, all$RNA_RPC31_pre3)
#
#names(all_HTS1) <- gsub(x = names(all_HTS1), pattern = "all.", replacement = "")
#names(all_RPC31) <- gsub(x = names(all_RPC31), pattern = "all.", replacement = "")
#
#all_HTS1[is.na(all_HTS1)] <- 0
#all_RPC31[is.na(all_RPC31)] <- 0
#
#head(all_HTS1)
#head(all_RPC31)

#create dataframes that mpralm can use, pool seq runs, filter with a 32 count cut off for DNA pre samples
#filter for either left or right > 32 read counts and merge with barcode assignment dataframe
#his4 <- merge(his4_seq1, his4_moreseq, by="barcode", all=TRUE)
#his4[is.na(his4)] <- 0
#head(his4)

#all_his4 <- data.frame(his4$barcode)
#colnames(all_his4)[colnames(all_his4)=="his4.barcode"] <- "barcode"

#all_his4$IVT_preL <- his4$IVT_preL + his4$IVT_preL2
#all_his4$IVT_postL <- his4$IVT_postL + his4$IVT_postL2
#all_his4$IVT_preR <- his4$IVT_preR + his4$IVT_preR2
#all_his4$IVT_postR <- his4$IVT_postR + his4$IVT_postR2
#all_his4$RNA_preL <- his4$RNA_preL + his4$RNA_preL2
#all_his4$RNA_postL <- his4$RNA_postL + his4$RNA_postL2
#all_his4$RNA_preR <- his4$RNA_preR + his4$RNA_preR2
#all_his4$RNA_postR <- his4$RNA_postR + his4$RNA_postR2
#
#head(all_his4)

all_his4 <- his4_seq1

his4_32 <- filter(all_his4,
                  all_his4$IVT_preL > 32 | all_his4$IVT_preR > 32)

#filter for either left or right > 32 read counts and merge with barcode assignment dataframe
names(all_HTS1) <- gsub(x = names(all_HTS1), pattern = ".count.txt", replacement = "")
names(all_RPC31) <- gsub(x = names(all_RPC31), pattern = ".count.txt", replacement = "")

all_HTS1_32 <- filter(all_HTS1, 
                      all_HTS1$IVT_HTS1_pre1 > 32 | all_HTS1$IVT_HTS1_pre3 > 32)

all_RPC31_32 <- filter(all_RPC31,
                       all_RPC31$IVT_RPC31_pre1 > 32 | all_RPC31$IVT_RPC31_pre3 > 32)

#prepare His4post vs. HTS1post and His4post vs RPC31post
#his4post_32 <- his4_32[c(1,3,5,7,9)]
#names(his4post_32) <- c("barcode", "IVT_post1", "IVT_post3", "RNA_post1", "RNA_post3")
#HTS1post_32 <- all_HTS1_32[c(1,2,4,6,8)]
#RPC31post_32 <- all_RPC31_32[c(1,2,4,6,8)]

hisvHTS1 <- merge(his4_32, all_HTS1_32, by="barcode")
hisvRPC31 <- merge(his4_32, all_RPC31_32, by="barcode")

hisvHTS1_guides <- merge(hisvHTS1, grna.assign.barcode.grna.good, by="barcode")
hisvRPC31_guides <- merge(hisvRPC31, grna.assign.barcode.grna.good, by="barcode")

head(hisvHTS1_guides)
head(hisvRPC31_guides)

#Preparing the DNA, RNA, and Element ID (eid) dataframes for His4 post vs HTS1 post mpralm analysis
dna <- hisvHTS1_guides
dna <-data.frame(row.names=dna$barcode, dna$IVT_post1, 
                 dna$IVT_post3, 
                 dna$IVT_HTS1_post1, 
                 dna$IVT_HTS1_post3)
names(dna) <- c("post1", "post3", "HTSpost1", "HTSpost3")

rna <- hisvHTS1_guides
rna <-data.frame(row.names=rna$barcode, rna$RNA_post1, 
                 rna$RNA_post3, 
                 rna$RNA_HTS1_post1, 
                 rna$RNA_HTS1_post3)
names(rna) <- c("post1", "post3", "HTSpost1", "HTSpost3")

eid <- as.character(hisvHTS1_guides$guide)

#mpra dataset construction, this will take the three above dataframes and create an mpra matrix
mpraset <- MPRASet(DNA = dna, RNA = rna, eid = eid, eseq = NULL, barcode = NULL)
mpraset

#mpra data set analysis, this block of code may be expected to take ~10ish mins to run
design <- data.frame(intcpt = 1,
                     post = grepl("HTS", colnames(mpraset)),
                     right = grepl("3", colnames(mpraset)))

mpralm_fit <- mpralm(object = mpraset, design = design,
                     aggregate = "sum", normalize = TRUE,
                     model_type = "indep_groups", plot = TRUE)
toptab <- topTable(mpralm_fit, coef = 2, number = Inf)
head(toptab)

#This example code will make modifications to the output analysis dataframe. 
#You'll need to replace "gene_name_mpralm" variable with the name you've chosen for your dataframe

his4postvHTS1post_sum_mpralm <- toptab 
setDT(his4postvHTS1post_sum_mpralm, keep.rownames = TRUE)[]
names(his4postvHTS1post_sum_mpralm)[names(his4postvHTS1post_sum_mpralm) == "rn"] <- "Guide"
his4postvHTS1post_sum_mpralm <- merge(his4postvHTS1post_sum_mpralm, guide.good.targets, by="Guide")
his4postvHTS1post_sum_mpralm$gene <- sgd[match(his4postvHTS1post_sum_mpralm$Yorf1, sgd$name), "gene"]
his4postvHTS1post_sum_mpralm$desc <- sgd[match(his4postvHTS1post_sum_mpralm$Yorf1, sgd$name), "desc"]
head(his4postvHTS1post_sum_mpralm)

#_______________
#Preparing the DNA, RNA, and Element ID (eid) dataframes for His4 post vs RPC31 post mpralm analysis
dna <- hisvRPC31_guides
dna <-data.frame(row.names=dna$barcode, dna$IVT_post1, 
                 dna$IVT_post3, 
                 dna$IVT_RPC31_post1, 
                 dna$IVT_RPC31_post3)
names(dna) <- c("post1", "post3", "RPCpost1", "RPCpost3")

rna <- hisvRPC31_guides
rna <-data.frame(row.names=rna$barcode, rna$RNA_post1, 
                 rna$RNA_post3, 
                 rna$RNA_RPC31_post1, 
                 rna$RNA_RPC31_post3)
names(rna) <- c("post1", "post3", "RPCpost1", "RPCpost3")

eid <- as.character(hisvRPC31_guides$guide)

#mpra dataset construction, this will take the three above dataframes and create an mpra matrix
mpraset <- MPRASet(DNA = dna, RNA = rna, eid = eid, eseq = NULL, barcode = NULL)
mpraset

#mpra data set analysis, this block of code may be expected to take ~10ish mins to run
design <- data.frame(intcpt = 1,
                     post = grepl("RPC", colnames(mpraset)),
                     right = grepl("3", colnames(mpraset)))

mpralm_fit <- mpralm(object = mpraset, design = design,
                     aggregate = "sum", normalize = TRUE,
                     model_type = "indep_groups", plot = TRUE)
toptab <- topTable(mpralm_fit, coef = 2, number = Inf)
head(toptab)

#This example code will make modifications to the output analysis dataframe. 
#You'll need to replace "gene_name_mpralm" variable with the name you've chosen for your dataframe

his4postvRPC31post_sum_mpralm <- toptab 
setDT(his4postvRPC31post_sum_mpralm, keep.rownames = TRUE)[]
names(his4postvRPC31post_sum_mpralm)[names(his4postvRPC31post_sum_mpralm) == "rn"] <- "Guide"
his4postvRPC31post_sum_mpralm <- merge(his4postvRPC31post_sum_mpralm, guide.good.targets, by="Guide")
his4postvRPC31post_sum_mpralm$gene <- sgd[match(his4postvRPC31post_sum_mpralm$Yorf1, sgd$name), "gene"]
his4postvRPC31post_sum_mpralm$desc <- sgd[match(his4postvRPC31post_sum_mpralm$Yorf1, sgd$name), "desc"]
head(his4postvRPC31post_sum_mpralm)

#_________________


#saving mpralm analysis
write.table(his4postvHTS1post_sum_mpralm, "~/CiBER_seq_package/all_raw_fastq/epistasis_3AT_HTS1_RPC31/his4postvHTS1post_sum_mpralm.txt", sep="\t")
write.table(his4postvRPC31post_sum_mpralm, "~/CiBER_seq_package/all_raw_fastq/epistasis_3AT_HTS1_RPC31/his4postvRPC31post_sum_mpralm.txt", sep="\t")





