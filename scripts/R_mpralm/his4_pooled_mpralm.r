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
# options(stringsAsFactors=FALSE)
# if (!exists("his4_seq1")) {
#   his4_seq1 <- read.delim("~/CiBER_seq_package/all_raw_fastq/HIS4_PGK1_3AT/all_his4_seq1_counts.txt",
#                                 stringsAsFactors=FALSE)
# }
# 
# options(stringsAsFactors=FALSE)
# if (!exists("pgk1_seq1")) {
#   pgk1_seq1 <- read.delim("~/CiBER_seq_package/all_raw_fastq/HIS4_PGK1_3AT/all_pgk1_seq1_counts.txt",
#                                 stringsAsFactors=FALSE)
# }

# Load in raw count data and relevant gRNA barcode assignment and gRNA assignment of target files
options(stringsAsFactors=FALSE)
if (!exists("his4_moreseq")) {
  all_his4 <- read.delim("~/CiBER_seq_package/all_raw_fastq/HIS4_PGK1_pooled/all_his4_pooled_counts.txt",
                                stringsAsFactors=FALSE)
}

options(stringsAsFactors=FALSE)
if (!exists("pgk1_moreseq")) {
  all_pgk1 <- read.delim("~/CiBER_seq_package/all_raw_fastq/HIS4_PGK1_pooled/all_pgk1_pooled_counts.txt",
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


grna.assign.barcode.grna.good$guide <- as.character(grna.assign.barcode.grna.good$guide)
grna.assign.barcode.grna.good$guide[grna.assign.barcode.grna.good$guide == "No_gRNA"] <- paste("No_gRNA", seq(1:788), sep=" ")

# #Incorporate a lookup table strategy to assign long column names to easier variale names
# old_pgk1_vars <- c("barcode", "IVT_3AT_L_S21_L008_R1_001pgk1.count.txt", "IVT_3AT_R_S22_L008_R1_001pgk1.count.txt", "IVT_PostL_S19_L008_R1_001pgk1.count.txt",
#                    "IVT_PostR_S20_L008_R1_001pgk1.count.txt", "IVT_PreL_S17_L008_R1_001pgk1.count.txt", "IVT_PreR_S18_L008_R1_001pgk1.count.txt",
#                    "RNA_3AT_L_S27_L008_R1_001pgk1.count.txt", "RNA_3AT_R_S28_L008_R1_001pgk1.count.txt", "RNA_PostL_S25_L008_R1_001pgk1.count.txt",
#                    "RNA_PostR_S26_L008_R1_001pgk1.count.txt", "RNA_PreL_S23_L008_R1_001pgk1.count.txt", "RNA_PreR_S24_L008_R1_001pgk1.count.txt")
# 
# old_his4_vars <- c("barcode", "IVT_3AT_L_S21_L008_R1_001his4.count.txt", "IVT_3AT_R_S22_L008_R1_001his4.count.txt", "IVT_PostL_S19_L008_R1_001his4.count.txt",
#                    "IVT_PostR_S20_L008_R1_001his4.count.txt", "IVT_PreL_S17_L008_R1_001his4.count.txt", "IVT_PreR_S18_L008_R1_001his4.count.txt",
#                    "RNA_3AT_L_S27_L008_R1_001his4.count.txt", "RNA_3AT_R_S28_L008_R1_001his4.count.txt", "RNA_PostL_S25_L008_R1_001his4.count.txt",
#                    "RNA_PostR_S26_L008_R1_001his4.count.txt", "RNA_PreL_S23_L008_R1_001his4.count.txt", "RNA_PreR_S24_L008_R1_001his4.count.txt")
# 
# newvars <- c("barcode", "IVT_3AT_L", "IVT_3AT_R", "IVT_postL", "IVT_postR", "IVT_preL", "IVT_preR", "RNA_3AT_L", "RNA_3AT_R", "RNA_postL", "RNA_postR", "RNA_preL", "RNA_preR")
# 
# old_pgk1_moreseq <- c("barcode", "IVT_PostL_S52_L008_R1_001pgk1.count.txt", "IVT_PostR_S53_L008_R1_001pgk1.count.txt",
#                       "IVT_PreL_S50_L008_R1_001pgk1.count.txt", "IVT_PreR_S51_L008_R1_001pgk1.count.txt", "RNA_PostL_S56_L008_R1_001pgk1.count.txt",
#                       "RNA_PostR_S57_L008_R1_001pgk1.count.txt", "RNA_PreL_S54_L008_R1_001pgk1.count.txt",
#                       "RNA_PreR_S55_L008_R1_001pgk1.count.txt")
# 
# old_his4_moreseq <- c("barcode", "IVT_PostL_S52_L008_R1_001his4.count.txt", "IVT_PostR_S53_L008_R1_001his4.count.txt",
#                       "IVT_PreL_S50_L008_R1_001his4.count.txt", "IVT_PreR_S51_L008_R1_001his4.count.txt", "RNA_PostL_S56_L008_R1_001his4.count.txt",
#                       "RNA_PostR_S57_L008_R1_001his4.count.txt", "RNA_PreL_S54_L008_R1_001his4.count.txt",
#                       "RNA_PreR_S55_L008_R1_001his4.count.txt")
# 
# newvars2 <- c("barcode", "IVT_postL2", "IVT_postR2", "IVT_preL2", "IVT_preR2", "RNA_postL2", "RNA_postR2", "RNA_preL2", "RNA_preR2")
# 
# lookup1 = data.frame(old_pgk1_vars, old_his4_vars, newvars)
# lookup2 = data.frame(old_pgk1_moreseq, old_his4_moreseq, newvars2)
# 
# names(his4_seq1) <- lookup1[match(names(his4_seq1), lookup1$old_his4_vars),"newvars"]
# names(pgk1_seq1) <- lookup1[match(names(pgk1_seq1), lookup1$old_pgk1_vars),"newvars"]
# 
# names(his4_moreseq) <- lookup2[match(names(his4_moreseq), lookup2$old_his4_moreseq),"newvars2"]
# names(pgk1_moreseq) <- lookup2[match(names(pgk1_moreseq), lookup2$old_pgk1_moreseq),"newvars2"]
# 
# head(his4_seq1)
# head(pgk1_seq1)
# head(his4_moreseq)
# head(pgk1_moreseq)

#create dataframes that mpralm can use, pool seq runs, filter with a 32 count cut off for DNA pre samples
#filter for either left or right > 32 read counts and merge with barcode assignment dataframe
# pgk1 <- merge(pgk1_seq1, pgk1_moreseq, by="barcode", all=TRUE)
# his4 <- merge(his4_seq1, his4_moreseq, by="barcode", all=TRUE)
# pgk1[is.na(pgk1)] <- 0
# his4[is.na(his4)] <- 0
# head(pgk1)
# head(his4)
# 
# all_pgk1 <- data.frame(pgk1$barcode)
# all_his4 <- data.frame(his4$barcode)
# colnames(all_pgk1)[colnames(all_pgk1)=="pgk1.barcode"] <- "barcode"
# colnames(all_his4)[colnames(all_his4)=="his4.barcode"] <- "barcode"
# 
# all_pgk1$IVT_preL <- pgk1$IVT_preL + pgk1$IVT_preL2
# all_pgk1$IVT_postL <- pgk1$IVT_postL + pgk1$IVT_postL2
# all_pgk1$IVT_preR <- pgk1$IVT_preR + pgk1$IVT_preR2
# all_pgk1$IVT_postR <- pgk1$IVT_postR + pgk1$IVT_postR2
# all_pgk1$RNA_preL <- pgk1$RNA_preL + pgk1$RNA_preL2
# all_pgk1$RNA_postL <- pgk1$RNA_postL + pgk1$RNA_postL2
# all_pgk1$RNA_preR <- pgk1$RNA_preR + pgk1$RNA_preR2
# all_pgk1$RNA_postR <- pgk1$RNA_postR + pgk1$RNA_postR2
# 
# all_his4$IVT_preL <- his4$IVT_preL + his4$IVT_preL2
# all_his4$IVT_postL <- his4$IVT_postL + his4$IVT_postL2
# all_his4$IVT_preR <- his4$IVT_preR + his4$IVT_preR2
# all_his4$IVT_postR <- his4$IVT_postR + his4$IVT_postR2
# all_his4$RNA_preL <- his4$RNA_preL + his4$RNA_preL2
# all_his4$RNA_postL <- his4$RNA_postL + his4$RNA_postL2
# all_his4$RNA_preR <- his4$RNA_preR + his4$RNA_preR2
# all_his4$RNA_postR <- his4$RNA_postR + his4$RNA_postR2

names(all_his4) <- gsub(x = names(all_his4), pattern = "his4.count.txt", replacement = "")
names(all_his4) <- gsub(x = names(all_his4), pattern = "PH_", replacement = "")

names(all_pgk1) <- gsub(x = names(all_pgk1), pattern = "pgk1.count.txt", replacement = "")
names(all_pgk1) <- gsub(x = names(all_pgk1), pattern = "PH_", replacement = "")

head(all_pgk1)
head(all_his4)

pgk1_32 <- filter(all_pgk1, 
                  all_pgk1$IVT_preL > 32 | all_pgk1$IVT_preR > 32)

his4_32 <- filter(all_his4,
                  all_his4$IVT_preL > 32 | all_his4$IVT_preR > 32)

pgk1_32_guides <- merge(pgk1_32, grna.assign.barcode.grna.good, by="barcode")

his4_32_guides <- merge(his4_32, grna.assign.barcode.grna.good, by="barcode")

#Preparing the DNA, RNA, and Element ID (eid) dataframes for His4 mpralm analysis
dna <- his4_32_guides
dna <-data.frame(row.names=dna$barcode, dna$IVT_preL, 
                 dna$IVT_preR, 
                 dna$IVT_postL, 
                 dna$IVT_postR)
names(dna) <- c("preL", "preR", "postL", "postR")

rna <- his4_32_guides
rna <-data.frame(row.names=rna$barcode, rna$RNA_preL, 
                 rna$RNA_preR, 
                 rna$RNA_postL, 
                 rna$RNA_postR)
names(rna) <- c("preL", "preR", "postL", "postR")

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
dna <-data.frame(row.names=dna$barcode, dna$IVT_preL, 
                 dna$IVT_preR, 
                 dna$IVT_postL, 
                 dna$IVT_postR)
names(dna) <- c("preL", "preR", "postL", "postR")

rna <- pgk1_32_guides
rna <-data.frame(row.names=rna$barcode, rna$RNA_preL, 
                 rna$RNA_preR, 
                 rna$RNA_postL, 
                 rna$RNA_postR)
names(rna) <- c("preL", "preR", "postL", "postR")

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
write.table(pgk1_sum_mpralm, "~/CiBER_seq_package/all_raw_fastq/HIS4_PGK1_pooled/pgk1_pooled_sum_mpralm.txt", sep="\t")
write.table(his4_sum_mpralm, "~/CiBER_seq_package/all_raw_fastq/HIS4_PGK1_pooled/his4_pooled_sum_mpralm.txt", sep="\t")
