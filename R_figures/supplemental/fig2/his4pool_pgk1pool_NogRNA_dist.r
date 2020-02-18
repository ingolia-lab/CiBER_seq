#!/usr/bin/env Rscript

if (!requireNamespace("dplyr", quietly = TRUE))
  install.packages("dplyr")
library(dplyr)

mpralmpath <- "~/CiBER_seq_package/all_raw_fasta_gz/"
GOpath <- "~/CiBER_seq_package/scripts/GO_analysis_files/GO_annotation_lists/"

pgk1_pooled <- read.delim(paste(mpralmpath, "HIS4_PGK1_pooled/pgk1_pooled_sum_mpralm", ".txt", sep=""),
                          stringsAsFactors=FALSE, header = TRUE)
his4_pooled <- read.delim(paste(mpralmpath, "HIS4_PGK1_pooled/his4_pooled_sum_mpralm", ".txt", sep=""),
                          stringsAsFactors=FALSE, header = TRUE)
pgk1_pooled_sub <- pgk1_pooled[c(1,2,6,9)] 
his4_pooled_sub <- his4_pooled[c(1,2,6,9)] 
names(pgk1_pooled_sub) <- c("Guide", "pgk1_pooled_logFC", "pgk1_pooled_adj.P.val", "pgk1_pooled_Yorf1")
names(his4_pooled_sub) <- c("Guide", "his4_pooled_logFC", "his4_pooled_adj.P.val", "his4_pooled_Yorf1")

his4no_grna <- filter(his4_pooled, grepl("No_gRNA+", his4_pooled_sub$Guide))
pgk1no_grna <- filter(pgk1_pooled, grepl("No_gRNA+", pgk1_pooled_sub$Guide))

hist(his4no_grna$P.Value, breaks = 100)
hist(pgk1no_grna$P.Value, breaks = 100)

hist(his4no_grna$adj.P.Val)


