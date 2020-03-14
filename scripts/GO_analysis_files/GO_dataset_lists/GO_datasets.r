#!/usr/bin/env Rscript

if (!requireNamespace("dplyr", quietly = TRUE))
  install.packages("dplyr")
library(dplyr)

mpralmpath <- "~/CiBER_seq_package/all_raw_fasta_gz/"
options(stringsAsFactors=FALSE)
for (dir_file in c("HIS4_PGK1_3AT/his4_postv3AT_sum_mpralm", "HIS4_PGK1_3AT/pgk1_postv3AT_sum_mpralm",
                   "HIS4_PGK1_pooled/his4_pooled_sum_mpralm", "HIS4_PGK1_pooled/pgk1_pooled_sum_mpralm",
                   "HTS1_RPC31/HTS1_sum_mpralm", "HTS1_RPC31/RPC31_sum_mpralm",
                   "GCN4_CDS_UTR/cds_sum_mpralm", "GCN4_CDS_UTR/utr_sum_mpralm")){
  y <- read.delim(paste(mpralmpath, dir_file, ".txt", sep=""),
                  stringsAsFactors=FALSE, header = TRUE)
  temp1 <- gsub(".*/", "", dir_file)
  temp3 <- temp2 <- gsub("_sum_mpralm", "", temp1)
  assign(paste(temp2, sep=""), y)
  m <- y
  m[m == "Neg_ctrl"] <- NA
  m <- m[complete.cases(m[,9]),]
  m_up <- filter(m, m$logFC < -0.5 & m$adj.P.Val < 0.05)
  m_dwn <- filter(m, m$logFC > 0.5 & m$adj.P.Val < 0.05)
  assign(paste(temp3, "sig_up", sep=""), m_up)
  assign(paste(temp3, "sig_dwn", sep=""), m_dwn)
  write.table(m$Yorf1, paste("~/CiBER_seq_package/scripts/GO_analysis_files/GO_dataset_lists/", temp3, "all.txt", sep=""), 
  row.names = FALSE, col.names = FALSE, quote = FALSE)  
  write.table(m_up$Yorf1, paste("~/CiBER_seq_package/scripts/GO_analysis_files/GO_dataset_lists/", temp3, "sig_up.txt", sep=""), 
  row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(m_dwn$Yorf1, paste("~/CiBER_seq_package/scripts/GO_analysis_files/GO_dataset_lists/", temp3, "sig_dwn.txt", sep=""), 
  row.names = FALSE, col.names = FALSE, quote = FALSE)
}

