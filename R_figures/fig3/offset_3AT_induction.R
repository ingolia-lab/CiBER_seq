#offset 3AT treatment
no_gRNA <- filter(his4_postv3AT, grepl("No_gRNA",his4_postv3AT$Guide))
median(no_gRNA$logFC)
mean(no_gRNA$logFC)

options(stringsAsFactors=FALSE)
if (!exists("pgk1_postv3AT")) {
  pgk1_postv3AT <- read.delim("~/CiBER_seq_package/all_raw_fasta_gz/HIS4_PGK1_3AT/pgk1_postv3AT_sum_mpralm.txt",
                                   stringsAsFactors=FALSE)
}

no_gRNA_pgk1 <- filter(pgk1_postv3AT, grepl("No_gRNA",pgk1_postv3AT$Guide))
median(no_gRNA_pgk1$logFC)
mean(no_gRNA_pgk1$logFC)

-0.3975933 -0.02257629



