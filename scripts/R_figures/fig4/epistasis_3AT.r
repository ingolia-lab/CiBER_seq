#!/usr/bin/env Rscript

#3AT epistasis analysis
# In the first phase of our experiment, we separate guides that increase His4 transcription (pre vs post)
# In the second phase, we identify guides that increase activation further 
#(knockdown makes cells hypersensitive to 3AT)
# as well as guides that prevent additional activation (either have already activated or are epistatic)

if (!requireNamespace("dplyr", quietly = TRUE))
  install.packages("dplyr")
library(dplyr)

if (!file.exists("SGD_features.tab")) {
  sgd <- download.file('https://downloads.yeastgenome.org/curation/chromosomal_feature/SGD_features.tab', destfile="SGD_features.tab")
}
sgd <- read.delim("SGD_features.tab", header=FALSE, quote="",
                  col.names=c("sgdid", "type", "qual", "name", "gene", "alias",
                              "parent", "sgdid2", "chrom", "start", "end",
                              "strand", "genpos", "cver", "sver", "desc"))
sgd_sub <- sgd[c(4,5,16)]
names(sgd_sub) <- c("name", "gene", "desc")

mpralmpath <- "~/CiBER_seq_package/all_raw_fasta_gz/"
GOpath <- "~/CiBER_seq_package/scripts/GO_analysis_files/GO_annotation_lists/"
options(stringsAsFactors=FALSE)
for (annot in c("aabiosynthesis_handcur", "polIIIsubunit_handcur", "transcription_handcur", 
                "translationcontrol_handcur", "tRNAaasynthetase_handcur", "tRNAprocessing_handcur", "glycolysis_handcur",
                "TRiC_CCT_handcur", "proteasome_handcur", "actin_handcur")){
  x <- read.table(paste(GOpath, annot, ".txt", sep=""),
                  stringsAsFactors=FALSE, header = TRUE)
  assign(paste(annot, sep=""), x)
  for (dir_file in c("HTS1_RPC31/HTS1_sum_mpralm", "HTS1_RPC31/RPC31_sum_mpralm", 
                     "HIS4_PGK1_pooled/his4_pooled_sum_mpralm", "HIS4_PGK1_3AT/his4_postv3AT_sum_mpralm", 
                     "HIS4_PGK1_3AT/his4_prev3AT_sum_mpralm")){
    y <- read.delim(paste(mpralmpath, dir_file, ".txt", sep=""),
                    stringsAsFactors=FALSE, header = TRUE)
    temp1 <- gsub(".*/", "", dir_file)
    temp3 <- temp2 <- gsub("_sum_mpralm", "", temp1)
    assign(paste(temp2, sep=""), y)
    m <- subset(y, (y$Yorf1 %in% x$name))    
    assign(paste(annot, temp3, sep=""), m)
  }
}
  

# To filter activators over the rest
his4_pooled_sig_act <- filter(his4_pooled, his4_pooled$adj.P.Val < 0.05 & his4_pooled$logFC < 0)
his4_pooled_nonact <- filter(his4_pooled, his4_pooled$adj.P.Val > 0.05 | his4_pooled$logFC > 0)
  
# To filter those that already induce from those that block induction
his4_postv3AT_lowind <- filter(his4_postv3AT, his4_postv3AT$adj.P.Val < 0.05 & his4_postv3AT$logFC > 0)




#treat the His4 pre vs post as the stage 1 condition
#treat postv3AT, His4 post vs HTS1(post), and RPC31(prevpost) as stage 2



