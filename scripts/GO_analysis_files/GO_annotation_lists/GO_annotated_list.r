#!/usr/bin/env Rscript

#creating hand_curated annotated gene lists

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

actin_extra <- c("ACT1", "COF1", "ARP2", "ARP3", "ARC35", "ARC40", "ARC19", "ARC18", 
"ARC15", "SCD5", "CAP1", "ARP8", "YSC84", "ARK1", "ENT1", "GEA1", 
"YIH1", "ARP9", "AIM7", "ACF2", "ARP6", "SYP1")

ER_traff <- c("SEC65", "SRP14", "SRP21", "SRP54", "SRP68", 
         "SRP72", "SRP102", "SEC62", "SEC61", "SRP101", "PHO88", "SEC72", "SEC63", "SEC66")

mitochondrial <- filter(his4_pooled, grepl("mitochond",his4_pooled$desc))
mitochondrial <- unique(mitochondrial$gene)

ribosomal_subs <- c("ASC1", "RPS0A", "RPS0B", "RPS0A", "RPS1A", "RPS1B", "RPS2", 
"RPS3", "RPS4A", "RPS4B", "RPS5", "RPS6A", "RPS6B", "RPS7A", "RPS7B", "RPS8A", "RPS8B", 
"RPS9A", "RPS9B", "RPS10A", "RPS10B", "RPS11A", "RPS11B", "RPS12", "RPS13", "RPS14A", 
"RPS14B", "RPS15", "RPS16A", "RPS16B", "RPS17A", "RPS17B", "RPS18A", "RPS18B", 
"RPS19A", "RPS19B", "RPS20", "RPS21A", "RPS21B", "RPS22A", "RPL2A", "RPL2B", "RPL3", 
"RPL4A", "RPL4B", "RPL5", "RPL6A", "RPL6B", "RPL7A", "RPL7B", "RPL8A", "RPL8B", 
"RPL9A", "RPL9B", "RPL10", "RPL11A", "RPL11B", "RPL12A", "RPL12B", "RPL13A", 
"RPL13B", "RPL14A", "RPL14B", "RPL15A", "RPL15B", "RPL16A", "RPL16B", "RPL17A", 
"RPL17B", "RPL18A", "RPL18B", "RPL19A", "RPL19B", "RPL20A", "RPL20B", "RPL21A", "RPL21B", 
"RPL22A", "RPL22B")


actin_extra_handcur <- subset(sgd, (sgd$gene %in% actin_extra))
srp_handcur <- subset(sgd, (sgd$gene %in% srp))
ribosomal_subs_handcur <- subset(sgd, (sgd$gene %in% ribosomal_subs))

ER_traff_handcur <- subset(sgd, (sgd$gene %in% ER_traff))

mitochondrial_handcur <- filter(sgd, grepl("mitochond",sgd$desc))
mitochondrial_handcur <- filter(mitochondrial_handcur, !grepl("Dubious",mitochondrial_handcur$qual))
mitochondrial_handcur <- filter(mitochondrial_handcur, !grepl("Uncharacterized",mitochondrial_handcur$qual))

write.table(actin_extra_handcur, 
            "~/CiBER_seq_package/scripts/GO_analysis_files/GO_annotation_lists/actin_extra_handcur.txt", 
            sep="\t")
write.table(srp_handcur, 
            "~/CiBER_seq_package/scripts/GO_analysis_files/GO_annotation_lists/srp_handcur.txt", 
            sep="\t")
write.table(ribosomal_subs_handcur, 
            "~/CiBER_seq_package/scripts/GO_analysis_files/GO_annotation_lists/ribosomal_subs_handcur.txt", 
            sep="\t")

write.table(ER_traff_handcur, 
            "~/CiBER_seq_package/scripts/GO_analysis_files/GO_annotation_lists/ER_traff_handcur.txt", 
            sep="\t")
write.table(mitochondrial_handcur, 
            "~/CiBER_seq_package/scripts/GO_analysis_files/GO_annotation_lists/mitochondrial_handcur.txt", 
            sep="\t")

#Gene Lists
#RPN3, RPN7, RPN15, RPN12, RPN11, RPN5, RPN6, RPN9, RPN8, RPN10, RPN1, RPN2, RPN13
#RPT1, RPT2, RPT3, RPT4, RPT5, RPT6, 
#PRE1, PRE2, PRE3, PRE4, PRE5, PRE6, PRE7, PRE8, PRE9, PRE10,
#PUP1, PUP2, PUP3,
#SCL1, UMP1

#TCP1, CCT2, CCT3, CCT4, CCT5, CCT6, CCT7, CCT8

#ACT1, COF1, ARP2, ARP3, ARC35, ARC40, ARC19, ARC18, ARC15


