#!/usr/bin/env Rscript

datadir = "~/CiBER_seq_package/all_raw_fastq/bc_validation/"

options(stringsAsFactors=FALSE)
if (!exists("DNA_RNA_BC1")) {
  DNA_RNA_BC1 <- read.delim(paste(datadir, "DNA_RNA.BC1.tsv", sep = ""), stringsAsFactors=FALSE, header = FALSE)
}

if (!exists("DNA_RNA_BC2")) {
  DNA_RNA_BC2 <- read.delim(paste(datadir, "DNA_RNA.BC2.tsv", sep = ""), stringsAsFactors=FALSE, header = FALSE)
}

if (!exists("DNA_RNA_BC3")) {
  DNA_RNA_BC3 <- read.delim(paste(datadir, "DNA_RNA.BC3.tsv", sep = ""), stringsAsFactors=FALSE, header = FALSE)
}

if (!exists("DNA_RNA_BC4")) {
  DNA_RNA_BC4 <- read.delim(paste(datadir, "DNA_RNA.BC4.tsv", sep = ""), stringsAsFactors=FALSE, header = FALSE)
}

if (!exists("DNA_RNA_BC5")) {
  DNA_RNA_BC5 <- read.delim(paste(datadir, "DNA_RNA.BC5.tsv", sep = ""), stringsAsFactors=FALSE, header = FALSE)
}

if (!exists("DNA_RNA_BC6")) {
  DNA_RNA_BC6 <- read.delim(paste(datadir, "DNA_RNA.BC6.tsv", sep = ""), stringsAsFactors=FALSE, header = FALSE)
}

head(DNA_RNA_BC1)

#loading packages and libraries
if (!requireNamespace("dplyr", quietly = TRUE))
  install.packages("dplyr")

if (!requireNamespace("data.table", quietly = TRUE))
  install.packages("data.table")

library(dplyr)
library(data.table)


#filter for either left or right > 32 read counts and merge with barcode assignment dataframe
DNA_RNA_BC1_32 <- filter(DNA_RNA_BC1,
                  DNA_RNA_BC1$V2 > 32 & DNA_RNA_BC1$V3 > 32)

DNA_RNA_BC2_32 <- filter(DNA_RNA_BC2,
                  DNA_RNA_BC2$V2 > 32 & DNA_RNA_BC2$V3 > 32)

DNA_RNA_BC3_32 <- filter(DNA_RNA_BC3,
                  DNA_RNA_BC3$V2 > 32 & DNA_RNA_BC3$V3 > 32)

DNA_RNA_BC4_32 <- filter(DNA_RNA_BC4,
                  DNA_RNA_BC4$V2 > 32 & DNA_RNA_BC4$V3 > 32)

DNA_RNA_BC5_32 <- filter(DNA_RNA_BC5,
                  DNA_RNA_BC5$V2 > 32 & DNA_RNA_BC5$V3 > 32)

DNA_RNA_BC6_32 <- filter(DNA_RNA_BC6,
                  DNA_RNA_BC6$V2 > 32 & DNA_RNA_BC6$V3 > 32)

head(DNA_RNA_BC1_32)

#calculate RNA:DNA ratios

DNA_RNA_BC1_32$RNA_DNA_ratio <- (DNA_RNA_BC1_32$V4 + DNA_RNA_BC1_32$V5)/(DNA_RNA_BC1_32$V2 + DNA_RNA_BC1_32$V3)
DNA_RNA_BC1_32$cond <- "a.0nM"
head(DNA_RNA_BC1_32)

DNA_RNA_BC2_32$RNA_DNA_ratio <- (DNA_RNA_BC2_32$V4 + DNA_RNA_BC2_32$V5)/(DNA_RNA_BC2_32$V2 + DNA_RNA_BC2_32$V3)
DNA_RNA_BC2_32$cond <- "b.16nM"
head(DNA_RNA_BC2_32)

DNA_RNA_BC3_32$RNA_DNA_ratio <- (DNA_RNA_BC3_32$V4 + DNA_RNA_BC3_32$V5)/(DNA_RNA_BC3_32$V2 + DNA_RNA_BC3_32$V3)
DNA_RNA_BC3_32$cond <- "c.128nM"
head(DNA_RNA_BC3_32)

DNA_RNA_BC4_32$RNA_DNA_ratio <- (DNA_RNA_BC4_32$V4 + DNA_RNA_BC4_32$V5)/(DNA_RNA_BC4_32$V2 + DNA_RNA_BC4_32$V3)
DNA_RNA_BC4_32$cond <- "a.No_gRNA"
head(DNA_RNA_BC4_32)

DNA_RNA_BC5_32$RNA_DNA_ratio <- (DNA_RNA_BC5_32$V4 + DNA_RNA_BC5_32$V5)/(DNA_RNA_BC5_32$V2 + DNA_RNA_BC5_32$V3)
DNA_RNA_BC5_32$cond <- "CGA_SC"
head(DNA_RNA_BC5_32)

DNA_RNA_BC6_32$RNA_DNA_ratio <- (DNA_RNA_BC6_32$V4 + DNA_RNA_BC6_32$V5)/(DNA_RNA_BC6_32$V2 + DNA_RNA_BC6_32$V3)
DNA_RNA_BC6_32$cond <- "b.ADH1_gRNA"
head(DNA_RNA_BC6_32)

est_RNA_DNA_ratios <- append(DNA_RNA_BC1_32$RNA_DNA_ratio, DNA_RNA_BC2_32$RNA_DNA_ratio)
est_RNA_DNA_ratios <- append(est_RNA_DNA_ratios, DNA_RNA_BC3_32$RNA_DNA_ratio)
est_conditions <- append(DNA_RNA_BC1_32$cond, DNA_RNA_BC2_32$cond)
est_conditions <- append(est_conditions, DNA_RNA_BC3_32$cond)

gRNA_RNA_DNA_ratios <- append(DNA_RNA_BC4_32$RNA_DNA_ratio, DNA_RNA_BC6_32$RNA_DNA_ratio)
gRNA_conditions <- append(DNA_RNA_BC4_32$cond, DNA_RNA_BC6_32$cond)

est <- data.frame(est_RNA_DNA_ratios, est_conditions)
names(est) <- c("est_RNA_DNA_ratios", "est_conditions")
gRNA <- data.frame(gRNA_RNA_DNA_ratios, gRNA_conditions)
names(gRNA) <- c("gRNA_RNA_DNA_ratios", "gRNA_conditions")

est$est_conditions <- as.factor(est$est_conditions)
gRNA$gRNA_conditions <- as.factor(gRNA$gRNA_conditions)

head(est)
tail(est)
head(gRNA)
tail(gRNA)

library(ggplot2)

l10Decades <- function(xs) {
  log10(c(sapply(xs, function(x) { seq(2,9)*(10**x) })))
}

# Violin with blox plot inside for estradiol titration
estra <- ggplot(est, aes(x=est_conditions, y=log10(est_RNA_DNA_ratios + 1), fill=est_conditions)) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +  
  geom_violin(trim = TRUE, scale = "area", width=1.1)+
  geom_boxplot(width=0.04, fill="white", outlier.shape = NA)+
  labs(title="Plot of RNA_DNA ratio by estradiol concentration",x="Estradiol (nM)", y = "log 10 RNA_DNA ratio")

pdf("est_barcode.pdf", useDingbats=FALSE, width=8, height=8, colormodel="rgb")
estra + theme_classic() + annotation_logticks(base = 10, sides = "l", scaled = TRUE,
  colour = "black", size = 0.5, linetype = 1, alpha = 1)
dev.off()

# Violin with blox plot inside for tet titration
tet <- ggplot(gRNA, aes(x=gRNA_conditions, y=log10(gRNA_RNA_DNA_ratios + 1), fill=gRNA_conditions)) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  geom_violin(trim = TRUE, scale = "area", width=1.1)+
  geom_boxplot(width=0.04, fill="white", outlier.shape = NA)+
  labs(title="Plot of RNA_DNA ratio by tet concentration",x="Anhydrotetracycline (ng/mL)", y = "log 10 RNA_DNA ratio")

pdf("tet_barcode.pdf", useDingbats=FALSE, width=8, height=8, colormodel="rgb")
tet + theme_classic() + annotation_logticks(base = 10, sides = "l", scaled = TRUE,
  colour = "black", size = 0.5, linetype = 1, alpha = 1)
dev.off()


#saving modified datasets
#write.table(cds_sum_mpralm, "~/CiBER_seq_package/all_raw_fasta_gz/GCN4_CDS_UTR/cds_sum_mpralm.txt", sep="\t")
#write.table(utr_sum_mpralm, "~/CiBER_seq_package/all_raw_fasta_gz/GCN4_CDS_UTR/utr_sum_mpralm.txt", sep="\t")

print(length(DNA_RNA_BC1_32$RNA_DNA_ratio))
print(length(DNA_RNA_BC2_32$RNA_DNA_ratio))
print(length(DNA_RNA_BC3_32$RNA_DNA_ratio))
print(length(DNA_RNA_BC4_32$RNA_DNA_ratio))
print(length(DNA_RNA_BC6_32$RNA_DNA_ratio))
