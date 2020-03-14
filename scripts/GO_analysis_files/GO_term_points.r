#!/usr/bin/env Rscript

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

mpralmpath <- "/mnt/ingolialab/rmuller1/test_sandbox/gz_to_mpralm/rehearsal_dinner/all_raw_fasta_gz/"
GOpath <- "/mnt/ingolialab/rmuller1/test_sandbox/gz_to_mpralm/rehearsal_dinner/scripts/GO_analysis_files/GO_annotation_lists/"
sgd_yorf <- data.frame(sgd$sgdid, sgd$name)
head(sgd_yorf)
options(stringsAsFactors=FALSE)
for (annot in c("cellular_amino_acid_metabolic_process", "nucleic_acid_metabolic_process",
                "RNA_processing", "tRNA_transcription", "cellular_component_biogenesis",
                "proteasomal_ubiquitin-independent_protein_catabolic_process",
                "RNP_complex_biogenesis", "vesicle_fusion_with_Golgi_apparatus", "gene_expression",
                "ribosome_biogenesis", "tRNA_aminoacylation", "glycolytic_process",
                "RNA_metabolic_process", "tRNA_metabolic_process", "nuclear_pore_distribution",
                "protein_sumoylation", "regulation_of_translation")){
  x <- read.delim(paste(GOpath, annot, ".txt", sep=""),
                  stringsAsFactors=FALSE, header = FALSE)
  x$sgd.sgdid <- gsub("SGD:", "", x$V1)
  z <- merge(x, sgd_yorf, by = "sgd.sgdid")
  colnames(z)[colnames(z)=="sgd.name"] <- "Yorf1"
  assign(paste(annot, sep=""), z)
  for (dir_file in c("HIS4_PGK1_3AT/his4_postv3AT_sum_mpralm", "HIS4_PGK1_3AT/pgk1_postv3AT_sum_mpralm",
                     "HIS4_PGK1_pooled/his4_pooled_sum_mpralm", "HIS4_PGK1_pooled/pgk1_pooled_sum_mpralm",
                     "HTS1_RPC31/HTS1_sum_mpralm", "HTS1_RPC31/RPC31_sum_mpralm",
                     "GCN4_CDS_UTR/cds_sum_mpralm", "GCN4_CDS_UTR/utr_sum_mpralm")){
    y <- read.delim(paste(mpralmpath, dir_file, ".txt", sep=""),
                    stringsAsFactors=FALSE, header = TRUE)
    temp1 <- gsub(".*/", "", dir_file)
    temp3 <- temp2 <- gsub("_sum_mpralm", "", temp1)
    assign(paste(temp2, sep=""), y)
    m <- merge(z, y, by = "Yorf1")
    assign(paste(annot, temp3, sep=""), m)
  }
}

## Hybrid PDF-PNG scatter plots
## based on https://jonathanchang.org/blog/how-to-partially-rasterize-a-figure-plotted-with-r/
## by jonathan.chang@monash.edu
## We use par(mar=...) and plot(...) rather than plot.window(mar=...) and points(...)

#install.packages("png")
library(png)
library(scales)
options(stringsAsFactors=FALSE)

#His4 post vs 3AT treatment_________________

## Set up values to plot
his4_postv3AT_insig <- filter(his4_postv3AT, his4_postv3AT$adj.P.Val > 0.05)
his4_postv3AT$LFC_plot <- -his4_postv3AT$logFC
his4_postv3AT$sig_plot <- -log10(pmax(10^-20, his4_postv3AT$adj.P.Val))

pdf("his4_postv3AT_volcano.pdf", useDingbats=FALSE, width=5, height=5, colormodel="rgb")
## Plot with type="n" to skip drawing the points
plot(x = his4_postv3AT$LFC_plot, y = his4_postv3AT$sig_plot,
     type="n",
     xlim=c(-6, 6), ylim=c(0, 20),
     axes=FALSE,
     xlab="Log2 Fold Before vs After 3AT Treatment", 
     ylab = expression('Statistical Significance -log10(Q)'), 
     col="gray52")
points(-cellular_amino_acid_metabolic_processhis4_postv3AT$logFC,
       -log10(pmax(10^-20, cellular_amino_acid_metabolic_processhis4_postv3AT$adj.P.Val)),
       type="n",
       col="violet"
)
points(-tRNA_transcriptionhis4_postv3AT$logFC,
       -log10(pmax(10^-20, tRNA_transcriptionhis4_postv3AT$adj.P.Val)),
       type="n",
       col="palegreen3"
)
points(-tRNA_metabolic_processhis4_postv3AT$logFC,
       -log10(pmax(10^-20, tRNA_metabolic_processhis4_postv3AT$adj.P.Val)),
       type="n",
       col="orange1"
)
points(-tRNA_aminoacylationhis4_postv3AT$logFC,
       -log10(pmax(10^-20, tRNA_aminoacylationhis4_postv3AT$adj.P.Val)),
       type="n",
       col="orangered1"
)
rect(-6,0,6,-log10(0.05), col = alpha("gray52", 0.3))
abline(h = 1.3, col="black")
abline(v = -0.5, col="black")
abline(v = 0.5, col="black")

#rect(-6,0,6,-log10(0.05), col=rgb(0.5,0.5,0.5,1/4, alpha=0.3))

## Capture the plot window coordinates
coords <- par("usr")
gx <- grconvertX(coords[1:2], "user", "inches")
gy <- grconvertY(coords[3:4], "user", "inches")
width <- max(gx) - min(gx)
height <- max(gy) - min(gy)

## Plot points to a PNG file
pointsfile <- "his4_postv3AT_volcano.png"
png(pointsfile, width=width, height=height, units="in", res=1500, bg="transparent")
par(mar=c(0,0,0,0))
## Exactly the same plot as above, but no type="n" and instead specify points, color, etc.
plot(x = his4_postv3AT$LFC_plot, y = his4_postv3AT$sig_plot,
     pch=16, cex=0.9,
     xlim=c(-6, 6), ylim=c(0, 20),
     axes=FALSE,
     xlab="Log2 Fold Before vs After 3AT Treatment", 
     ylab = expression('Statistical Significance -log10(Q)'), 
     col = alpha("gray52", 1.0))
points(-cellular_amino_acid_metabolic_processhis4_postv3AT$logFC,
       -log10(pmax(10^-20, cellular_amino_acid_metabolic_processhis4_postv3AT$adj.P.Val)),
       pch=16, cex=0.9,
       col = alpha("mediumorchid2", 1.0)
)
points(-tRNA_transcriptionhis4_postv3AT$logFC,
       -log10(pmax(10^-20, tRNA_transcriptionhis4_postv3AT$adj.P.Val)),
       pch=16, cex=0.9,
       col = alpha("palegreen4", 1.0)
)
points(-tRNA_metabolic_processhis4_postv3AT$logFC,
       -log10(pmax(10^-20, tRNA_metabolic_processhis4_postv3AT$adj.P.Val)),
       pch=16, cex=0.9,
       col = alpha("orange1", 1.0)
)
points(-tRNA_aminoacylationhis4_postv3AT$logFC,
       -log10(pmax(10^-20, tRNA_aminoacylationhis4_postv3AT$adj.P.Val)),
       pch=16, cex=0.9,
       col = alpha("orangered1", 1.0)
)
rect(-6,0,6,-log10(0.05), col = alpha("gray52", 0.3))
abline(h = 1.3, col="black")
abline(v = -0.5, col="black")
abline(v = 0.5, col="black")
dev.off()

## Read in the PNG file and draw as a raster image
panel <- readPNG(pointsfile, native = TRUE)
rasterImage(panel, coords[1], coords[3], coords[2], coords[4])

## Draw axes, title, etc.
axis(side = 1, at = c(-6, -4, -2, 0, 2, 4, 6), labels=c("1/64", "1/16", "1/4", "1", "4", "16", "64"))
axis(side = 2, at = c(0, 5, 10, 15, 20), labels=c("0", "5", "10", "15", ">20"), col = alpha("orangered1", 1.0))
title(main="P(His4) 3AT Treatment")
dev.off()
