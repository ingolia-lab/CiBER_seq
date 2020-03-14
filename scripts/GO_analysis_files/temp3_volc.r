#!/usr/bin/env Rscript

if (!requireNamespace("dplyr", quietly = TRUE))
  install.packages("dplyr")
library(dplyr)

mpralmpath <- "~/CiBER_seq_package/all_raw_fasta_gz/"
GOpath <- "~/CiBER_seq_package/scripts/GO_analysis_files/GO_annotation_lists/"
options(stringsAsFactors=FALSE)
for (annot in c("aabiosynthesis_handcur", "polIIIsubunit_handcur", "transcription_handcur", 
"translationcontrol_handcur", "tRNAaasynthetase_handcur", "tRNAprocessing_handcur", "glycolysis_handcur")){
  x <- read.table(paste(GOpath, annot, ".txt", sep=""),
                  stringsAsFactors=FALSE, header = TRUE)
  assign(paste(annot, sep=""), x)
  for (dir_file in c("HIS4_PGK1_3AT/his4_postv3AT_sum_mpralm", "HIS4_PGK1_3AT/pgk1_postv3AT_sum_mpralm",
                     "HIS4_PGK1_pooled/his4_pooled_sum_mpralm", "HIS4_PGK1_pooled/pgk1_pooled_sum_mpralm",
                     "HTS1_RPC31/HTS1_sum_mpralm", "HTS1_RPC31/RPC31_sum_mpralm",
                     "GCN4_CDS_UTR/cds_sum_mpralm", "GCN4_CDS_UTR/utr_sum_mpralm")){
    y <- read.delim(paste(mpralmpath, dir_file, ".txt", sep=""),
                    stringsAsFactors=FALSE, header = TRUE)
    temp1 <- gsub(".*/", "", dir_file)
    temp3 <- temp2 <- gsub("_sum_mpralm", "", temp1)
    assign(paste(temp2, sep=""), y)
    m <- subset(y, (y$Yorf1 %in% x$name))    
    assign(paste(annot, temp3, sep=""), m)
  }
}

his4grna <- filter(his4_pooled, his4_pooled$gene == "HIS4")

head(aabiosynthesis_handcur)
head(his4_pooled)
head(aabiosynthesis_handcurhis4_pooled)

head(transcription_handcurhis4_pooled)
#write.table(transcription_handcurhis4_pooled, "/mnt/ingolialab/rmuller1/test_sandbox/gz_to_mpralm/rehearsal_dinner/scripts/GO_analysis_files/transcription_handcurhis4_pooled.txt")

## Hybrid PDF-PNG scatter plots
## based on https://jonathanchang.org/blog/how-to-partially-rasterize-a-figure-plotted-with-r/
## by jonathan.chang@monash.edu
## We use par(mar=...) and plot(...) rather than plot.window(mar=...) and points(...)

#install.packages("png")
library(png)
library(scales)
options(stringsAsFactors=FALSE)

#His4 pooled before vs after guide tet induction treatment_________________

## Set up values to plot
his4_pooled_insig <- filter(his4_pooled, his4_pooled$adj.P.Val > 0.05)
temp1 <- pmax(-4, his4_pooled$logFC)
temp2 <- pmin(4, temp1)
his4_pooled$LFC_plot <- -temp2
his4_pooled$sig_plot <- -log10(pmax(10^-15, his4_pooled$adj.P.Val))

pdf("his4_pool_activators.pdf", useDingbats=FALSE, width=5, height=5, colormodel="rgb")
## Plot with type="n" to skip drawing the points
plot(x = his4_pooled$LFC_plot, y = his4_pooled$sig_plot,
     type="n",
     xlim=c(-4, 4), ylim=c(0, 15),
     axes=FALSE,
     xlab="Log2 Fold Before vs After Guide tet Induction", 
     ylab = expression('Statistical Significance -log10(Q)'), 
     col="gray52")
points(-pmin(4, aabiosynthesis_handcurhis4_pooled$logFC),
       -log10(pmax(10^-15, aabiosynthesis_handcurhis4_pooled$adj.P.Val)),
       type="n",
       col="royalblue1"
)
points(-pmin(4, tRNAaasynthetase_handcurhis4_pooled$logFC),
       -log10(pmax(10^-15, tRNAaasynthetase_handcurhis4_pooled$adj.P.Val)),
       type="n",
       col="seagreen4"
)
points(-pmin(4, translationcontrol_handcurhis4_pooled$logFC),
       -log10(pmax(10^-15, translationcontrol_handcurhis4_pooled$adj.P.Val)),
       type="n",
       col="paleturquoise"
)
points(-pmin(4, polIIIsubunit_handcurhis4_pooled$logFC),
       -log10(pmax(10^-15, polIIIsubunit_handcurhis4_pooled$adj.P.Val)),
       type="n",
       col="violet"
)
points(-pmin(4, tRNAprocessing_handcurhis4_pooled$logFC),
       -log10(pmax(10^-15, tRNAprocessing_handcurhis4_pooled$adj.P.Val)),
       type="n",
       col="red4"
)
points(-pmin(4, his4grna$logFC),
       -log10(pmax(10^-15, his4grna$adj.P.Val)),
       type="n",
       col="#b58900"
)
#rect(-6,0,6,-log10(0.05), col = alpha("gray52", 0.3))
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
pointsfile <- "his4_pool_activators.png"
png(pointsfile, width=width, height=height, units="in", res=1500, bg="transparent")
par(mar=c(0,0,0,0))
## Exactly the same plot as above, but no type="n" and instead specify points, color, etc.
plot(x = his4_pooled$LFC_plot, y = his4_pooled$sig_plot,
     pch=16, cex=0.9,
     xlim=c(-4, 4), ylim=c(0, 15),
     axes=FALSE,
     xlab="Log2 Fold Before vs After Guide tet Induction",
     ylab = expression('Statistical Significance -log10(Q)'),
     col="gray52")
points(-pmin(4, aabiosynthesis_handcurhis4_pooled$logFC),
       -log10(pmax(10^-15, aabiosynthesis_handcurhis4_pooled$adj.P.Val)),
       pch=16, cex=0.9,
       col="royalblue1"
)
points(-pmin(4, tRNAaasynthetase_handcurhis4_pooled$logFC),
       -log10(pmax(10^-15, tRNAaasynthetase_handcurhis4_pooled$adj.P.Val)),
       pch=16, cex=0.9,
       col="seagreen4"
)
points(-pmin(4, translationcontrol_handcurhis4_pooled$logFC),
       -log10(pmax(10^-15, translationcontrol_handcurhis4_pooled$adj.P.Val)),
       pch=16, cex=0.9,
       col="paleturquoise"
)
points(-pmin(4, polIIIsubunit_handcurhis4_pooled$logFC),
       -log10(pmax(10^-15, polIIIsubunit_handcurhis4_pooled$adj.P.Val)),
       pch=16, cex=0.9,
       col="violet"
)
points(-pmin(4, tRNAprocessing_handcurhis4_pooled$logFC),
       -log10(pmax(10^-15, tRNAprocessing_handcurhis4_pooled$adj.P.Val)),
       pch=16, cex=0.9,
       col="red4"
)
points(-pmin(4, his4grna$logFC),
       -log10(pmax(10^-15, his4grna$adj.P.Val)),
       pch=16, cex=0.9,
       col="#b58900"
)
#rect(-6,0,6,-log10(0.05), col = alpha("gray52", 0.3))
abline(h = 1.3, col="black")
abline(v = -0.5, col="black")
abline(v = 0.5, col="black")
dev.off()

## Read in the PNG file and draw as a raster image
panel <- readPNG(pointsfile, native = TRUE)
rasterImage(panel, coords[1], coords[3], coords[2], coords[4])

## Draw axes, title, etc.
axis(side = 1, at = c(-4, -2, 0, 2, 4), labels=c("<1/16", "1/4", "1", "4", ">16"))
axis(side = 2, at = c(0, 5, 10, 15), labels=c("0", "5", "10", ">15"))
title(main="P(His4) Guide Induction")
dev.off()
