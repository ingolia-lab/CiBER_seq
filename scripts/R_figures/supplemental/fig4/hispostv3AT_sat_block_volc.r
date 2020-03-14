#!/usr/bin/env Rscript

if (!requireNamespace("dplyr", quietly = TRUE))
  install.packages("dplyr")
library(dplyr)

mpralmpath <- "~/CiBER_seq_package/all_raw_fasta_gz/"
GOpath <- "~/CiBER_seq_package/scripts/GO_analysis_files/GO_annotation_lists/"
options(stringsAsFactors=FALSE)
for (annot in c("aabiosynthesis_handcur", "polIIIsubunit_handcur", "transcription_handcur", 
                "translationcontrol_handcur", "tRNAaasynthetase_handcur", "tRNAprocessing_handcur", 
                "actin_handcur", "proteasome_handcur", "ribosomal_subs_handcur", 
                "ER_traff_handcur", "mitochondrial_handcur")){
  x <- read.table(paste(GOpath, annot, ".txt", sep=""),
                  stringsAsFactors=FALSE, header = TRUE)
  assign(paste(annot, sep=""), x)
  for (dir_file in c("HIS4_PGK1_3AT/his4_postv3AT_sum_mpralm")){
    y <- read.delim(paste(mpralmpath, dir_file, ".txt", sep=""),
                    stringsAsFactors=FALSE, header = TRUE)
    temp1 <- gsub(".*/", "", dir_file)
    temp3 <- temp2 <- gsub("_sum_mpralm", "", temp1)
    assign(paste(temp2, sep=""), y)
    m <- subset(y, (y$Yorf1 %in% x$name))    
    assign(paste(annot, temp3, sep=""), m)
  }
}

his4grna <- filter(his4_postv3AT, his4_postv3AT$gene == "HIS4")

head(aabiosynthesis_handcur)
head(his4_postv3AT)
head(aabiosynthesis_handcurhis4_postv3AT)

head(transcription_handcurhis4_postv3AT)
#write.table(transcription_handcurhis4_postv3AT, "/mnt/ingolialab/rmuller1/test_sandbox/gz_to_mpralm/rehearsal_dinner/scripts/GO_analysis_files/transcription_handcurhis4_postv3AT.txt")

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
his4_postv3AT_insig <- filter(his4_postv3AT, his4_postv3AT$adj.P.Val > 0.05)
temp1 <- pmax(-4, his4_postv3AT$logFC)
temp2 <- pmin(4, temp1)
his4_postv3AT$LFC_plot <- -temp2
his4_postv3AT$sig_plot <- -log10(pmax(10^-15, his4_postv3AT$adj.P.Val))

pdf("his4_postv3AT_satblock.pdf", useDingbats=FALSE, width=5, height=5, colormodel="rgb")
## Plot with type="n" to skip drawing the points
plot(x = his4_postv3AT$LFC_plot, y = his4_postv3AT$sig_plot,
     type="n",
     xlim=c(-4, 4), ylim=c(0, 15),
     axes=FALSE,
     xlab="Log2 Fold Before vs After Guide tet Induction", 
     ylab = expression('Statistical Significance -log10(Q)'), 
     col="gray75")
points(-pmin(4, aabiosynthesis_handcurhis4_postv3AT$logFC),
       -log10(pmax(10^-15, aabiosynthesis_handcurhis4_postv3AT$adj.P.Val)),
       type="n",
       col="mediumaquamarine"
)
points(-pmin(4, tRNAaasynthetase_handcurhis4_postv3AT$logFC),
       -log10(pmax(10^-15, tRNAaasynthetase_handcurhis4_postv3AT$adj.P.Val)),
       type="n",
       col="mediumaquamarine"
)
points(-pmin(4, translationcontrol_handcurhis4_postv3AT$logFC),
       -log10(pmax(10^-15, translationcontrol_handcurhis4_postv3AT$adj.P.Val)),
       type="n",
       col="mediumaquamarine"
)
points(-pmin(4, polIIIsubunit_handcurhis4_postv3AT$logFC),
       -log10(pmax(10^-15, polIIIsubunit_handcurhis4_postv3AT$adj.P.Val)),
       type="n",
       col="mediumaquamarine"
)
points(-pmin(4, tRNAprocessing_handcurhis4_postv3AT$logFC),
       -log10(pmax(10^-15, tRNAprocessing_handcurhis4_postv3AT$adj.P.Val)),
       type="n",
       col="mediumaquamarine"
)
points(-pmin(4, his4grna$logFC),
       -log10(pmax(10^-15, his4grna$adj.P.Val)),
       type="n",
       col="gray13"
)
points(-pmin(4, actin_handcurhis4_postv3AT$logFC),
       -log10(pmax(10^-15, actin_handcurhis4_postv3AT$adj.P.Val)),
       type="n",
       col="maroon"
)
#points(-pmin(4, ribosomal_subs_handcurhis4_postv3AT$logFC),
#       -log10(pmax(10^-15, ribosomal_subs_handcurhis4_postv3AT$adj.P.Val)),
#       type="n",
#       col="lightpink3"
#)

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
pointsfile <- "his4_postv3AT_satblock.png"
png(pointsfile, width=width, height=height, units="in", res=1500, bg="transparent")
par(mar=c(0,0,0,0))
## Exactly the same plot as above, but no type="n" and instead specify points, color, etc.
plot(x = his4_postv3AT$LFC_plot, y = his4_postv3AT$sig_plot,
     pch=16, cex=0.9,
     xlim=c(-4, 4), ylim=c(0, 15),
     axes=FALSE,
     xlab="Log2 Fold Before vs After Guide tet Induction",
     ylab = expression('Statistical Significance -log10(Q)'),
     col="gray75")
points(-pmin(4, aabiosynthesis_handcurhis4_postv3AT$logFC),
       -log10(pmax(10^-15, aabiosynthesis_handcurhis4_postv3AT$adj.P.Val)),
       pch=16, cex=0.9,
       col="mediumaquamarine"
)
points(-pmin(4, tRNAaasynthetase_handcurhis4_postv3AT$logFC),
       -log10(pmax(10^-15, tRNAaasynthetase_handcurhis4_postv3AT$adj.P.Val)),
       pch=16, cex=0.9,
       col="mediumaquamarine"
)
points(-pmin(4, translationcontrol_handcurhis4_postv3AT$logFC),
       -log10(pmax(10^-15, translationcontrol_handcurhis4_postv3AT$adj.P.Val)),
       pch=16, cex=0.9,
       col="mediumaquamarine"
)
points(-pmin(4, polIIIsubunit_handcurhis4_postv3AT$logFC),
       -log10(pmax(10^-15, polIIIsubunit_handcurhis4_postv3AT$adj.P.Val)),
       pch=16, cex=0.9,
       col="mediumaquamarine"
)
points(-pmin(4, tRNAprocessing_handcurhis4_postv3AT$logFC),
       -log10(pmax(10^-15, tRNAprocessing_handcurhis4_postv3AT$adj.P.Val)),
       pch=16, cex=0.9,
       col="mediumaquamarine"
)
points(-pmin(4, his4grna$logFC),
       -log10(pmax(10^-15, his4grna$adj.P.Val)),
       pch=16, cex=0.9,
       col="gray13"
)
points(-pmin(4, actin_handcurhis4_postv3AT$logFC),
       -log10(pmax(10^-15, actin_handcurhis4_postv3AT$adj.P.Val)),
       pch=16, cex=0.9,
       col="maroon"
)
#points(-pmin(4, ribosomal_subs_handcurhis4_postv3AT$logFC),
#       -log10(pmax(10^-15, ribosomal_subs_handcurhis4_postv3AT$adj.P.Val)),
#       pch=16, cex=0.9,
#       col="lightpink3"
#)
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
