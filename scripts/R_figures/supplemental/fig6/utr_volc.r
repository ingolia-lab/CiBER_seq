#!/usr/bin/env Rscript

if (!requireNamespace("dplyr", quietly = TRUE))
  install.packages("dplyr")
library(dplyr)

mpralmpath <- "~/CiBER_seq_package/all_raw_fastq/"
GOpath <- "~/CiBER_github/CiBER_seq/scripts/GO_analysis_files/GO_annotation_lists/"
options(stringsAsFactors=FALSE)
for (annot in c("aabiosynthesis_handcur", "polIIIsubunit_handcur", "transcription_handcur", 
                "translationcontrol_handcur", "tRNAaasynthetase_handcur", "tRNAprocessing_handcur", "glycolysis_handcur")){
  x <- read.table(paste(GOpath, annot, ".txt", sep=""),
                  stringsAsFactors=FALSE, header = TRUE)
  assign(paste(annot, sep=""), x)
  for (dir_file in c("GCN4_CDS_UTR/utr_sum_mpralm")){
    y <- read.delim(paste(mpralmpath, dir_file, ".txt", sep=""),
                    stringsAsFactors=FALSE, header = TRUE)
    temp1 <- gsub(".*/", "", dir_file)
    temp3 <- temp2 <- gsub("_sum_mpralm", "", temp1)
    assign(paste(temp2, sep=""), y)
    m <- subset(y, (y$Yorf1 %in% x$name))    
    assign(paste(annot, temp3, sep=""), m)
  }
}

gal1grna <- filter(utr, utr$gene == "GAL1")

head(aabiosynthesis_handcur)
head(utr)
head(aabiosynthesis_handcurutr)

head(transcription_handcurutr)
#write.table(transcription_handcurutr, "/mnt/ingolialab/rmuller1/test_sandbox/gz_to_mpralm/rehearsal_dinner/scripts/GO_analysis_files/transcription_handcurutr.txt")

## Hybrid PDF-PNG scatter plots
## based on https://jonathanchang.org/blog/how-to-partially-rasterize-a-figure-plotted-with-r/
## by jonathan.chang@monash.edu
## We use par(mar=...) and plot(...) rather than plot.window(mar=...) and points(...)

#install.packages("png")
library(png)
library(scales)
options(stringsAsFactors=FALSE)

#GCN4 UTR-TF before vs after guide tet induction treatment_________________

## Set up values to plot
utr_insig <- filter(utr, utr$adj.P.Val > 0.05)
temp1 <- pmax(-6, utr$logFC)
temp2 <- pmin(6, temp1)
utr$LFC_plot <- -temp2
utr$sig_plot <- -log10(pmax(10^-10, utr$adj.P.Val))

pdf("GCN4_cds_TF.pdf", useDingbats=FALSE, width=5, height=5, colormodel="rgb")
## Plot with type="n" to skip drawing the points
plot(x = utr$LFC_plot, y = utr$sig_plot,
     type="n",
     xlim=c(-6, 6), ylim=c(0, 10),
     axes=FALSE,
     xlab="Log2 Fold Before vs After Guide tet Induction", 
     ylab = expression('Statistical Significance -log10(Q)'), 
     col="gray75")
points(-pmin(6, aabiosynthesis_handcurutr$logFC),
       -log10(pmax(10^-10, aabiosynthesis_handcurutr$adj.P.Val)),
       type="n",
       col="royalblue1"
)
points(-pmin(6, tRNAaasynthetase_handcurutr$logFC),
       -log10(pmax(10^-10, tRNAaasynthetase_handcurutr$adj.P.Val)),
       type="n",
       col="seagreen4"
)
points(-pmin(6, translationcontrol_handcurutr$logFC),
       -log10(pmax(10^-10, translationcontrol_handcurutr$adj.P.Val)),
       type="n",
       col="paleturquoise"
)
points(-pmin(6, polIIIsubunit_handcurutr$logFC),
       -log10(pmax(10^-10, polIIIsubunit_handcurutr$adj.P.Val)),
       type="n",
       col="violet"
)
points(-pmin(6, tRNAprocessing_handcurutr$logFC),
       -log10(pmax(10^-10, tRNAprocessing_handcurutr$adj.P.Val)),
       type="n",
       col="red4"
)
points(-pmin(6, gal1grna$logFC),
       -log10(pmax(10^-10, gal1grna$adj.P.Val)),
       type="n",
       col="gray13"
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
pointsfile <- "cds_TF.png"
png(pointsfile, width=width, height=height, units="in", res=1500, bg="transparent")
par(mar=c(0,0,0,0))
## Exactly the same plot as above, but no type="n" and instead specify points, color, etc.
plot(x = utr$LFC_plot, y = utr$sig_plot,
     pch=16, cex=0.9,
     xlim=c(-6, 6), ylim=c(0, 10),
     axes=FALSE,
     xlab="Log2 Fold Before vs After Guide tet Induction",
     ylab = expression('Statistical Significance -log10(Q)'),
     col="gray75")
points(-pmin(6, aabiosynthesis_handcurutr$logFC),
       -log10(pmax(10^-10, aabiosynthesis_handcurutr$adj.P.Val)),
       pch=16, cex=0.9,
       col="royalblue1"
)
points(-pmin(6, tRNAaasynthetase_handcurutr$logFC),
       -log10(pmax(10^-10, tRNAaasynthetase_handcurutr$adj.P.Val)),
       pch=16, cex=0.9,
       col="seagreen4"
)
points(-pmin(6, translationcontrol_handcurutr$logFC),
       -log10(pmax(10^-10, translationcontrol_handcurutr$adj.P.Val)),
       pch=16, cex=0.9,
       col="paleturquoise"
)
points(-pmin(6, polIIIsubunit_handcurutr$logFC),
       -log10(pmax(10^-10, polIIIsubunit_handcurutr$adj.P.Val)),
       pch=16, cex=0.9,
       col="violet"
)
points(-pmin(6, tRNAprocessing_handcurutr$logFC),
       -log10(pmax(10^-10, tRNAprocessing_handcurutr$adj.P.Val)),
       pch=16, cex=0.9,
       col="red4"
)
points(-pmin(6, gal1grna$logFC),
       -log10(pmax(10^-10, gal1grna$adj.P.Val)),
       pch=16, cex=0.9,
       col="gray13"
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
axis(side = 1, at = c(-6, -4, -2, 0, 2, 4, 6), 
     labels=c("<1/64", "1/16", "1/4", "1", "4", "16", ">64"))
axis(side = 2, at = c(0, 5, 10), labels=c("0", "5", ">10"))
title(main="P(Synthetic) Guide Induction")
dev.off()
