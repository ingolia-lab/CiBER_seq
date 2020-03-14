#!/usr/bin/env Rscript

if (!requireNamespace("dplyr", quietly = TRUE))
  install.packages("dplyr")
library(dplyr)

mpralmpath <- "~/CiBER_seq_package/all_raw_fasta_gz/"
GOpath <- "~/CiBER_seq_package/scripts/GO_analysis_files/GO_annotation_lists/"
options(stringsAsFactors=FALSE)
for (annot in c("aabiosynthesis_handcur", "polIIIsubunit_handcur", "transcription_handcur", 
                "translationcontrol_handcur", "tRNAaasynthetase_handcur", 
                "tRNAprocessing_handcur", "glycolysis_handcur")){
  x <- read.table(paste(GOpath, annot, ".txt", sep=""),
                  stringsAsFactors=FALSE, header = TRUE)
  assign(paste(annot, sep=""), x)
  for (dir_file in c("HIS4_PGK1_pooled/pgk1_pooled_sum_mpralm")){
    y <- read.delim(paste(mpralmpath, dir_file, ".txt", sep=""),
                    stringsAsFactors=FALSE, header = TRUE)
    temp1 <- gsub(".*/", "", dir_file)
    temp3 <- temp2 <- gsub("_sum_mpralm", "", temp1)
    assign(paste(temp2, sep=""), y)
    m <- subset(y, (y$Yorf1 %in% x$name))    
    assign(paste(annot, temp3, sep=""), m)
  }
}

pgk1grna <- filter(pgk1_pooled, pgk1_pooled$gene == "PGK1")

head(aabiosynthesis_handcur)
head(pgk1_pooled)
head(aabiosynthesis_handcurpgk1_pooled)

head(transcription_handcurpgk1_pooled)
#write.table(transcription_handcurpgk1_pooled, "/mnt/ingolialab/rmuller1/test_sandbox/gz_to_mpralm/rehearsal_dinner/scripts/GO_analysis_files/transcription_handcurpgk1_pooled.txt")

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
pgk1_pooled_insig <- filter(pgk1_pooled, pgk1_pooled$adj.P.Val > 0.05)
temp1 <- pmax(-4, pgk1_pooled$logFC)
temp2 <- pmin(4, temp1)
pgk1_pooled$LFC_plot <- -temp2
pgk1_pooled$sig_plot <- -log10(pmax(10^-15, pgk1_pooled$adj.P.Val))

pdf("pgk1_pool_activators.pdf", useDingbats=FALSE, width=5, height=5, colormodel="rgb")
## Plot with type="n" to skip drawing the points
plot(x = pgk1_pooled$LFC_plot, y = pgk1_pooled$sig_plot,
     type="n",
     xlim=c(-4, 4), ylim=c(0, 15),
     axes=FALSE,
     xlab="Log2 Fold Before vs After Guide tet Induction", 
     ylab = expression('Statistical Significance -log10(Q)'), 
     col="gray75")
points(-pmin(4, glycolysis_handcurpgk1_pooled$logFC),
       -log10(pmax(10^-15, glycolysis_handcurpgk1_pooled$adj.P.Val)),
       type="n",
       col="orange"
)
points(-pmin(4, transcription_handcurpgk1_pooled$logFC),
       -log10(pmax(10^-15, transcription_handcurpgk1_pooled$adj.P.Val)),
       type="n",
       col="red4"
)
points(-pmin(4, pgk1grna$logFC),
       -log10(pmax(10^-15, pgk1grna$adj.P.Val)),
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
pointsfile <- "pgk1_pool_activators.png"
png(pointsfile, width=width, height=height, units="in", res=1500, bg="transparent")
par(mar=c(0,0,0,0))
## Exactly the same plot as above, but no type="n" and instead specify points, color, etc.
plot(x = pgk1_pooled$LFC_plot, y = pgk1_pooled$sig_plot,
     pch=16, cex=0.9,
     xlim=c(-4, 4), ylim=c(0, 15),
     axes=FALSE,
     xlab="Log2 Fold Before vs After Guide tet Induction",
     ylab = expression('Statistical Significance -log10(Q)'),
     col="gray75")
points(-pmin(4, glycolysis_handcurpgk1_pooled$logFC),
       -log10(pmax(10^-15, glycolysis_handcurpgk1_pooled$adj.P.Val)),
       pch=16, cex=0.9,
       col="orange"
)
points(-pmin(4, transcription_handcurpgk1_pooled$logFC),
       -log10(pmax(10^-15, transcription_handcurpgk1_pooled$adj.P.Val)),
       pch=16, cex=0.9,
       col="red4"
)
points(-pmin(4, pgk1grna$logFC),
       -log10(pmax(10^-15, pgk1grna$adj.P.Val)),
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
axis(side = 1, at = c(-4, -2, 0, 2, 4), labels=c("<1/16", "1/4", "1", "4", ">16"))
axis(side = 2, at = c(0, 5, 10, 15), labels=c("0", "5", "10", ">15"))
title(main="P(PGK1) Guide Induction")
dev.off()
