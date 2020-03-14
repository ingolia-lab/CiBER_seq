#!/usr/bin/env Rscript

if (!requireNamespace("dplyr", quietly = TRUE))
  install.packages("dplyr")
library(dplyr)

mpralmpath <- "~/CiBER_seq_package/all_raw_fasta_gz/"
GOpath <- "~/CiBER_seq_package/scripts/GO_analysis_files/GO_annotation_lists/"
options(stringsAsFactors=FALSE)
for (annot in c("aabiosynthesis_handcur", "polIIIsubunit_handcur", "transcription_handcur", 
"translationcontrol_handcur", "tRNAaasynthetase_handcur", "tRNAprocessing_handcur", "glycolysis_handcur", "sumoylation_handcur", "nuc_pore_handcur")){
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

head(nuc_pore_handcur)
head(cds)
head(nuc_pore_handcurcds)

#write.table(transcription_handcurhis4_pooled, "/mnt/ingolialab/rmuller1/test_sandbox/gz_to_mpralm/rehearsal_dinner/scripts/GO_analysis_files/transcription_handcurhis4_pooled.txt")

## Hybrid PDF-PNG scatter plots
## based on https://jonathanchang.org/blog/how-to-partially-rasterize-a-figure-plotted-with-r/
## by jonathan.chang@monash.edu
## We use par(mar=...) and plot(...) rather than plot.window(mar=...) and points(...)

#install.packages("png")
library(png)
library(scales)
options(stringsAsFactors=FALSE)

#cds before vs after guide tet induction treatment_________________

## Set up values to plot
cds_insig <- filter(cds, cds$adj.P.Val > 0.05)
temp1 <- pmax(-4, cds$logFC)
temp2 <- pmin(4, temp1)
cds$LFC_plot <- -temp2
cds$sig_plot <- -log10(pmax(10^-15, cds$adj.P.Val))

pdf("cds_activators.pdf", useDingbats=FALSE, width=5, height=5, colormodel="rgb")
## Plot with type="n" to skip drawing the points
plot(x = cds$LFC_plot, y = cds$sig_plot,
     type="n",
     xlim=c(-4, 4), ylim=c(0, 15),
     axes=FALSE,
     xlab="Log2 Fold Before vs After Guide tet Induction", 
     ylab = expression('Statistical Significance -log10(Q)'), 
     col="gray52")
points(-pmin(4, nuc_pore_handcurcds$logFC),
       -log10(pmax(10^-15, nuc_pore_handcurcds$adj.P.Val)),
       type="n",
       col="royalblue1"
)
points(-pmin(4, sumoylation_handcurcds$logFC),
       -log10(pmax(10^-15, sumoylation_handcurcds$adj.P.Val)),
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
pointsfile <- "cds_activators.png"
png(pointsfile, width=width, height=height, units="in", res=1500, bg="transparent")
par(mar=c(0,0,0,0))
## Exactly the same plot as above, but no type="n" and instead specify points, color, etc.
plot(x = cds$LFC_plot, y = cds$sig_plot,
     pch=16, cex=0.9,
     xlim=c(-4, 4), ylim=c(0, 15),
     axes=FALSE,
     xlab="Log2 Fold Before vs After Guide tet Induction",
     ylab = expression('Statistical Significance -log10(Q)'),
     col="gray52")
points(-pmin(4, nuc_pore_handcurcds$logFC),
       -log10(pmax(10^-15, nuc_pore_handcurcds$adj.P.Val)),
       pch=16, cex=0.9,
       col="royalblue1"
)
points(-pmin(4, sumoylation_handcurcds$logFC),
       -log10(pmax(10^-15, sumoylation_handcurcds$adj.P.Val)),
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
title(main="TF-GCN4 CDS Guide Induction")
dev.off()
