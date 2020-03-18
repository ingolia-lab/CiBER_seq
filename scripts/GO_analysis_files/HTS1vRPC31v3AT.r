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
  for (dir_file in c("HTS1_RPC31/HTS1_sum_mpralm", "HTS1_RPC31/RPC31_sum_mpralm", "HIS4_PGK1_3AT/his4_postv3AT_sum_mpralm")){
    y <- read.delim(paste(mpralmpath, dir_file, ".txt", sep=""),
                    stringsAsFactors=FALSE, header = TRUE)
    temp1 <- gsub(".*/", "", dir_file)
    temp3 <- temp2 <- gsub("_sum_mpralm", "", temp1)
    assign(paste(temp2, sep=""), y)
    m <- subset(y, (y$Yorf1 %in% x$name))    
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

#HTS1 vs RPC31 logFC plot 

## Set up values to plot
HTS1_sub <- HTS1[c(1,2,6,9)] 
RPC31_sub <- RPC31[c(1,2,6,9)] 
thrAT_sub <- his4_postv3AT[c(1,2,6,9)]
names(HTS1_sub) <- c("Guide", "HTS1_logFC", "HTS1_adj.P.val", "HTS1_Yorf1")
names(RPC31_sub) <- c("Guide", "RPC31_logFC", "RPC31_adj.P.val", "RPC31_Yorf1")
names(thrAT_sub) <- c("Guide", "thrAT_logFC", "thrAT_adj.P.val", "thrAT_Yorf1")

head(HTS1_sub)
head(RPC31_sub)
head(thrAT_sub)

HTS1vRPC31 <- merge(HTS1_sub, RPC31_sub, by = "Guide")
head(HTS1vRPC31)

HTS1vRPC31v3AT <- merge(HTS1vRPC31, thrAT_sub, by = "Guide")
head(HTS1vRPC31v3AT)

HTS1_sig <- filter(HTS1vRPC31v3AT, HTS1vRPC31v3AT$HTS1_adj.P.val < 0.05)
RPC31_sig <- filter(HTS1vRPC31v3AT, HTS1vRPC31v3AT$RPC31_adj.P.val < 0.05)
both_sig <- filter(HTS1vRPC31v3AT, HTS1vRPC31v3AT$RPC31_adj.P.val < 0.05 & HTS1vRPC31v3AT$HTS1_adj.P.val < 0.05)
all_sig <- filter(HTS1vRPC31v3AT, HTS1vRPC31v3AT$RPC31_adj.P.val < 0.05 & HTS1vRPC31v3AT$HTS1_adj.P.val < 0.05 & HTS1vRPC31v3AT$thrAT_adj.P.val < 0.05)

pdf("HTS1vRPC31v3AT.pdf", useDingbats=FALSE, width=5, height=5, colormodel="rgb")
## Plot with type="n" to skip drawing the points
plot(x = -HTS1vRPC31v3AT$HTS1_logFC, y = -HTS1vRPC31v3AT$RPC31_logFC,
     type="n",
     axes=FALSE,
     xlab="HTS1 Log Fold Change", 
     ylab = "RPC31 Log Fold Change", 
     col="gray52")
points(-HTS1_sig$HTS1_logFC,
       -HTS1_sig$RPC31_logFC,
       type="n",
       col="firebrick"
)
points(-RPC31_sig$HTS1_logFC,
       -RPC31_sig$RPC31_logFC,
       type="n",
       col="royalblue"
)
points(-both_sig$HTS1_logFC,
       -both_sig$RPC31_logFC,
       type="n",
       col="purple3"
)
points(-all_sig$HTS1_logFC,
       -all_sig$RPC31_logFC,
       type="n",
       col="black"
)
abline(h = 0, col="black")
abline(v = 0, col="black")

## Capture the plot window coordinates
coords <- par("usr")
gx <- grconvertX(coords[1:2], "user", "inches")
gy <- grconvertY(coords[3:4], "user", "inches")
width <- max(gx) - min(gx)
height <- max(gy) - min(gy)

## Plot points to a PNG file
pointsfile <- "HTS1vRPC31.png"
png(pointsfile, width=width, height=height, units="in", res=1500, bg="transparent")
par(mar=c(0,0,0,0))
## Exactly the same plot as above, but no type="n" and instead specify points, color, etc.
plot(x = -HTS1vRPC31v3AT$HTS1_logFC, y = -HTS1vRPC31v3AT$RPC31_logFC,
     pch=16, cex=0.9,
     axes=FALSE,
     xlab="HTS1 Log Fold Change", 
     ylab = "RPC31 Log Fold Change", 
     col="gray52")
points(-HTS1_sig$HTS1_logFC,
       -HTS1_sig$RPC31_logFC,
       pch=16, cex=0.9,
       col="firebrick"
)
points(-RPC31_sig$HTS1_logFC,
       -RPC31_sig$RPC31_logFC,
       pch=16, cex=0.9,
       col="royalblue"
)
points(-both_sig$HTS1_logFC,
       -both_sig$RPC31_logFC,
       pch=16, cex=0.9,
       col="purple3"
)
points(-all_sig$HTS1_logFC,
       -all_sig$RPC31_logFC,
       pch=16, cex=0.9,
       col="black"
)
dev.off()
## Read in the PNG file and draw as a raster image
panel <- readPNG(pointsfile, native = TRUE)
rasterImage(panel, coords[1], coords[3], coords[2], coords[4])

## Draw axes, title, etc.
axis(side = 1, at = c(-6, -4, -2, 0, 2, 4, 6), labels=c("<1/64", "1/16", "1/4", "1", "4", "16", ">64"))
axis(side = 2, at = c(-6, -4, -2, 0, 2, 4, 6), labels=c("<1/64", "1/16", "1/4", "1", "4", "16", ">64"))
title(main="HTS1 vs RPC31 LogFC")
dev.off()