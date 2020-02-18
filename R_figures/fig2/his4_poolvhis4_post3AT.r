#!/usr/bin/env Rscript

if (!requireNamespace("dplyr", quietly = TRUE))
  install.packages("dplyr")
library(dplyr)

mpralmpath <- "/mnt/ingolialab/rmuller1/test_sandbox/gz_to_mpralm/rehearsal_dinner/all_raw_fasta_gz/"
GOpath <- "/mnt/ingolialab/rmuller1/test_sandbox/gz_to_mpralm/rehearsal_dinner/scripts/GO_analysis_files/GO_annotation_lists/"
options(stringsAsFactors=FALSE)
for (annot in c("aabiosynthesis_handcur", "polIIIsubunit_handcur", "transcription_handcur", 
                "translationcontrol_handcur", "tRNAaasynthetase_handcur", "tRNAprocessing_handcur", "glycolysis_handcur")){
  x <- read.table(paste(GOpath, annot, ".txt", sep=""),
                  stringsAsFactors=FALSE, header = TRUE)
  assign(paste(annot, sep=""), x)
  for (dir_file in c("HIS4_PGK1_3AT/his4_postv3AT_sum_mpralm", "HIS4_PGK1_pooled/his4_pooled_sum_mpralm")){
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

#His_pooled vs his4_3AT logFC plot 

## Set up values to plot
hispool_sub <- his4_pooled[c(1,2,6,9)] 
his3AT_sub <- his4_postv3AT[c(1,2,6,9)] 
names(hispool_sub) <- c("Guide", "hispool_logFC", "hispool_adj.P.val", "hispool_Yorf1")
names(his3AT_sub) <- c("Guide", "his3AT_logFC", "his3AT_adj.P.val", "his3AT_Yorf1")

head(hispool_sub)
head(his3AT_sub)
     
hispoolv3AT <- merge(hispool_sub, his3AT_sub, by = "Guide")

head(hispoolv3AT)

hispool_sig <- filter(hispoolv3AT, hispoolv3AT$hispool_adj.P.val < 0.05)
his3AT_sig <- filter(hispoolv3AT, hispoolv3AT$his3AT_adj.P.val < 0.05)
both_sig <- filter(hispoolv3AT, hispoolv3AT$his3AT_adj.P.val < 0.05 & hispoolv3AT$hispool_adj.P.val < 0.05)

pdf("hispoolv3AT.pdf", useDingbats=FALSE, width=5, height=5, colormodel="rgb")
## Plot with type="n" to skip drawing the points
plot(x = -hispoolv3AT$hispool_logFC, y = -hispoolv3AT$his3AT_logFC,
     type="n",
     axes=FALSE,
     xlab="hispool Log Fold Change", 
     ylab = "his3AT Log Fold Change", 
     col="gray52")
points(-hispool_sig$hispool_logFC,
       -hispool_sig$his3AT_logFC,
       type="n",
       col="firebrick"
)
points(-his3AT_sig$hispool_logFC,
       -his3AT_sig$his3AT_logFC,
       type="n",
       col="royalblue"
)
points(-both_sig$hispool_logFC,
       -both_sig$his3AT_logFC,
       type="n",
       col="purple3"
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
pointsfile <- "his4_pooled_volcano.png"
png(pointsfile, width=width, height=height, units="in", res=1500, bg="transparent")
par(mar=c(0,0,0,0))
## Exactly the same plot as above, but no type="n" and instead specify points, color, etc.
plot(x = -hispoolv3AT$hispool_logFC, y = -hispoolv3AT$his3AT_logFC,
     pch=16, cex=0.9,
     axes=FALSE,
     xlab="hispool Log Fold Change", 
     ylab = "his3AT Log Fold Change", 
     col="gray52")
points(-hispool_sig$hispool_logFC,
       -hispool_sig$his3AT_logFC,
       pch=16, cex=0.9,
       col="firebrick"
)
points(-his3AT_sig$hispool_logFC,
       -his3AT_sig$his3AT_logFC,
       pch=16, cex=0.9,
       col="royalblue"
)
points(-both_sig$hispool_logFC,
       -both_sig$his3AT_logFC,
       pch=16, cex=0.9,
       col="purple3"
)
dev.off()
## Read in the PNG file and draw as a raster image
panel <- readPNG(pointsfile, native = TRUE)
rasterImage(panel, coords[1], coords[3], coords[2], coords[4])

## Draw axes, title, etc.
axis(side = 1, at = c(-6, -4, -2, 0, 2, 4, 6), labels=c("<1/64", "1/16", "1/4", "1", "4", "16", ">64"))
axis(side = 2, at = c(-6, -4, -2, 0, 2, 4, 6), labels=c("<1/64", "1/16", "1/4", "1", "4", "16", ">64"))
title(main="hispool vs his3AT LogFC")
dev.off()
