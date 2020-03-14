#!/usr/bin/env Rscript

if (!requireNamespace("dplyr", quietly = TRUE))
  install.packages("dplyr")
library(dplyr)

mpralmpath <- "~/CiBER_seq_package/all_raw_fasta_gz/"
GOpath <- "~/CiBER_seq_package/scripts/GO_analysis_files/GO_annotation_lists/"

his4_postv3AT <- read.delim(paste(mpralmpath, "HIS4_PGK1_3AT/his4_postv3AT_sum_mpralm", ".txt", sep=""),
                  stringsAsFactors=FALSE, header = TRUE)
his4_pooled <- read.delim(paste(mpralmpath, "HIS4_PGK1_pooled/his4_pooled_sum_mpralm", ".txt", sep=""),
                stringsAsFactors=FALSE, header = TRUE)
his4_postv3AT_sub <- his4_postv3AT[c(1,2,6,9)] 
his4_pooled_sub <- his4_pooled[c(1,2,6,9)] 
names(his4_postv3AT_sub) <- c("Guide", "his4_postv3AT_logFC", "his4_postv3AT_adj.P.val", "his4_postv3AT_Yorf1")
names(his4_pooled_sub) <- c("Guide", "his4_pooled_logFC", "his4_pooled_adj.P.val", "his4_pooled_Yorf1")

his4_poolv3AT <- merge(his4_postv3AT_sub, his4_pooled_sub, by = "Guide")
his4_poolv3AT_sig <- filter(his4_poolv3AT, 
                            his4_poolv3AT$his4_postv3AT_adj.P.val < 0.05 | his4_poolv3AT$his4_pooled_adj.P.val < 0.05)
his4_poolv3AT_insig <- filter(his4_poolv3AT, 
                            his4_poolv3AT$his4_postv3AT_adj.P.val > 0.05 & his4_poolv3AT$his4_pooled_adj.P.val > 0.05)

his4grna <- filter(his4_poolv3AT, his4_poolv3AT$his4_postv3AT_Yorf1 == "YCL030C")

head(his4_poolv3AT_sig)
head(his4_poolv3AT_insig)

options(stringsAsFactors=FALSE)
for (annot in c("aabiosynthesis_handcur", "polIIIsubunit_handcur", "transcription_handcur", 
                "translationcontrol_handcur", "tRNAaasynthetase_handcur", "tRNAprocessing_handcur", "glycolysis_handcur")){
  x <- read.table(paste(GOpath, annot, ".txt", sep=""),
                  stringsAsFactors=FALSE, header = TRUE)
  m <- subset(his4_poolv3AT_sig, (his4_poolv3AT_sig$his4_pooled_Yorf1 %in% x$name)) 
  assign(paste(annot, "his4_poolv3AT_sig", sep=""), m)
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

pdf("new_hispoolv3AT.pdf", useDingbats=FALSE, width=8.5, height=8.5, colormodel="rgb")
## Plot with type="n" to skip drawing the points
#plot(x = -his4_poolv3AT_insig$his4_pooled_logFC, y = -his4_poolv3AT_insig$his4_postv3AT_logFC,
#     type="n",
#     axes=FALSE,
#     xlab="his4 pooled Log Fold Change",
#     ylab = "his3AT Log Fold Change",
#     col="gray75"
#)
plot(x = -his4_poolv3AT_sig$his4_pooled_logFC, y = -his4_poolv3AT_sig$his4_postv3AT_logFC,
     type="n",
     xlim = c(-4.5,4.5),
     ylim = c(-4.5,4.5), 
     axes=FALSE,
     xlab="his4 pooled Log Fold Change", 
     ylab = "his3AT Log Fold Change", 
     col="gray75"
     )
points(x = -aabiosynthesis_handcurhis4_poolv3AT_sig$his4_pooled_logFC, 
       y = -aabiosynthesis_handcurhis4_poolv3AT_sig$his4_postv3AT_logFC,
       type="n",
       col="royalblue1"
)
points(x = -polIIIsubunit_handcurhis4_poolv3AT_sig$his4_pooled_logFC, 
       y = -polIIIsubunit_handcurhis4_poolv3AT_sig$his4_postv3AT_logFC,
       type="n",
       col="violet"
)
points(x = -tRNAprocessing_handcurhis4_poolv3AT_sig$his4_pooled_logFC, 
       y = -tRNAprocessing_handcurhis4_poolv3AT_sig$his4_postv3AT_logFC,
       type="n",
       col="red4"
)
points(x = -tRNAaasynthetase_handcurhis4_poolv3AT_sig$his4_pooled_logFC, 
       y = -tRNAaasynthetase_handcurhis4_poolv3AT_sig$his4_postv3AT_logFC,
       type="n",
       col="seagreen4"
)
points(x = -translationcontrol_handcurhis4_poolv3AT_sig$his4_pooled_logFC, 
       y = -translationcontrol_handcurhis4_poolv3AT_sig$his4_postv3AT_logFC,
       type="n",
       col="paleturquoise"
)
points(x = -his4grna$his4_pooled_logFC,
       y = -his4grna$his4_postv3AT_logFC,
       type="n",
       col="gray13"
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
pointsfile <- "his4_pooled_vs3AT.png"
png(pointsfile, width=width, height=height, units="in", res=1500, bg="transparent")
par(mar=c(0,0,0,0))
## Exactly the same plot as above, but no type="n" and instead specify points, color, etc.
#plot(x = -his4_poolv3AT_insig$his4_pooled_logFC, y = -his4_poolv3AT_insig$his4_postv3AT_logFC,
#     pch=16, cex=0.9,
#     axes=FALSE,
#     xlab="his4 pooled Log Fold Change",
#     ylab = "his3AT Log Fold Change",
#     col="gray75"
#)
plot(x = -his4_poolv3AT_sig$his4_pooled_logFC, y = -his4_poolv3AT_sig$his4_postv3AT_logFC,
     pch=16, cex=0.9,
     xlab="his4 pooled Log Fold Change",
     ylab = "his3AT Log Fold Change",
     xlim = c(-4.5,4.5),
     ylim = c(-4.5,4.5), 
     axes=FALSE,
     col="gray75"
     )
points(x = -aabiosynthesis_handcurhis4_poolv3AT_sig$his4_pooled_logFC,
       y = -aabiosynthesis_handcurhis4_poolv3AT_sig$his4_postv3AT_logFC,
       pch=16, cex=0.9,
       col="royalblue1"
)
points(x = -polIIIsubunit_handcurhis4_poolv3AT_sig$his4_pooled_logFC,
       y = -polIIIsubunit_handcurhis4_poolv3AT_sig$his4_postv3AT_logFC,
       pch=16, cex=0.9,
       col="violet"
)
points(x = -tRNAprocessing_handcurhis4_poolv3AT_sig$his4_pooled_logFC,
       y = -tRNAprocessing_handcurhis4_poolv3AT_sig$his4_postv3AT_logFC,
       pch=16, cex=0.9,
       col="red4"
)
points(x = -tRNAaasynthetase_handcurhis4_poolv3AT_sig$his4_pooled_logFC,
       y = -tRNAaasynthetase_handcurhis4_poolv3AT_sig$his4_postv3AT_logFC,
       pch=16, cex=0.9,
       col="seagreen4"
)
points(x = -translationcontrol_handcurhis4_poolv3AT_sig$his4_pooled_logFC,
       y = -translationcontrol_handcurhis4_poolv3AT_sig$his4_postv3AT_logFC,
       pch=16, cex=0.9,
       col="paleturquoise"
)
points(x = -his4grna$his4_pooled_logFC,
       y = -his4grna$his4_postv3AT_logFC,
       pch=16, cex=0.9,
       col="gray13"
)
dev.off()
## Read in the PNG file and draw as a raster image
panel <- readPNG(pointsfile, native = TRUE)
rasterImage(panel, coords[1], coords[3], coords[2], coords[4])

## Draw axes, title, etc.
axis(side = 1, at = c(-4, -3, -2, -1, 0, 1, 2, 3, 4), labels=c("<1/16", "1/8", "1/4", "1/2", "1", "2", "4", "8", ">16"))
axis(side = 2, at = c(-4, -3, -2, -1, 0, 1, 2, 3, 4), labels=c("<1/16", "1/8", "1/4", "1/2", "1", "2", "4", "8", ">16"))
title(main="hispool vs his3AT LogFC")
dev.off()
