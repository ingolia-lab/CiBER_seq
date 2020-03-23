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
sgd_sub <- sgd[c(4,5,16)]
names(sgd_sub) <- c("name", "gene", "desc")

mpralmpath <- "~/CiBER_seq_package/all_raw_fastq/"
GOpath <- "~/CiBER_github/CiBER_seq/scripts/GO_analysis_files/GO_annotation_lists/"
options(stringsAsFactors=FALSE)
for (annot in c("aabiosynthesis_handcur", "polIIIsubunit_handcur", "transcription_handcur", 
                "translationcontrol_handcur", "tRNAaasynthetase_handcur", "tRNAprocessing_handcur",
                "TRiC_CCT_handcur", "proteasome_handcur", "actin_handcur", 
                "ER_traff_handcur")){
  x <- read.table(paste(GOpath, annot, ".txt", sep=""),
                  stringsAsFactors=FALSE, header = TRUE)
  assign(paste(annot, sep=""), x)
  for (dir_file in c("epistasis_3AT_HTS1_RPC31/his4_postv3AT_sum_mpralm", "epistasis_3AT_HTS1_RPC31/his4postvHTS1post_sum_mpralm", 
                     "epistasis_3AT_HTS1_RPC31/his4postvRPC31post_sum_mpralm", "HIS4_PGK1_pooled/his4_pooled_sum_mpralm")){
    y <- read.delim(paste(mpralmpath, dir_file, ".txt", sep=""),
                    stringsAsFactors=FALSE, header = TRUE)
    temp1 <- gsub(".*/", "", dir_file)
    temp3 <- temp2 <- gsub("_sum_mpralm", "", temp1)
    assign(paste(temp2, sep=""), y)
    m <- subset(y, (y$Yorf1 %in% x$name))    
    assign(paste(annot, temp3, sep=""), m)
  }
}

#merge 3AT, HTS1, RPC31
postvthrAT_sub <- his4_postv3AT[c(1,2,6,9,23,24)]

postvHTS1_sub <- data.frame(his4postvHTS1post$Guide, his4postvHTS1post$logFC, 
                            his4postvHTS1post$adj.P.Val, his4postvHTS1post$Yorf1)
names(postvHTS1_sub) <- c("Guide", "HTS1_logFC", "HTS1_adj.P.Val", "HTS1_Yorf1")

postvRPC31_sub <- data.frame(his4postvRPC31post$Guide, his4postvRPC31post$logFC, 
                             his4postvRPC31post$adj.P.Val, his4postvRPC31post$Yorf1)
names(postvRPC31_sub) <- c("Guide", "RPC31_logFC", "RPC31_adj.P.Val", "RPC31_Yorf1")

all <- merge(postvthrAT_sub, postvHTS1_sub, by="Guide")
all <- merge(all, postvRPC31_sub, by="Guide")
a <- cor(-all$logFC, all$HTS1_logFC)
all <- filter(all, all$adj.P.Val < 0.05 | all$HTS1_adj.P.Val < 0.05)

ER_traff <- subset(all, (all$HTS1_Yorf1 %in% ER_traff_handcur$name))
TRiC <- subset(all, (all$HTS1_Yorf1 %in% TRiC_CCT_handcur$name))
prot <- subset(all, (all$HTS1_Yorf1 %in% proteasome_handcur$name))
actin <- subset(all, (all$HTS1_Yorf1 %in% actin_handcur$name))
aabiosyn <- subset(all, (all$HTS1_Yorf1 %in% aabiosynthesis_handcur$name))
polIII <- subset(all, (all$HTS1_Yorf1 %in% polIIIsubunit_handcur$name))
transcr <- subset(all, (all$HTS1_Yorf1 %in% transcription_handcur$name))
transl <- subset(all, (all$HTS1_Yorf1 %in% translationcontrol_handcur$name))
tRNAsynth <- subset(all, (all$HTS1_Yorf1 %in% tRNAaasynthetase_handcur$name))
tRNAproc <- subset(all, (all$HTS1_Yorf1 %in% tRNAprocessing_handcur$name))

pcl5grna <- filter(all, all$HTS1_Yorf1 == "YHR071W")
yih1grna <- filter(all, all$HTS1_Yorf1 == "YCR059C")
gcn4grna <- filter(all, all$HTS1_Yorf1 == "YEL009C")

## Hybrid PDF-PNG scatter plots
## based on https://jonathanchang.org/blog/how-to-partially-rasterize-a-figure-plotted-with-r/
## by jonathan.chang@monash.edu
## We use par(mar=...) and plot(...) rather than plot.window(mar=...) and points(...)

#install.packages("png")
library(png)
library(scales)
options(stringsAsFactors=FALSE)

pdf("epistasis_3AT_HTS1.pdf", useDingbats=FALSE, width=5, height=5, colormodel="rgb")
## Plot with type="n" to skip drawing the points
plot(x = -all$logFC, y = all$HTS1_logFC,
     type="n",
     axes=FALSE,
     xlim = c(-6,6),
     ylim = c(-6,6), 
     xlab="postv3AT Fold Change", 
     ylab = "postvHTS1 Fold Change", 
     col="gray75")
points(-aabiosyn$logFC,
       aabiosyn$HTS1_logFC,
       type="n",
       col="mediumaquamarine"
)
points(-transl$logFC,
       transl$HTS1_logFC,
       type="n",
       col="mediumaquamarine"
)
points(-tRNAsynth$logFC,
       tRNAsynth$HTS1_logFC,
       type="n",
       col="mediumaquamarine"
)
points(-polIII$logFC,
       polIII$HTS1_logFC,
       type="n",
       col="mediumaquamarine"
)
points(-tRNAproc$logFC,
       tRNAproc$HTS1_logFC,
       type="n",
       col="mediumaquamarine"
)
#points(-TRiC$logFC,
#       TRiC$HTS1_logFC,
#       type="n",
#       col="tan4"
#)
#points(-prot$logFC,
#       prot$HTS1_logFC,
#       type="n",
#       col="orange1"
#)
points(-actin$logFC,
       actin$HTS1_logFC,
       type="n",
       col="maroon"
)
#points(-ER_traff$logFC,
#       ER_traff$HTS1_logFC,
#       type="n",
#       col="maroon"
#)
points(-pcl5grna$logFC,
       pcl5grna$HTS1_logFC,
       type="n",
       col="gray5"
)
points(-yih1grna$logFC,
       yih1grna$HTS1_logFC,
       type="n",
       col="gray5"
)
points(-gcn4grna$logFC,
       gcn4grna$HTS1_logFC,
       type="n",
       col="gray5"
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
pointsfile <- "epistasis3ATvsHTS1.png"
png(pointsfile, width=width, height=height, units="in", res=1500, bg="transparent")
par(mar=c(0,0,0,0))
## Exactly the same plot as above, but no type="n" and instead specify points, color, etc.
plot(x = -all$logFC, y = all$HTS1_logFC,
     pch=16, cex=0.9,
     axes=FALSE,
     xlim = c(-6,6),
     ylim = c(-6,6), 
     xlab="3AT Log Fold Change", 
     ylab = "RPC31 Log Fold Change", 
     col="gray75")
points(-aabiosyn$logFC,
       aabiosyn$HTS1_logFC,
       pch=16, cex=0.9,
       col="mediumaquamarine"
)
points(-transl$logFC,
       transl$HTS1_logFC,
       pch=16, cex=0.9,
       col="mediumaquamarine"
)
points(-tRNAsynth$logFC,
       tRNAsynth$HTS1_logFC,
       pch=16, cex=0.9,
       col="mediumaquamarine"
)
points(-polIII$logFC,
       polIII$HTS1_logFC,
       pch=16, cex=0.9,
       col="mediumaquamarine"
)
points(-tRNAproc$logFC,
       tRNAproc$HTS1_logFC,
       pch=16, cex=0.9,
       col="mediumaquamarine"
)
#points(-TRiC$logFC,
#       TRiC$HTS1_logFC,
#       pch=16, cex=0.9,
#       col="tan4"
#)
#points(-prot$logFC,
#       prot$HTS1_logFC,
#       pch=16, cex=0.9,
#       col="orange1"
#)
points(-actin$logFC,
       actin$HTS1_logFC,
       pch=16, cex=0.9,
       col="maroon"
)
#points(-ER_traff$logFC,
#       ER_traff$HTS1_logFC,
#       pch=16, cex=0.9,
#       col="maroon"
#)
points(-pcl5grna$logFC,
       pcl5grna$HTS1_logFC,
       pch=16, cex=0.9,
       col="gray5"
)
points(-yih1grna$logFC,
       yih1grna$HTS1_logFC,
       pch=16, cex=0.9,
       col="gray5"
)
points(-gcn4grna$logFC,
       gcn4grna$HTS1_logFC,
       pch=16, cex=0.9,
       col="gray5"
)
dev.off()
## Read in the PNG file and draw as a raster image
panel <- readPNG(pointsfile, native = TRUE)
rasterImage(panel, coords[1], coords[3], coords[2], coords[4])

## Draw axes, title, etc.
axis(side = 1, at = c(-6, -4, -2, 0, 2, 4, 6), labels=c("<1/64", "1/16", "1/4", "1", "4", "16", ">64"))
axis(side = 2, at = c(-6, -4, -2, 0, 2, 4, 6), labels=c("<1/64", "1/16", "1/4", "1", "4", "16", ">64"))
title(main="postv3AT vs postvHTS1 LogFC")
dev.off()

#__________________________
# postv3AT vs postvRPC31

all <- merge(postvthrAT_sub, postvHTS1_sub, by="Guide")
all <- merge(all, postvRPC31_sub, by="Guide")
b <- cor(-all$logFC, all$RPC31_logFC)
all <- filter(all, all$adj.P.Val < 0.05 | all$RPC31_adj.P.Val < 0.05)

ER_traff <- subset(all, (all$HTS1_Yorf1 %in% ER_traff_handcur$name))
TRiC <- subset(all, (all$HTS1_Yorf1 %in% TRiC_CCT_handcur$name))
prot <- subset(all, (all$HTS1_Yorf1 %in% proteasome_handcur$name))
actin <- subset(all, (all$HTS1_Yorf1 %in% actin_handcur$name))
aabiosyn <- subset(all, (all$HTS1_Yorf1 %in% aabiosynthesis_handcur$name))
polIII <- subset(all, (all$HTS1_Yorf1 %in% polIIIsubunit_handcur$name))
transcr <- subset(all, (all$HTS1_Yorf1 %in% transcription_handcur$name))
transl <- subset(all, (all$HTS1_Yorf1 %in% translationcontrol_handcur$name))
tRNAsynth <- subset(all, (all$HTS1_Yorf1 %in% tRNAaasynthetase_handcur$name))
tRNAproc <- subset(all, (all$HTS1_Yorf1 %in% tRNAprocessing_handcur$name))

pcl5grna <- filter(all, all$HTS1_Yorf1 == "YHR071W")
yih1grna <- filter(all, all$HTS1_Yorf1 == "YCR059C")
gcn4grna <- filter(all, all$HTS1_Yorf1 == "YEL009C")

pdf("epistasis_3AT_RPC31.pdf", useDingbats=FALSE, width=5, height=5, colormodel="rgb")
## Plot with type="n" to skip drawing the points
plot(x = -all$logFC, y = all$RPC31_logFC,
     type="n",
     axes=FALSE,
     xlim = c(-6,6),
     ylim = c(-6,6), 
     xlab="postv3AT Fold Change", 
     ylab = "postvRPC31 Fold Change", 
     col="gray75")
points(-aabiosyn$logFC,
       aabiosyn$RPC31_logFC,
       type="n",
       col="mediumaquamarine"
)
points(-transl$logFC,
       transl$RPC31_logFC,
       type="n",
       col="mediumaquamarine"
)
points(-tRNAsynth$logFC,
       tRNAsynth$RPC31_logFC,
       type="n",
       col="mediumaquamarine"
)
points(-polIII$logFC,
       polIII$RPC31_logFC,
       type="n",
       col="mediumaquamarine"
)
points(-tRNAproc$logFC,
       tRNAproc$RPC31_logFC,
       type="n",
       col="mediumaquamarine"
)
#points(-TRiC$logFC,
#       TRiC$RPC31_logFC,
#       type="n",
#       col="tan4"
#)
#points(-prot$logFC,
#       prot$RPC31_logFC,
#       type="n",
#       col="orange1"
#)
points(-actin$logFC,
       actin$RPC31_logFC,
       type="n",
       col="maroon"
)
#points(-ER_traff$logFC,
#       ER_traff$RPC31_logFC,
#       type="n",
#       col="maroon"
#)
points(-pcl5grna$logFC,
       pcl5grna$RPC31_logFC,
       type="n",
       col="gray5"
)
points(-yih1grna$logFC,
       yih1grna$RPC31_logFC,
       type="n",
       col="gray5"
)
points(-gcn4grna$logFC,
       gcn4grna$RPC31_logFC,
       type="n",
       col="gray5"
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
pointsfile <- "epistasis3ATvsRPC31.png"
png(pointsfile, width=width, height=height, units="in", res=1500, bg="transparent")
par(mar=c(0,0,0,0))
## Exactly the same plot as above, but no type="n" and instead specify points, color, etc.
plot(x = -all$logFC, y = all$RPC31_logFC,
     pch=16, cex=0.9,
     axes=FALSE,
     xlim = c(-6,6),
     ylim = c(-6,6), 
     xlab="RPC31 Log Fold Change", 
     ylab = "RPC31 Log Fold Change", 
     col="gray75")
points(-aabiosyn$logFC,
       aabiosyn$RPC31_logFC,
       pch=16, cex=0.9,
       col="mediumaquamarine"
)
points(-transl$logFC,
       transl$RPC31_logFC,
       pch=16, cex=0.9,
       col="mediumaquamarine"
)
points(-tRNAsynth$logFC,
       tRNAsynth$RPC31_logFC,
       pch=16, cex=0.9,
       col="mediumaquamarine"
)
points(-polIII$logFC,
       polIII$RPC31_logFC,
       pch=16, cex=0.9,
       col="mediumaquamarine"
)
points(-tRNAproc$logFC,
       tRNAproc$RPC31_logFC,
       pch=16, cex=0.9,
       col="mediumaquamarine"
)
#points(-TRiC$logFC,
#       TRiC$RPC31_logFC,
#       pch=16, cex=0.9,
#       col="tan4"
#)
#points(-prot$logFC,
#       prot$RPC31_logFC,
#       pch=16, cex=0.9,
#       col="orange1"
#)
points(-actin$logFC,
       actin$RPC31_logFC,
       pch=16, cex=0.9,
       col="maroon"
)
#points(-ER_traff$logFC,
#       ER_traff$RPC31_logFC,
#       pch=16, cex=0.9,
#       col="maroon"
#)
points(-pcl5grna$logFC,
       pcl5grna$RPC31_logFC,
       pch=16, cex=0.9,
       col="gray5"
)
points(-yih1grna$logFC,
       yih1grna$RPC31_logFC,
       pch=16, cex=0.9,
       col="gray5"
)
points(-gcn4grna$logFC,
       gcn4grna$RPC31_logFC,
       pch=16, cex=0.9,
       col="gray5"
)
dev.off()
## Read in the PNG file and draw as a raster image
panel <- readPNG(pointsfile, native = TRUE)
rasterImage(panel, coords[1], coords[3], coords[2], coords[4])

## Draw axes, title, etc.
axis(side = 1, at = c(-6, -4, -2, 0, 2, 4, 6), labels=c("<1/64", "1/16", "1/4", "1", "4", "16", ">64"))
axis(side = 2, at = c(-6, -4, -2, 0, 2, 4, 6), labels=c("<1/64", "1/16", "1/4", "1", "4", "16", ">64"))
title(main="postv3AT vs postvRPC31 LogFC")
dev.off()

#______________________ postvHTS1 vs postvRPC31

all <- merge(postvthrAT_sub, postvHTS1_sub, by="Guide")
all <- merge(all, postvRPC31_sub, by="Guide")
c <- cor(all$HTS1_logFC, all$RPC31_logFC)
all <- filter(all, all$HTS1_adj.P.Val < 0.05 | all$RPC31_adj.P.Val < 0.05)

ER_traff <- subset(all, (all$HTS1_Yorf1 %in% ER_traff_handcur$name))
TRiC <- subset(all, (all$HTS1_Yorf1 %in% TRiC_CCT_handcur$name))
prot <- subset(all, (all$HTS1_Yorf1 %in% proteasome_handcur$name))
actin <- subset(all, (all$HTS1_Yorf1 %in% actin_handcur$name))
aabiosyn <- subset(all, (all$HTS1_Yorf1 %in% aabiosynthesis_handcur$name))
polIII <- subset(all, (all$HTS1_Yorf1 %in% polIIIsubunit_handcur$name))
transcr <- subset(all, (all$HTS1_Yorf1 %in% transcription_handcur$name))
transl <- subset(all, (all$HTS1_Yorf1 %in% translationcontrol_handcur$name))
tRNAsynth <- subset(all, (all$HTS1_Yorf1 %in% tRNAaasynthetase_handcur$name))
tRNAproc <- subset(all, (all$HTS1_Yorf1 %in% tRNAprocessing_handcur$name))

pcl5grna <- filter(all, all$HTS1_Yorf1 == "YHR071W")
yih1grna <- filter(all, all$HTS1_Yorf1 == "YCR059C")
gcn4grna <- filter(all, all$HTS1_Yorf1 == "YEL009C")

pdf("epistasis_HTS1_RPC31.pdf", useDingbats=FALSE, width=5, height=5, colormodel="rgb")
## Plot with type="n" to skip drawing the points
plot(x = all$HTS1_logFC, y = all$RPC31_logFC,
     type="n",
     axes=FALSE,
     xlim = c(-6,6),
     ylim = c(-6,6), 
     xlab="postvHTS1 Fold Change", 
     ylab = "postvRPC31 Fold Change", 
     col="gray75")
points(aabiosyn$HTS1_logFC,
       aabiosyn$RPC31_logFC,
       type="n",
       col="mediumaquamarine"
)
points(transl$HTS1_logFC,
       transl$RPC31_logFC,
       type="n",
       col="mediumaquamarine"
)
points(tRNAsynth$HTS1_logFC,
       tRNAsynth$RPC31_logFC,
       type="n",
       col="mediumaquamarine"
)
points(polIII$HTS1_logFC,
       polIII$RPC31_logFC,
       type="n",
       col="mediumaquamarine"
)
points(tRNAproc$HTS1_logFC,
       tRNAproc$RPC31_logFC,
       type="n",
       col="mediumaquamarine"
)
#points(TRiC$HTS1_logFC,
#       TRiC$RPC31_logFC,
#       type="n",
#       col="tan4"
#)
#points(prot$HTS1_logFC,
#       prot$RPC31_logFC,
#       type="n",
#       col="orange1"
#)
points(actin$HTS1_logFC,
       actin$RPC31_logFC,
       type="n",
       col="maroon"
)
#points(ER_traff$HTS1_logFC,
#       ER_traff$RPC31_logFC,
#       type="n",
#       col="maroon"
#)
points(pcl5grna$HTS1_logFC,
       pcl5grna$RPC31_logFC,
       type="n",
       col="gray5"
)
points(yih1grna$HTS1_logFC,
       yih1grna$RPC31_logFC,
       type="n",
       col="gray5"
)
points(gcn4grna$HTS1_logFC,
       gcn4grna$RPC31_logFC,
       type="n",
       col="gray5"
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
pointsfile <- "epistasisHTS1vsRPC31.png"
png(pointsfile, width=width, height=height, units="in", res=1500, bg="transparent")
par(mar=c(0,0,0,0))
## Exactly the same plot as above, but no type="n" and instead specify points, color, etc.
plot(x = all$HTS1_logFC, y = all$RPC31_logFC,
     pch=16, cex=0.9,
     axes=FALSE,
     xlim = c(-6,6),
     ylim = c(-6,6), 
     xlab="HTS1 Log Fold Change", 
     ylab = "RPC31 Log Fold Change", 
     col="gray75")
points(aabiosyn$HTS1_logFC,
       aabiosyn$RPC31_logFC,
       pch=16, cex=0.9,
       col="mediumaquamarine"
)
points(transl$HTS1_logFC,
       transl$RPC31_logFC,
       pch=16, cex=0.9,
       col="mediumaquamarine"
)
points(tRNAsynth$HTS1_logFC,
       tRNAsynth$RPC31_logFC,
       pch=16, cex=0.9,
       col="mediumaquamarine"
)
points(polIII$HTS1_logFC,
       polIII$RPC31_logFC,
       pch=16, cex=0.9,
       col="mediumaquamarine"
)
points(tRNAproc$HTS1_logFC,
       tRNAproc$RPC31_logFC,
       pch=16, cex=0.9,
       col="mediumaquamarine"
)
#points(TRiC$HTS1_logFC,
#       TRiC$RPC31_logFC,
#       pch=16, cex=0.9,
#       col="tan4"
#)
#points(prot$HTS1_logFC,
#       prot$RPC31_logFC,
#       pch=16, cex=0.9,
#       col="orange1"
#)
points(actin$HTS1_logFC,
       actin$RPC31_logFC,
       pch=16, cex=0.9,
       col="maroon"
)
#points(ER_traff$HTS1_logFC,
#       ER_traff$RPC31_logFC,
#       pch=16, cex=0.9,
#       col="maroon"
#)
points(pcl5grna$HTS1_logFC,
       pcl5grna$RPC31_logFC,
       pch=16, cex=0.9,
       col="gray5"
)
points(yih1grna$HTS1_logFC,
       yih1grna$RPC31_logFC,
       pch=16, cex=0.9,
       col="gray5"
)
points(gcn4grna$HTS1_logFC,
       gcn4grna$RPC31_logFC,
       pch=16, cex=0.9,
       col="gray5"
)
dev.off()
## Read in the PNG file and draw as a raster image
panel <- readPNG(pointsfile, native = TRUE)
rasterImage(panel, coords[1], coords[3], coords[2], coords[4])

## Draw axes, title, etc.
axis(side = 1, at = c(-6, -4, -2, 0, 2, 4, 6), labels=c("<1/64", "1/16", "1/4", "1", "4", "16", ">64"))
axis(side = 2, at = c(-6, -4, -2, 0, 2, 4, 6), labels=c("<1/64", "1/16", "1/4", "1", "4", "16", ">64"))
title(main="postvHTS1 vs postvRPC31 LogFC")
dev.off()

a
b
c

#temp <- merge(postvthrAT_sub, postvHTS1_sub, by="Guide")
#temp <- filter(temp, temp$adj.P.Val < 0.05 | temp$HTS1_adj.P.Val < 0.05)
#plot(-temp$logFC, temp$HTS1_logFC, pch=20, col="blue")

#temp2 <- merge(postvthrAT_sub, postvRPC31_sub, by="Guide")
#temp2 <- filter(temp2, temp2$adj.P.Val < 0.05 | temp2$RPC31_adj.P.Val < 0.05)
#plot(-temp2$logFC, temp2$RPC31_logFC, pch=20, col="blue")
#temp2$diff <- temp2$logFC + temp2$RPC31_logFC

#_________________________________
#filter HTS1 and RCP31 that block activation post_vs_HTS1post and post_vs_RPC31post
HTS1_block <- all
RPC31_block <- all

#filter out simple activators seen in (LFC in HIS4 pre v post > 0.5)
'%notin%' <- Negate('%in%')
his4_pooled_act <- filter(his4_pooled, his4_pooled$adj.P.Val < 0.05 & his4_pooled$logFC < -0.5)
HTS1_act_no_ind <- subset(all, (all$HTS1_Yorf1 %in% his4_pooled_act$Yorf1))
HTS1_block <- subset(all, (all$HTS1_Yorf1 %notin% his4_pooled_act$Yorf1))

#GO analysis on those that block activation
GO_HTS1_block <- HTS1_block
GO_HTS1_block <- filter(GO_HTS1_block, 
                        GO_HTS1_block$HTS1_adj.P.Val < 0.05 & GO_HTS1_block$HTS1_logFC < -0.5)
GO_HTS1_block <- GO_HTS1_block[!is.na(GO_HTS1_block[,4]),]
write.table(as.factor(GO_HTS1_block$Yorf1), "~/GO_HTS1_block.txt", sep="\t", 
            row.names = FALSE, col.names = FALSE)


