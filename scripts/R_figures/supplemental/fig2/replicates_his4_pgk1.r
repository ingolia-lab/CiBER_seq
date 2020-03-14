#graph replicates RNA and DNA

mpralmpath <- "~/CiBER_seq_package/all_raw_fasta_gz/"

his4_pool_reads <- read.delim(paste(mpralmpath, 
                                "HIS4_PGK1_pooled/more_seq/all_his4_moreseq_counts", 
                                ".txt", sep=""),
                          stringsAsFactors=FALSE, header = TRUE)

pgk1_pool_reads <- read.delim(paste(mpralmpath, 
                                "HIS4_PGK1_pooled/more_seq/all_pgk1_moreseq_counts", 
                                ".txt", sep=""),
                          stringsAsFactors=FALSE, header = TRUE)

if (!requireNamespace("dplyr", quietly = TRUE))
  install.packages("dplyr")
library(dplyr)

if (!requireNamespace("scales", quietly = TRUE))
  install.packages("scales")
library(scales)

if (!requireNamespace("png", quietly = TRUE))
  install.packages("png")
library(png)

his4_pool_reads <- filter(his4_pool_reads, 
                          his4_pool_reads$IVT_PreL_S50_L008_R1_001his4.count.txt > 32 & his4_pool_reads$IVT_PreR_S51_L008_R1_001his4.count.txt > 32)
pgk1_pool_reads <- filter(pgk1_pool_reads, 
                          pgk1_pool_reads$IVT_PreL_S50_L008_R1_001pgk1.count.txt > 32 & pgk1_pool_reads$IVT_PreR_S51_L008_R1_001pgk1.count.txt > 32)

#his4 IVT pre replicates__________________________________________________
pdf("his4_IVT_pre_replicates.pdf", useDingbats=FALSE, width=5, height=5, colormodel="rgb")
par(pty="s")
plot(log10(his4_pool_reads$IVT_PreL_S50_L008_R1_001his4.count.txt + 1), 
     log10(his4_pool_reads$IVT_PreR_S51_L008_R1_001his4.count.txt + 1), pch=20, 
     col = alpha("black", 0.1), axes = FALSE, 
     type="n")
## Capture the plot window coordinates
coords <- par("usr")
gx <- grconvertX(coords[1:2], "user", "inches")
gy <- grconvertY(coords[3:4], "user", "inches")
width <- max(gx) - min(gx)
height <- max(gy) - min(gy)

pointsfile <- "his4_IVT_pre_replicates.png"
png(pointsfile, width=width, height=height, units="in", res=1500, bg="transparent")
par(mar=c(0,0,0,0))

plot(log10(his4_pool_reads$IVT_PreL_S50_L008_R1_001his4.count.txt + 1), 
     log10(his4_pool_reads$IVT_PreR_S51_L008_R1_001his4.count.txt + 1), pch=20, 
     col = alpha("black", 0.1), axes = FALSE)
dev.off()

## Read in the PNG file and draw as a raster image
panel <- readPNG(pointsfile, native = TRUE)
rasterImage(panel, coords[1], coords[3], coords[2], coords[4])

axis(side = 1, at = c(log10(1), log10(2), log10(3), log10(4), 
                      log10(5), log10(6), log10(7), log10(8), 
                      log10(9), log10(10), log10(20), log10(30), log10(40), 
                      log10(50), log10(60), log10(70), log10(80), 
                      log10(90), log10(100), log10(200), log10(300), log10(400), 
                      log10(500), log10(600), log10(700), log10(800), 
                      log10(900), log10(1000), log10(2000), log10(3000), log10(4000), log10(5000), log10(6000), log10(7000), log10(8000), log10(9000), log10(10000)), 
     labels=c("1", "2", "3", "4", "5", "6", "7", "8", "9",
              "10", "20", "30", "40", "50", "60", "70", "80", "90", 
              "100", "200", "300", "400", "500", "600", "700", "800", "900", "1000", 
              "2000", "3000", "4000", "5000", "6000", "7000", "8000", "9000", "10000"))
axis(side = 2, at = c(log10(1), log10(2), log10(3), log10(4), 
                      log10(5), log10(6), log10(7), log10(8), 
                      log10(9), log10(10), log10(20), log10(30), log10(40), 
                      log10(50), log10(60), log10(70), log10(80), 
                      log10(90), log10(100), log10(200), log10(300), log10(400), 
                      log10(500), log10(600), log10(700), log10(800), 
                      log10(900), log10(1000), log10(2000), log10(3000), log10(4000), log10(5000), log10(6000), log10(7000), log10(8000), log10(9000), log10(10000)), 
     labels=c("1", "2", "3", "4", "5", "6", "7", "8", "9",
              "10", "20", "30", "40", "50", "60", "70", "80", "90", 
              "100", "200", "300", "400", "500", "600", "700", "800", "900", "1000", 
              "2000", "3000", "4000", "5000", "6000", "7000", "8000", "9000", "10000"))
title(main = "his4_IVT_pre_replicates")
dev.off()

#his4 IVT post replicates__________________________________________________
pdf("his4_IVT_post_replicates.pdf", useDingbats=FALSE, width=5, height=5, colormodel="rgb")
par(pty="s")
plot(log10(his4_pool_reads$IVT_PostL_S52_L008_R1_001his4.count.txt + 1), 
     log10(his4_pool_reads$IVT_PostR_S53_L008_R1_001his4.count.txt + 1), pch=20, 
     col = alpha("black", 0.1), axes = FALSE, 
     type="n")
## Capture the plot window coordinates
coords <- par("usr")
gx <- grconvertX(coords[1:2], "user", "inches")
gy <- grconvertY(coords[3:4], "user", "inches")
width <- max(gx) - min(gx)
height <- max(gy) - min(gy)

pointsfile <- "his4_IVT_post_replicates.png"
png(pointsfile, width=width, height=height, units="in", res=1500, bg="transparent")
par(mar=c(0,0,0,0))

plot(log10(his4_pool_reads$IVT_PostL_S52_L008_R1_001his4.count.txt + 1), 
     log10(his4_pool_reads$IVT_PostR_S53_L008_R1_001his4.count.txt + 1), pch=20, 
     col = alpha("black", 0.1), axes = FALSE)
dev.off()

## Read in the PNG file and draw as a raster image
panel <- readPNG(pointsfile, native = TRUE)
rasterImage(panel, coords[1], coords[3], coords[2], coords[4])

axis(side = 1, at = c(log10(1), log10(2), log10(3), log10(4), 
                      log10(5), log10(6), log10(7), log10(8), 
                      log10(9), log10(10), log10(20), log10(30), log10(40), 
                      log10(50), log10(60), log10(70), log10(80), 
                      log10(90), log10(100), log10(200), log10(300), log10(400), 
                      log10(500), log10(600), log10(700), log10(800), 
                      log10(900), log10(1000), log10(2000), log10(3000), log10(4000), log10(5000), log10(6000), log10(7000), log10(8000), log10(9000), log10(10000)), 
     labels=c("1", "2", "3", "4", "5", "6", "7", "8", "9",
              "10", "20", "30", "40", "50", "60", "70", "80", "90", 
              "100", "200", "300", "400", "500", "600", "700", "800", "900", "1000", 
              "2000", "3000", "4000", "5000", "6000", "7000", "8000", "9000", "10000"))
axis(side = 2, at = c(log10(1), log10(2), log10(3), log10(4), 
                      log10(5), log10(6), log10(7), log10(8), 
                      log10(9), log10(10), log10(20), log10(30), log10(40), 
                      log10(50), log10(60), log10(70), log10(80), 
                      log10(90), log10(100), log10(200), log10(300), log10(400), 
                      log10(500), log10(600), log10(700), log10(800), 
                      log10(900), log10(1000), log10(2000), log10(3000), log10(4000), log10(5000), log10(6000), log10(7000), log10(8000), log10(9000), log10(10000)), 
     labels=c("1", "2", "3", "4", "5", "6", "7", "8", "9",
              "10", "20", "30", "40", "50", "60", "70", "80", "90", 
              "100", "200", "300", "400", "500", "600", "700", "800", "900", "1000", 
              "2000", "3000", "4000", "5000", "6000", "7000", "8000", "9000", "10000"))
title(main = "his4_IVT_post_replicates")
dev.off()

#pgk1 IVT pre replicates__________________________________________________
pdf("pgk1_IVT_pre_replicates.pdf", useDingbats=FALSE, width=5, height=5, colormodel="rgb")
par(pty="s")
plot(log10(pgk1_pool_reads$IVT_PreL_S50_L008_R1_001pgk1.count.txt + 1), 
     log10(pgk1_pool_reads$IVT_PreR_S51_L008_R1_001pgk1.count.txt + 1), pch=20, 
     col = alpha("black", 0.1), axes = FALSE, 
     type="n")
## Capture the plot window coordinates
coords <- par("usr")
gx <- grconvertX(coords[1:2], "user", "inches")
gy <- grconvertY(coords[3:4], "user", "inches")
width <- max(gx) - min(gx)
height <- max(gy) - min(gy)

pointsfile <- "pgk1_IVT_pre_replicates.png"
png(pointsfile, width=width, height=height, units="in", res=1500, bg="transparent")
par(mar=c(0,0,0,0))

plot(log10(pgk1_pool_reads$IVT_PreL_S50_L008_R1_001pgk1.count.txt + 1), 
     log10(pgk1_pool_reads$IVT_PreR_S51_L008_R1_001pgk1.count.txt + 1), pch=20, 
     col = alpha("black", 0.1), axes = FALSE)
dev.off()

## Read in the PNG file and draw as a raster image
panel <- readPNG(pointsfile, native = TRUE)
rasterImage(panel, coords[1], coords[3], coords[2], coords[4])

axis(side = 1, at = c(log10(1), log10(2), log10(3), log10(4), 
                      log10(5), log10(6), log10(7), log10(8), 
                      log10(9), log10(10), log10(20), log10(30), log10(40), 
                      log10(50), log10(60), log10(70), log10(80), 
                      log10(90), log10(100), log10(200), log10(300), log10(400), 
                      log10(500), log10(600), log10(700), log10(800), 
                      log10(900), log10(1000), log10(2000), log10(3000), log10(4000), log10(5000), log10(6000), log10(7000), log10(8000), log10(9000), log10(10000)), 
     labels=c("1", "2", "3", "4", "5", "6", "7", "8", "9",
              "10", "20", "30", "40", "50", "60", "70", "80", "90", 
              "100", "200", "300", "400", "500", "600", "700", "800", "900", "1000", 
              "2000", "3000", "4000", "5000", "6000", "7000", "8000", "9000", "10000"))
axis(side = 2, at = c(log10(1), log10(2), log10(3), log10(4), 
                      log10(5), log10(6), log10(7), log10(8), 
                      log10(9), log10(10), log10(20), log10(30), log10(40), 
                      log10(50), log10(60), log10(70), log10(80), 
                      log10(90), log10(100), log10(200), log10(300), log10(400), 
                      log10(500), log10(600), log10(700), log10(800), 
                      log10(900), log10(1000), log10(2000), log10(3000), log10(4000), log10(5000), log10(6000), log10(7000), log10(8000), log10(9000), log10(10000)), 
     labels=c("1", "2", "3", "4", "5", "6", "7", "8", "9",
              "10", "20", "30", "40", "50", "60", "70", "80", "90", 
              "100", "200", "300", "400", "500", "600", "700", "800", "900", "1000", 
              "2000", "3000", "4000", "5000", "6000", "7000", "8000", "9000", "10000"))
title(main = "pgk1_IVT_pre_replicates")
dev.off()


#pgk1 IVT post replicates__________________________________________________
pdf("pgk1_IVT_post_replicates.pdf", useDingbats=FALSE, width=5, height=5, colormodel="rgb")
par(pty="s")
plot(log10(pgk1_pool_reads$IVT_PostL_S52_L008_R1_001pgk1.count.txt + 1), 
     log10(pgk1_pool_reads$IVT_PostR_S53_L008_R1_001pgk1.count.txt + 1), pch=20, 
     col = alpha("black", 0.1), axes = FALSE, 
     type="n")
## Capture the plot window coordinates
coords <- par("usr")
gx <- grconvertX(coords[1:2], "user", "inches")
gy <- grconvertY(coords[3:4], "user", "inches")
width <- max(gx) - min(gx)
height <- max(gy) - min(gy)

pointsfile <- "pgk1_IVT_post_replicates.png"
png(pointsfile, width=width, height=height, units="in", res=1500, bg="transparent")
par(mar=c(0,0,0,0))

plot(log10(pgk1_pool_reads$IVT_PostL_S52_L008_R1_001pgk1.count.txt + 1), 
     log10(pgk1_pool_reads$IVT_PostR_S53_L008_R1_001pgk1.count.txt + 1), pch=20, 
     col = alpha("black", 0.1), axes = FALSE)
dev.off()

## Read in the PNG file and draw as a raster image
panel <- readPNG(pointsfile, native = TRUE)
rasterImage(panel, coords[1], coords[3], coords[2], coords[4])

axis(side = 1, at = c(log10(1), log10(2), log10(3), log10(4), 
                      log10(5), log10(6), log10(7), log10(8), 
                      log10(9), log10(10), log10(20), log10(30), log10(40), 
                      log10(50), log10(60), log10(70), log10(80), 
                      log10(90), log10(100), log10(200), log10(300), log10(400), 
                      log10(500), log10(600), log10(700), log10(800), 
                      log10(900), log10(1000), log10(2000), log10(3000), log10(4000), log10(5000), log10(6000), log10(7000), log10(8000), log10(9000), log10(10000)), 
     labels=c("1", "2", "3", "4", "5", "6", "7", "8", "9",
              "10", "20", "30", "40", "50", "60", "70", "80", "90", 
              "100", "200", "300", "400", "500", "600", "700", "800", "900", "1000", 
              "2000", "3000", "4000", "5000", "6000", "7000", "8000", "9000", "10000"))
axis(side = 2, at = c(log10(1), log10(2), log10(3), log10(4), 
                      log10(5), log10(6), log10(7), log10(8), 
                      log10(9), log10(10), log10(20), log10(30), log10(40), 
                      log10(50), log10(60), log10(70), log10(80), 
                      log10(90), log10(100), log10(200), log10(300), log10(400), 
                      log10(500), log10(600), log10(700), log10(800), 
                      log10(900), log10(1000), log10(2000), log10(3000), log10(4000), log10(5000), log10(6000), log10(7000), log10(8000), log10(9000), log10(10000)), 
     labels=c("1", "2", "3", "4", "5", "6", "7", "8", "9",
              "10", "20", "30", "40", "50", "60", "70", "80", "90", 
              "100", "200", "300", "400", "500", "600", "700", "800", "900", "1000", 
              "2000", "3000", "4000", "5000", "6000", "7000", "8000", "9000", "10000"))
title(main = "pgk1_IVT_post_replicates")
dev.off()

#his4 RNA pre replicates__________________________________________________
pdf("his4_RNA_pre_replicates.pdf", useDingbats=FALSE, width=5, height=5, colormodel="rgb")
par(pty="s")
plot(log10(his4_pool_reads$RNA_PreL_S54_L008_R1_001his4.count.txt + 1), 
     log10(his4_pool_reads$RNA_PreR_S55_L008_R1_001his4.count.txt + 1), pch=20, 
     col = alpha("black", 0.1), axes = FALSE, 
     type="n")
## Capture the plot window coordinates
coords <- par("usr")
gx <- grconvertX(coords[1:2], "user", "inches")
gy <- grconvertY(coords[3:4], "user", "inches")
width <- max(gx) - min(gx)
height <- max(gy) - min(gy)

pointsfile <- "his4_RNA_pre_replicates.png"
png(pointsfile, width=width, height=height, units="in", res=1500, bg="transparent")
par(mar=c(0,0,0,0))

plot(log10(his4_pool_reads$RNA_PreL_S54_L008_R1_001his4.count.txt + 1), 
     log10(his4_pool_reads$RNA_PreR_S55_L008_R1_001his4.count.txt + 1), pch=20, 
     col = alpha("black", 0.1), axes = FALSE)
dev.off()

## Read in the PNG file and draw as a raster image
panel <- readPNG(pointsfile, native = TRUE)
rasterImage(panel, coords[1], coords[3], coords[2], coords[4])

axis(side = 1, at = c(log10(1), log10(2), log10(3), log10(4), 
                      log10(5), log10(6), log10(7), log10(8), 
                      log10(9), log10(10), log10(20), log10(30), log10(40), 
                      log10(50), log10(60), log10(70), log10(80), 
                      log10(90), log10(100), log10(200), log10(300), log10(400), 
                      log10(500), log10(600), log10(700), log10(800), 
                      log10(900), log10(1000), log10(2000), log10(3000), log10(4000), log10(5000), log10(6000), log10(7000), log10(8000), log10(9000), log10(10000)), 
     labels=c("1", "2", "3", "4", "5", "6", "7", "8", "9",
              "10", "20", "30", "40", "50", "60", "70", "80", "90", 
              "100", "200", "300", "400", "500", "600", "700", "800", "900", "1000", 
              "2000", "3000", "4000", "5000", "6000", "7000", "8000", "9000", "10000"))
axis(side = 2, at = c(log10(1), log10(2), log10(3), log10(4), 
                      log10(5), log10(6), log10(7), log10(8), 
                      log10(9), log10(10), log10(20), log10(30), log10(40), 
                      log10(50), log10(60), log10(70), log10(80), 
                      log10(90), log10(100), log10(200), log10(300), log10(400), 
                      log10(500), log10(600), log10(700), log10(800), 
                      log10(900), log10(1000), log10(2000), log10(3000), log10(4000), log10(5000), log10(6000), log10(7000), log10(8000), log10(9000), log10(10000)), 
     labels=c("1", "2", "3", "4", "5", "6", "7", "8", "9",
              "10", "20", "30", "40", "50", "60", "70", "80", "90", 
              "100", "200", "300", "400", "500", "600", "700", "800", "900", "1000", 
              "2000", "3000", "4000", "5000", "6000", "7000", "8000", "9000", "10000"))
title(main = "his4_RNA_pre_replicates")
dev.off()

#his4 RNA post replicates__________________________________________________
pdf("his4_RNA_post_replicates.pdf", useDingbats=FALSE, width=5, height=5, colormodel="rgb")
par(pty="s")
plot(log10(his4_pool_reads$RNA_PostL_S56_L008_R1_001his4.count.txt + 1), 
     log10(his4_pool_reads$RNA_PostR_S57_L008_R1_001his4.count.txt + 1), pch=20, 
     col = alpha("black", 0.1), axes = FALSE, 
     type="n")
## Capture the plot window coordinates
coords <- par("usr")
gx <- grconvertX(coords[1:2], "user", "inches")
gy <- grconvertY(coords[3:4], "user", "inches")
width <- max(gx) - min(gx)
height <- max(gy) - min(gy)

pointsfile <- "his4_RNA_post_replicates.png"
png(pointsfile, width=width, height=height, units="in", res=1500, bg="transparent")
par(mar=c(0,0,0,0))

plot(log10(his4_pool_reads$RNA_PostL_S56_L008_R1_001his4.count.txt + 1), 
     log10(his4_pool_reads$RNA_PostR_S57_L008_R1_001his4.count.txt + 1), pch=20, 
     col = alpha("black", 0.1), axes = FALSE)
dev.off()

## Read in the PNG file and draw as a raster image
panel <- readPNG(pointsfile, native = TRUE)
rasterImage(panel, coords[1], coords[3], coords[2], coords[4])

axis(side = 1, at = c(log10(1), log10(2), log10(3), log10(4), 
                      log10(5), log10(6), log10(7), log10(8), 
                      log10(9), log10(10), log10(20), log10(30), log10(40), 
                      log10(50), log10(60), log10(70), log10(80), 
                      log10(90), log10(100), log10(200), log10(300), log10(400), 
                      log10(500), log10(600), log10(700), log10(800), 
                      log10(900), log10(1000), log10(2000), log10(3000), log10(4000), log10(5000), log10(6000), log10(7000), log10(8000), log10(9000), log10(10000)), 
     labels=c("1", "2", "3", "4", "5", "6", "7", "8", "9",
              "10", "20", "30", "40", "50", "60", "70", "80", "90", 
              "100", "200", "300", "400", "500", "600", "700", "800", "900", "1000", 
              "2000", "3000", "4000", "5000", "6000", "7000", "8000", "9000", "10000"))
axis(side = 2, at = c(log10(1), log10(2), log10(3), log10(4), 
                      log10(5), log10(6), log10(7), log10(8), 
                      log10(9), log10(10), log10(20), log10(30), log10(40), 
                      log10(50), log10(60), log10(70), log10(80), 
                      log10(90), log10(100), log10(200), log10(300), log10(400), 
                      log10(500), log10(600), log10(700), log10(800), 
                      log10(900), log10(1000), log10(2000), log10(3000), log10(4000), log10(5000), log10(6000), log10(7000), log10(8000), log10(9000), log10(10000)), 
     labels=c("1", "2", "3", "4", "5", "6", "7", "8", "9",
              "10", "20", "30", "40", "50", "60", "70", "80", "90", 
              "100", "200", "300", "400", "500", "600", "700", "800", "900", "1000", 
              "2000", "3000", "4000", "5000", "6000", "7000", "8000", "9000", "10000"))
title(main = "his4_RNA_post_replicates")
dev.off()

#pgk1 RNA pre replicates__________________________________________________
pdf("pgk1_RNA_pre_replicates.pdf", useDingbats=FALSE, width=5, height=5, colormodel="rgb")
par(pty="s")
plot(log10(pgk1_pool_reads$RNA_PreL_S54_L008_R1_001pgk1.count.txt + 1), 
     log10(pgk1_pool_reads$RNA_PreR_S55_L008_R1_001pgk1.count.txt + 1), pch=20, 
     col = alpha("black", 0.1), axes = FALSE, 
     type="n")
## Capture the plot window coordinates
coords <- par("usr")
gx <- grconvertX(coords[1:2], "user", "inches")
gy <- grconvertY(coords[3:4], "user", "inches")
width <- max(gx) - min(gx)
height <- max(gy) - min(gy)

pointsfile <- "pgk1_RNA_pre_replicates.png"
png(pointsfile, width=width, height=height, units="in", res=1500, bg="transparent")
par(mar=c(0,0,0,0))

plot(log10(pgk1_pool_reads$RNA_PreL_S54_L008_R1_001pgk1.count.txt + 1), 
     log10(pgk1_pool_reads$RNA_PreR_S55_L008_R1_001pgk1.count.txt + 1), pch=20, 
     col = alpha("black", 0.1), axes = FALSE)
dev.off()

## Read in the PNG file and draw as a raster image
panel <- readPNG(pointsfile, native = TRUE)
rasterImage(panel, coords[1], coords[3], coords[2], coords[4])

axis(side = 1, at = c(log10(1), log10(2), log10(3), log10(4), 
                      log10(5), log10(6), log10(7), log10(8), 
                      log10(9), log10(10), log10(20), log10(30), log10(40), 
                      log10(50), log10(60), log10(70), log10(80), 
                      log10(90), log10(100), log10(200), log10(300), log10(400), 
                      log10(500), log10(600), log10(700), log10(800), 
                      log10(900), log10(1000), log10(2000), log10(3000), log10(4000), log10(5000), log10(6000), log10(7000), log10(8000), log10(9000), log10(10000)), 
     labels=c("1", "2", "3", "4", "5", "6", "7", "8", "9",
              "10", "20", "30", "40", "50", "60", "70", "80", "90", 
              "100", "200", "300", "400", "500", "600", "700", "800", "900", "1000", 
              "2000", "3000", "4000", "5000", "6000", "7000", "8000", "9000", "10000"))
axis(side = 2, at = c(log10(1), log10(2), log10(3), log10(4), 
                      log10(5), log10(6), log10(7), log10(8), 
                      log10(9), log10(10), log10(20), log10(30), log10(40), 
                      log10(50), log10(60), log10(70), log10(80), 
                      log10(90), log10(100), log10(200), log10(300), log10(400), 
                      log10(500), log10(600), log10(700), log10(800), 
                      log10(900), log10(1000), log10(2000), log10(3000), log10(4000), log10(5000), log10(6000), log10(7000), log10(8000), log10(9000), log10(10000)), 
     labels=c("1", "2", "3", "4", "5", "6", "7", "8", "9",
              "10", "20", "30", "40", "50", "60", "70", "80", "90", 
              "100", "200", "300", "400", "500", "600", "700", "800", "900", "1000", 
              "2000", "3000", "4000", "5000", "6000", "7000", "8000", "9000", "10000"))
title(main = "pgk1_RNA_pre_replicates")
dev.off()


#pgk1 RNA post replicates__________________________________________________
pdf("pgk1_RNA_post_replicates.pdf", useDingbats=FALSE, width=5, height=5, colormodel="rgb")
par(pty="s")
plot(log10(pgk1_pool_reads$RNA_PostL_S56_L008_R1_001pgk1.count.txt + 1), 
     log10(pgk1_pool_reads$RNA_PostR_S57_L008_R1_001pgk1.count.txt + 1), pch=20, 
     col = alpha("black", 0.1), axes = FALSE, 
     type="n")
## Capture the plot window coordinates
coords <- par("usr")
gx <- grconvertX(coords[1:2], "user", "inches")
gy <- grconvertY(coords[3:4], "user", "inches")
width <- max(gx) - min(gx)
height <- max(gy) - min(gy)

pointsfile <- "pgk1_RNA_post_replicates.png"
png(pointsfile, width=width, height=height, units="in", res=1500, bg="transparent")
par(mar=c(0,0,0,0))

plot(log10(pgk1_pool_reads$RNA_PostL_S56_L008_R1_001pgk1.count.txt + 1), 
     log10(pgk1_pool_reads$RNA_PostR_S57_L008_R1_001pgk1.count.txt + 1), pch=20, 
     col = alpha("black", 0.1), axes = FALSE)
dev.off()

## Read in the PNG file and draw as a raster image
panel <- readPNG(pointsfile, native = TRUE)
rasterImage(panel, coords[1], coords[3], coords[2], coords[4])

axis(side = 1, at = c(log10(1), log10(2), log10(3), log10(4), 
                      log10(5), log10(6), log10(7), log10(8), 
                      log10(9), log10(10), log10(20), log10(30), log10(40), 
                      log10(50), log10(60), log10(70), log10(80), 
                      log10(90), log10(100), log10(200), log10(300), log10(400), 
                      log10(500), log10(600), log10(700), log10(800), 
                      log10(900), log10(1000), log10(2000), log10(3000), log10(4000), log10(5000), log10(6000), log10(7000), log10(8000), log10(9000), log10(10000)), 
     labels=c("1", "2", "3", "4", "5", "6", "7", "8", "9",
              "10", "20", "30", "40", "50", "60", "70", "80", "90", 
              "100", "200", "300", "400", "500", "600", "700", "800", "900", "1000", 
              "2000", "3000", "4000", "5000", "6000", "7000", "8000", "9000", "10000"))
axis(side = 2, at = c(log10(1), log10(2), log10(3), log10(4), 
                      log10(5), log10(6), log10(7), log10(8), 
                      log10(9), log10(10), log10(20), log10(30), log10(40), 
                      log10(50), log10(60), log10(70), log10(80), 
                      log10(90), log10(100), log10(200), log10(300), log10(400), 
                      log10(500), log10(600), log10(700), log10(800), 
                      log10(900), log10(1000), log10(2000), log10(3000), log10(4000), log10(5000), log10(6000), log10(7000), log10(8000), log10(9000), log10(10000)), 
     labels=c("1", "2", "3", "4", "5", "6", "7", "8", "9",
              "10", "20", "30", "40", "50", "60", "70", "80", "90", 
              "100", "200", "300", "400", "500", "600", "700", "800", "900", "1000", 
              "2000", "3000", "4000", "5000", "6000", "7000", "8000", "9000", "10000"))
title(main = "pgk1_RNA_post_replicates")
dev.off()


