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

GO_count <- data.frame("firstrow", 0, 0, 0)
names(GO_count) <- c("GOsample", "total","up", "dwn")
mpralmpath <- "/mnt/ingolialab/rmuller1/test_sandbox/gz_to_mpralm/rehearsal_dinner/all_raw_fasta_gz/"
GOpath <- "/mnt/ingolialab/rmuller1/test_sandbox/gz_to_mpralm/rehearsal_dinner/scripts/GO_analysis_files/GO_annotation_lists/"
sgd_yorf <- data.frame(sgd$sgdid, sgd$name)
head(sgd_yorf)
options(stringsAsFactors=FALSE)
for (annot in c("cellular_amino_acid_metabolic_process", "nucleic_acid_metabolic_process",
                "RNA_processing", "tRNA_transcription", "cellular_component_biogenesis",
                "proteasomal_ubiquitin_independent_protein_catabolic_process",
                "RNP_complex_biogenesis", "vesicle_fusion_with_Golgi_apparatus", "gene_expression",
                "ribosome_biogenesis", "tRNA_aminoacylation", "glycolytic_process",
                "RNA_metabolic_process", "tRNA_metabolic_process", "nuclear_pore_distribution",
                "protein_sumoylation", "regulation_of_translation")){
  x <- read.delim(paste(GOpath, annot, ".txt", sep=""),
                  stringsAsFactors=FALSE, header = FALSE)
  x <- data.frame(x$V1, x$V2)
  x$sgd.sgdid <- gsub("SGD:", "", x$x.V1)  
  z <- merge(x, sgd_yorf, by = "sgd.sgdid")
  colnames(z)[colnames(z)=="sgd.name"] <- "Yorf1"
  assign(paste(annot, sep=""), z)
  for (dir_file in c("HIS4_PGK1_3AT/his4_postv3AT_sum_mpralm", "HIS4_PGK1_3AT/pgk1_postv3AT_sum_mpralm",
                     "HIS4_PGK1_pooled/his4_pooled_sum_mpralm", "HIS4_PGK1_pooled/pgk1_pooled_sum_mpralm",
                     "HTS1_RPC31/HTS1_sum_mpralm", "HTS1_RPC31/RPC31_sum_mpralm",
                     "GCN4_CDS_UTR/cds_sum_mpralm", "GCN4_CDS_UTR/utr_sum_mpralm")){
    y <- read.delim(paste(mpralmpath, dir_file, ".txt", sep=""),
                    stringsAsFactors=FALSE, header = TRUE)
    temp1 <- gsub(".*/", "", dir_file)
    temp3 <- temp2 <- gsub("_sum_mpralm", "", temp1)
    assign(paste(temp2, sep=""), y)
    m <- merge(z, y, by = "Yorf1")
    m <- unique(m, incomparables = FALSE)
    assign(paste(annot, temp3, sep=""), m)
    m_up <- filter(m, m$logFC < -0.5 & m$adj.P.Val < 0.05)
    m_dwn <- filter(m, m$logFC > 0.5 & m$adj.P.Val < 0.05)
    assign(paste(annot, temp3, "sig_up", sep=""), m_up)
    assign(paste(annot, temp3, "sig_dwn", sep=""), m_dwn)
    y_up <- filter(y, y$logFC < -0.5 & y$adj.P.Val < 0.05)
    y_dwn <- filter(y, y$logFC > 0.5 & y$adj.P.Val < 0.05)
    assign(paste(temp3, "sig_up", sep=""), y_up)
    assign(paste(temp3, "sig_dwn", sep=""), y_dwn)
    full_mpralm <- data.frame(temp3, nrow(y), nrow(y_up), nrow(y_dwn))
    names(full_mpralm) <- c("GOsample", "total","up", "dwn")
    GO_mpralm <- data.frame(paste(annot, temp3, sep=""), nrow(m), nrow(m_up), nrow(m_dwn))
    names(GO_mpralm) <- c("GOsample", "total","up", "dwn")
    GO_count <- rbind(GO_count, full_mpralm)
    GO_count <- rbind(GO_count, GO_mpralm)
  }
}

GO_count <- unique(GO_count, incomparables = FALSE)
write.table(GO_count, "/mnt/ingolialab/rmuller1/test_sandbox/gz_to_mpralm/rehearsal_dinner/scripts/GO_analysis_files/GO_count.txt")

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
temp1 <- pmax(-6, his4_pooled$logFC)
temp2 <- pmin(6, temp1)
his4_pooled$LFC_plot <- -temp2
his4_pooled$sig_plot <- -log10(pmax(10^-20, his4_pooled$adj.P.Val))

pdf("his4_pooled_volcano.pdf", useDingbats=FALSE, width=5, height=5, colormodel="rgb")
## Plot with type="n" to skip drawing the points
plot(x = his4_pooled$LFC_plot, y = his4_pooled$sig_plot,
     type="n",
     xlim=c(-6, 6), ylim=c(0, 20),
     axes=FALSE,
     xlab="Log2 Fold Before vs After Guide tet Induction", 
     ylab = expression('Statistical Significance -log10(Q)'), 
     col="gray52")
points(-nucleic_acid_metabolic_processhis4_pooled$logFC,
       -log10(pmax(10^-20, nucleic_acid_metabolic_processhis4_pooled$adj.P.Val)),
       type="n",
       col="saddlebrown"
)
points(-RNA_metabolic_processhis4_pooled$logFC,
       -log10(pmax(10^-20, RNA_metabolic_processhis4_pooled$adj.P.Val)),
       type="n",
       col="darksalmon"
)
points(-cellular_amino_acid_metabolic_processhis4_pooled$logFC,
       -log10(pmax(10^-20, cellular_amino_acid_metabolic_processhis4_pooled$adj.P.Val)),
       type="n",
       col="violet"
)
points(-tRNA_transcriptionhis4_pooled$logFC,
       -log10(pmax(10^-20, tRNA_transcriptionhis4_pooled$adj.P.Val)),
       type="n",
       col="palegreen3"
)
points(-tRNA_metabolic_processhis4_pooled$logFC,
       -log10(pmax(10^-20, tRNA_metabolic_processhis4_pooled$adj.P.Val)),
       type="n",
       col="orange1"
)
points(-tRNA_aminoacylationhis4_pooled$logFC,
       -log10(pmax(10^-20, tRNA_aminoacylationhis4_pooled$adj.P.Val)),
       type="n",
       col="orangered1"
)
points(-regulation_of_translationhis4_pooled$logFC,
       -log10(pmax(10^-20, regulation_of_translationhis4_pooled$adj.P.Val)),
       type="n",
       col="royalblue1"
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
pointsfile <- "his4_pooled_volcano.png"
png(pointsfile, width=width, height=height, units="in", res=1500, bg="transparent")
par(mar=c(0,0,0,0))
## Exactly the same plot as above, but no type="n" and instead specify points, color, etc.
plot(x = his4_pooled$LFC_plot, y = his4_pooled$sig_plot,
     pch=16, cex=0.9,
     xlim=c(-6, 6), ylim=c(0, 20),
     axes=FALSE,
     xlab="Log2 Fold Before vs After Guide tet Induction", 
     ylab = expression('Statistical Significance -log10(Q)'), 
     col = alpha("gray52", 1.0))
points(-nucleic_acid_metabolic_processhis4_pooled$logFC,
       -log10(pmax(10^-20, nucleic_acid_metabolic_processhis4_pooled$adj.P.Val)),
       pch=16, cex=0.9,
       col="saddlebrown"
)
points(-RNA_metabolic_processhis4_pooled$logFC,
       -log10(pmax(10^-20, RNA_metabolic_processhis4_pooled$adj.P.Val)),
       pch=16, cex=0.9,
       col="darksalmon"
)
points(-cellular_amino_acid_metabolic_processhis4_pooled$logFC,
       -log10(pmax(10^-20, cellular_amino_acid_metabolic_processhis4_pooled$adj.P.Val)),
       pch=16, cex=0.9,
       col = alpha("mediumorchid2", 1.0)
)
points(-tRNA_transcriptionhis4_pooled$logFC,
       -log10(pmax(10^-20, tRNA_transcriptionhis4_pooled$adj.P.Val)),
       pch=16, cex=0.9,
       col = alpha("palegreen4", 1.0)
)
points(-tRNA_metabolic_processhis4_pooled$logFC,
       -log10(pmax(10^-20, tRNA_metabolic_processhis4_pooled$adj.P.Val)),
       pch=16, cex=0.9,
       col = alpha("orange1", 1.0)
)
points(-tRNA_aminoacylationhis4_pooled$logFC,
       -log10(pmax(10^-20, tRNA_aminoacylationhis4_pooled$adj.P.Val)),
       pch=16, cex=0.9,
       col = alpha("orangered1", 1.0)
)
points(-regulation_of_translationhis4_pooled$logFC,
       -log10(pmax(10^-20, regulation_of_translationhis4_pooled$adj.P.Val)),
       pch=16, cex=0.9,
       col="royalblue1"
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
axis(side = 1, at = c(-6, -4, -2, 0, 2, 4, 6), labels=c("<1/64", "1/16", "1/4", "1", "4", "16", ">64"))
axis(side = 2, at = c(0, 5, 10, 15, 20), labels=c("0", "5", "10", "15", ">20"))
title(main="P(His4) Guide Induction")
dev.off()



#Pgk1 pooled before vs after guide tet induction treatment_________________

## Set up values to plot
pgk1_pooled_insig <- filter(pgk1_pooled, pgk1_pooled$adj.P.Val > 0.05)
temp1 <- pmax(-6, pgk1_pooled$logFC)
temp2 <- pmin(6, temp1)
pgk1_pooled$LFC_plot <- -temp2
pgk1_pooled$sig_plot <- -log10(pmax(10^-20, pgk1_pooled$adj.P.Val))

pdf("pgk1_pooled_volcano.pdf", useDingbats=FALSE, width=5, height=5, colormodel="rgb")
## Plot with type="n" to skip drawing the points
plot(x = pgk1_pooled$LFC_plot, y = pgk1_pooled$sig_plot,
     type="n",
     xlim=c(-6, 6), ylim=c(0, 20),
     axes=FALSE,
     xlab="Log2 Fold Before vs After Guide tet Induction", 
     ylab = expression('Statistical Significance -log10(Q)'), 
     col="gray52")
points(-RNA_metabolic_processpgk1_pooled$logFC,
       -log10(pmax(10^-20, RNA_metabolic_processpgk1_pooled$adj.P.Val)),
       type="n",
       col="darksalmon"
)
points(-glycolytic_processpgk1_pooled$logFC,
       -log10(pmax(10^-20, glycolytic_processpgk1_pooled$adj.P.Val)),
       type="n",
       col="darkolivegreen1"
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
pointsfile <- "pgk1_pooled_volcano.png"
png(pointsfile, width=width, height=height, units="in", res=1500, bg="transparent")
par(mar=c(0,0,0,0))
## Exactly the same plot as above, but no type="n" and instead specify points, color, etc.
plot(x = pgk1_pooled$LFC_plot, y = pgk1_pooled$sig_plot,
     pch=16, cex=0.9,
     xlim=c(-6, 6), ylim=c(0, 20),
     axes=FALSE,
     xlab="Log2 Fold Before vs After Guide tet Induction", 
     ylab = expression('Statistical Significance -log10(Q)'), 
     col = alpha("gray52", 1.0))
points(-RNA_metabolic_processpgk1_pooled$logFC,
       -log10(pmax(10^-20, RNA_metabolic_processpgk1_pooled$adj.P.Val)),
       pch=16, cex=0.9,
       col="darksalmon"
)
points(-glycolytic_processpgk1_pooled$logFC,
       -log10(pmax(10^-20, glycolytic_processpgk1_pooled$adj.P.Val)),
       pch=16, cex=0.9,
       col="darkolivegreen1"
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
axis(side = 1, at = c(-6, -4, -2, 0, 2, 4, 6), labels=c("<1/64", "1/16", "1/4", "1", "4", "16", ">64"))
axis(side = 2, at = c(0, 5, 10, 15, 20), labels=c("0", "5", "10", "15", ">20"))
title(main="P(Pgk1) Guide Induction")
dev.off()



#His4 3AT treatment_________________

## Set up values to plot
his4_postv3AT_insig <- filter(his4_postv3AT, his4_postv3AT$adj.P.Val > 0.05)
temp1 <- pmax(-6, his4_postv3AT$logFC)
temp2 <- pmin(6, temp1)
his4_postv3AT$LFC_plot <- -temp2
his4_postv3AT$sig_plot <- -log10(pmax(10^-20, his4_postv3AT$adj.P.Val))

pdf("his4_postv3AT_volcano.pdf", useDingbats=FALSE, width=5, height=5, colormodel="rgb")
## Plot with type="n" to skip drawing the points
plot(x = his4_postv3AT$LFC_plot, y = his4_postv3AT$sig_plot,
     type="n",
     xlim=c(-6, 6), ylim=c(0, 20),
     axes=FALSE,
     xlab="Log2 Fold Before vs After Guide tet Induction", 
     ylab = expression('Statistical Significance -log10(Q)'), 
     col="gray52")
points(-cellular_amino_acid_metabolic_processhis4_postv3AT$logFC,
       -log10(pmax(10^-20, cellular_amino_acid_metabolic_processhis4_postv3AT$adj.P.Val)),
       type="n",
       col="violet"
)
points(-tRNA_transcriptionhis4_postv3AT$logFC,
       -log10(pmax(10^-20, tRNA_transcriptionhis4_postv3AT$adj.P.Val)),
       type="n",
       col="palegreen3"
)
points(-tRNA_metabolic_processhis4_postv3AT$logFC,
       -log10(pmax(10^-20, tRNA_metabolic_processhis4_postv3AT$adj.P.Val)),
       type="n",
       col="orange1"
)
points(-tRNA_aminoacylationhis4_postv3AT$logFC,
       -log10(pmax(10^-20, tRNA_aminoacylationhis4_postv3AT$adj.P.Val)),
       type="n",
       col="orangered1"
)
points(-regulation_of_translationhis4_postv3AT$logFC,
       -log10(pmax(10^-20, regulation_of_translationhis4_postv3AT$adj.P.Val)),
       type="n",
       col="royalblue1"
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
pointsfile <- "his4_postv3AT_volcano.png"
png(pointsfile, width=width, height=height, units="in", res=1500, bg="transparent")
par(mar=c(0,0,0,0))
## Exactly the same plot as above, but no type="n" and instead specify points, color, etc.
plot(x = his4_postv3AT$LFC_plot, y = his4_postv3AT$sig_plot,
     pch=16, cex=0.9,
     xlim=c(-6, 6), ylim=c(0, 20),
     axes=FALSE,
     xlab="Log2 Fold Before vs After Guide tet Induction", 
     ylab = expression('Statistical Significance -log10(Q)'), 
     col = alpha("gray52", 1.0))
points(-cellular_amino_acid_metabolic_processhis4_postv3AT$logFC,
       -log10(pmax(10^-20, cellular_amino_acid_metabolic_processhis4_postv3AT$adj.P.Val)),
       pch=16, cex=0.9,
       col = alpha("mediumorchid2", 1.0)
)
points(-tRNA_transcriptionhis4_postv3AT$logFC,
       -log10(pmax(10^-20, tRNA_transcriptionhis4_postv3AT$adj.P.Val)),
       pch=16, cex=0.9,
       col = alpha("palegreen4", 1.0)
)
points(-tRNA_metabolic_processhis4_postv3AT$logFC,
       -log10(pmax(10^-20, tRNA_metabolic_processhis4_postv3AT$adj.P.Val)),
       pch=16, cex=0.9,
       col = alpha("orange1", 1.0)
)
points(-tRNA_aminoacylationhis4_postv3AT$logFC,
       -log10(pmax(10^-20, tRNA_aminoacylationhis4_postv3AT$adj.P.Val)),
       pch=16, cex=0.9,
       col = alpha("orangered1", 1.0)
)
points(-regulation_of_translationhis4_postv3AT$logFC,
       -log10(pmax(10^-20, regulation_of_translationhis4_postv3AT$adj.P.Val)),
       pch=16, cex=0.9,
       col="royalblue1"
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
axis(side = 1, at = c(-6, -4, -2, 0, 2, 4, 6), labels=c("<1/64", "1/16", "1/4", "1", "4", "16", ">64"))
axis(side = 2, at = c(0, 5, 10, 15, 20), labels=c("0", "5", "10", "15", ">20"))
title(main="P(His4) 3AT Treatment")
dev.off()



#Pgk1 3AT treatment_________________

## Set up values to plot
pgk1_postv3AT_insig <- filter(pgk1_postv3AT, pgk1_postv3AT$adj.P.Val > 0.05)
temp1 <- pmax(-6, pgk1_postv3AT$logFC)
temp2 <- pmin(6, temp1)
pgk1_postv3AT$LFC_plot <- -temp2
pgk1_postv3AT$sig_plot <- -log10(pmax(10^-20, pgk1_postv3AT$adj.P.Val))

pdf("pgk1_postv3AT_volcano.pdf", useDingbats=FALSE, width=5, height=5, colormodel="rgb")
## Plot with type="n" to skip drawing the points
plot(x = pgk1_postv3AT$LFC_plot, y = pgk1_postv3AT$sig_plot,
     type="n",
     xlim=c(-6, 6), ylim=c(0, 20),
     axes=FALSE,
     xlab="Log2 Fold Before vs After Guide tet Induction", 
     ylab = expression('Statistical Significance -log10(Q)'), 
     col="gray52")
points(-cellular_amino_acid_metabolic_processpgk1_postv3AT$logFC,
       -log10(pmax(10^-20, cellular_amino_acid_metabolic_processpgk1_postv3AT$adj.P.Val)),
       type="n",
       col="violet"
)
points(-tRNA_transcriptionpgk1_postv3AT$logFC,
       -log10(pmax(10^-20, tRNA_transcriptionpgk1_postv3AT$adj.P.Val)),
       type="n",
       col="palegreen3"
)
points(-tRNA_metabolic_processpgk1_postv3AT$logFC,
       -log10(pmax(10^-20, tRNA_metabolic_processpgk1_postv3AT$adj.P.Val)),
       type="n",
       col="orange1"
)
points(-tRNA_aminoacylationpgk1_postv3AT$logFC,
       -log10(pmax(10^-20, tRNA_aminoacylationpgk1_postv3AT$adj.P.Val)),
       type="n",
       col="orangered1"
)
points(-regulation_of_translationpgk1_postv3AT$logFC,
       -log10(pmax(10^-20, regulation_of_translationpgk1_postv3AT$adj.P.Val)),
       type="n",
       col="royalblue1"
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
pointsfile <- "pgk1_postv3AT_volcano.png"
png(pointsfile, width=width, height=height, units="in", res=1500, bg="transparent")
par(mar=c(0,0,0,0))
## Exactly the same plot as above, but no type="n" and instead specify points, color, etc.
plot(x = pgk1_postv3AT$LFC_plot, y = pgk1_postv3AT$sig_plot,
     pch=16, cex=0.9,
     xlim=c(-6, 6), ylim=c(0, 20),
     axes=FALSE,
     xlab="Log2 Fold Before vs After Guide tet Induction", 
     ylab = expression('Statistical Significance -log10(Q)'), 
     col = alpha("gray52", 1.0))
points(-cellular_amino_acid_metabolic_processpgk1_postv3AT$logFC,
       -log10(pmax(10^-20, cellular_amino_acid_metabolic_processpgk1_postv3AT$adj.P.Val)),
       pch=16, cex=0.9,
       col = alpha("mediumorchid2", 1.0)
)
points(-tRNA_transcriptionpgk1_postv3AT$logFC,
       -log10(pmax(10^-20, tRNA_transcriptionpgk1_postv3AT$adj.P.Val)),
       pch=16, cex=0.9,
       col = alpha("palegreen4", 1.0)
)
points(-tRNA_metabolic_processpgk1_postv3AT$logFC,
       -log10(pmax(10^-20, tRNA_metabolic_processpgk1_postv3AT$adj.P.Val)),
       pch=16, cex=0.9,
       col = alpha("orange1", 1.0)
)
points(-tRNA_aminoacylationpgk1_postv3AT$logFC,
       -log10(pmax(10^-20, tRNA_aminoacylationpgk1_postv3AT$adj.P.Val)),
       pch=16, cex=0.9,
       col = alpha("orangered1", 1.0)
)
points(-regulation_of_translationpgk1_postv3AT$logFC,
       -log10(pmax(10^-20, regulation_of_translationpgk1_postv3AT$adj.P.Val)),
       pch=16, cex=0.9,
       col="royalblue1"
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
axis(side = 1, at = c(-6, -4, -2, 0, 2, 4, 6), labels=c("<1/64", "1/16", "1/4", "1", "4", "16", ">64"))
axis(side = 2, at = c(0, 5, 10, 15, 20), labels=c("0", "5", "10", "15", ">20"))
title(main="P(Pgk1) 3AT Treatment")
dev.off()



#HTS1 guide induction_________________

## Set up values to plot
HTS1_insig <- filter(HTS1, HTS1$adj.P.Val > 0.05)
temp1 <- pmax(-6, HTS1$logFC)
temp2 <- pmin(6, temp1)
HTS1$LFC_plot <- -temp2
HTS1$sig_plot <- -log10(pmax(10^-20, HTS1$adj.P.Val))

pdf("HTS1_volcano.pdf", useDingbats=FALSE, width=5, height=5, colormodel="rgb")
## Plot with type="n" to skip drawing the points
plot(x = HTS1$LFC_plot, y = HTS1$sig_plot,
     type="n",
     xlim=c(-6, 6), ylim=c(0, 20),
     axes=FALSE,
     xlab="Log2 Fold Before vs After Guide tet Induction", 
     ylab = expression('Statistical Significance -log10(Q)'), 
     col="gray52")
points(-ribosome_biogenesisHTS1$logFC,
       -log10(pmax(10^-20, ribosome_biogenesisHTS1$adj.P.Val)),
       type="n",
       col="orange1"
)
points(-proteasomal_ubiquitin_independent_protein_catabolic_processHTS1$logFC,
       -log10(pmax(10^-20, proteasomal_ubiquitin_independent_protein_catabolic_processHTS1$adj.P.Val)),
       type="n",
       col="violet"
)
points(-vesicle_fusion_with_Golgi_apparatusHTS1$logFC,
       -log10(pmax(10^-20, vesicle_fusion_with_Golgi_apparatusHTS1$adj.P.Val)),
       type="n",
       col="palegreen3"
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
pointsfile <- "HTS1_volcano.png"
png(pointsfile, width=width, height=height, units="in", res=1500, bg="transparent")
par(mar=c(0,0,0,0))
## Exactly the same plot as above, but no type="n" and instead specify points, color, etc.
plot(x = HTS1$LFC_plot, y = HTS1$sig_plot,
     pch=16, cex=0.9,
     xlim=c(-6, 6), ylim=c(0, 20),
     axes=FALSE,
     xlab="Log2 Fold Before vs After Guide tet Induction", 
     ylab = expression('Statistical Significance -log10(Q)'), 
     col = alpha("gray52", 1.0))
points(-ribosome_biogenesisHTS1$logFC,
       -log10(pmax(10^-20, ribosome_biogenesisHTS1$adj.P.Val)),
       pch=16, cex=0.9,
       col="orange1"
)
points(-proteasomal_ubiquitin_independent_protein_catabolic_processHTS1$logFC,
       -log10(pmax(10^-20, proteasomal_ubiquitin_independent_protein_catabolic_processHTS1$adj.P.Val)),
       pch=16, cex=0.9,
       col="violet"
)
points(-vesicle_fusion_with_Golgi_apparatusHTS1$logFC,
       -log10(pmax(10^-20, vesicle_fusion_with_Golgi_apparatusHTS1$adj.P.Val)),
       pch=16, cex=0.9,
       col="palegreen3"
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
axis(side = 1, at = c(-6, -4, -2, 0, 2, 4, 6), labels=c("<1/64", "1/16", "1/4", "1", "4", "16", ">64"))
axis(side = 2, at = c(0, 5, 10, 15, 20), labels=c("0", "5", "10", "15", ">20"))
title(main="P(His4) Guide Induction in HTS1 Knockdown background")
dev.off()



#RPC31 guide induction_________________

## Set up values to plot
RPC31_insig <- filter(RPC31, RPC31$adj.P.Val > 0.05)
temp1 <- pmax(-6, RPC31$logFC)
temp2 <- pmin(6, temp1)
RPC31$LFC_plot <- -temp2
RPC31$sig_plot <- -log10(pmax(10^-20, RPC31$adj.P.Val))

pdf("RPC31_volcano.pdf", useDingbats=FALSE, width=5, height=5, colormodel="rgb")
## Plot with type="n" to skip drawing the points
plot(x = RPC31$LFC_plot, y = RPC31$sig_plot,
     type="n",
     xlim=c(-6, 6), ylim=c(0, 20),
     axes=FALSE,
     xlab="Log2 Fold Before vs After Guide tet Induction", 
     ylab = expression('Statistical Significance -log10(Q)'), 
     col="gray52")
points(-ribosome_biogenesisRPC31$logFC,
       -log10(pmax(10^-20, ribosome_biogenesisRPC31$adj.P.Val)),
       type="n",
       col="orange1"
)
points(-proteasomal_ubiquitin_independent_protein_catabolic_processRPC31$logFC,
       -log10(pmax(10^-20, proteasomal_ubiquitin_independent_protein_catabolic_processRPC31$adj.P.Val)),
       type="n",
       col="violet"
)
points(-vesicle_fusion_with_Golgi_apparatusRPC31$logFC,
       -log10(pmax(10^-20, vesicle_fusion_with_Golgi_apparatusRPC31$adj.P.Val)),
       type="n",
       col="palegreen3"
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
pointsfile <- "RPC31_volcano.png"
png(pointsfile, width=width, height=height, units="in", res=1500, bg="transparent")
par(mar=c(0,0,0,0))
## Exactly the same plot as above, but no type="n" and instead specify points, color, etc.
plot(x = RPC31$LFC_plot, y = RPC31$sig_plot,
     pch=16, cex=0.9,
     xlim=c(-6, 6), ylim=c(0, 20),
     axes=FALSE,
     xlab="Log2 Fold Before vs After Guide tet Induction", 
     ylab = expression('Statistical Significance -log10(Q)'), 
     col = alpha("gray52", 1.0))
points(-ribosome_biogenesisRPC31$logFC,
       -log10(pmax(10^-20, ribosome_biogenesisRPC31$adj.P.Val)),
       pch=16, cex=0.9,
       col="orange1"
)
points(-proteasomal_ubiquitin_independent_protein_catabolic_processRPC31$logFC,
       -log10(pmax(10^-20, proteasomal_ubiquitin_independent_protein_catabolic_processRPC31$adj.P.Val)),
       pch=16, cex=0.9,
       col="violet"
)
points(-vesicle_fusion_with_Golgi_apparatusRPC31$logFC,
       -log10(pmax(10^-20, vesicle_fusion_with_Golgi_apparatusRPC31$adj.P.Val)),
       pch=16, cex=0.9,
       col="palegreen3"
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
axis(side = 1, at = c(-6, -4, -2, 0, 2, 4, 6), labels=c("<1/64", "1/16", "1/4", "1", "4", "16", ">64"))
axis(side = 2, at = c(0, 5, 10, 15, 20), labels=c("0", "5", "10", "15", ">20"))
title(main="P(His4) Guide Induction in RPC31 Knockdown background")
dev.off()


#TF-cds guide induction_________________

## Set up values to plot
cds_insig <- filter(cds, cds$adj.P.Val > 0.05)
temp1 <- pmax(-6, cds$logFC)
temp2 <- pmin(6, temp1)
cds$LFC_plot <- -temp2
cds$sig_plot <- -log10(pmax(10^-20, cds$adj.P.Val))

pdf("cds_volcano.pdf", useDingbats=FALSE, width=5, height=5, colormodel="rgb")
## Plot with type="n" to skip drawing the points
plot(x = cds$LFC_plot, y = cds$sig_plot,
     type="n",
     xlim=c(-6, 6), ylim=c(0, 20),
     axes=FALSE,
     xlab="Log2 Fold Before vs After Guide tet Induction", 
     ylab = expression('Statistical Significance -log10(Q)'), 
     col="gray52")
points(-RNA_metabolic_processcds$logFC,
       -log10(pmax(10^-20, RNA_metabolic_processcds$adj.P.Val)),
       type="n",
       col="orange1"
)
points(-nuclear_pore_distributioncds$logFC,
       -log10(pmax(10^-20, nuclear_pore_distributioncds$adj.P.Val)),
       type="n",
       col="violet"
)
points(-protein_sumoylationcds$logFC,
       -log10(pmax(10^-20, protein_sumoylationcds$adj.P.Val)),
       type="n",
       col="palegreen3"
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
pointsfile <- "cds_volcano.png"
png(pointsfile, width=width, height=height, units="in", res=1500, bg="transparent")
par(mar=c(0,0,0,0))
## Exactly the same plot as above, but no type="n" and instead specify points, color, etc.
plot(x = cds$LFC_plot, y = cds$sig_plot,
     pch=16, cex=0.9,
     xlim=c(-6, 6), ylim=c(0, 20),
     axes=FALSE,
     xlab="Log2 Fold Before vs After Guide tet Induction", 
     ylab = expression('Statistical Significance -log10(Q)'), 
     col = alpha("gray52", 1.0))
points(-RNA_metabolic_processcds$logFC,
       -log10(pmax(10^-20, RNA_metabolic_processcds$adj.P.Val)),
       pch=16, cex=0.9,
       col="orange1"
)
points(-nuclear_pore_distributioncds$logFC,
       -log10(pmax(10^-20, nuclear_pore_distributioncds$adj.P.Val)),
       pch=16, cex=0.9,
       col="violet"
)
points(-protein_sumoylationcds$logFC,
       -log10(pmax(10^-20, protein_sumoylationcds$adj.P.Val)),
       pch=16, cex=0.9,
       col="palegreen3"
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
axis(side = 1, at = c(-6, -4, -2, 0, 2, 4, 6), labels=c("<1/64", "1/16", "1/4", "1", "4", "16", ">64"))
axis(side = 2, at = c(0, 5, 10, 15, 20), labels=c("0", "5", "10", "15", ">20"))
title(main="P(Gal) Guide Induction in TF-CDS Knockdown background")
dev.off()


#TF-UTR before vs after guide tet induction treatment_________________

## Set up values to plot
utr_insig <- filter(utr, utr$adj.P.Val > 0.05)
temp1 <- pmax(-6, utr$logFC)
temp2 <- pmin(6, temp1)
utr$LFC_plot <- -temp2
utr$sig_plot <- -log10(pmax(10^-20, utr$adj.P.Val))

pdf("utr_volcano.pdf", useDingbats=FALSE, width=5, height=5, colormodel="rgb")
## Plot with type="n" to skip drawing the points
plot(x = utr$LFC_plot, y = utr$sig_plot,
     type="n",
     xlim=c(-6, 6), ylim=c(0, 20),
     axes=FALSE,
     xlab="Log2 Fold Before vs After Guide tet Induction", 
     ylab = expression('Statistical Significance -log10(Q)'), 
     col="gray52")
points(-cellular_amino_acid_metabolic_processutr$logFC,
       -log10(pmax(10^-20, cellular_amino_acid_metabolic_processutr$adj.P.Val)),
       type="n",
       col="violet"
)
points(-tRNA_transcriptionutr$logFC,
       -log10(pmax(10^-20, tRNA_transcriptionutr$adj.P.Val)),
       type="n",
       col="palegreen3"
)
points(-tRNA_metabolic_processutr$logFC,
       -log10(pmax(10^-20, tRNA_metabolic_processutr$adj.P.Val)),
       type="n",
       col="orange1"
)
points(-tRNA_aminoacylationutr$logFC,
       -log10(pmax(10^-20, tRNA_aminoacylationutr$adj.P.Val)),
       type="n",
       col="orangered1"
)
points(-regulation_of_translationutr$logFC,
       -log10(pmax(10^-20, regulation_of_translationutr$adj.P.Val)),
       type="n",
       col="royalblue1"
)
points(-RNA_metabolic_processutr$logFC,
       -log10(pmax(10^-20, RNA_metabolic_processutr$adj.P.Val)),
       type="n",
       col="goldenrod1"
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
pointsfile <- "utr_volcano.png"
png(pointsfile, width=width, height=height, units="in", res=1500, bg="transparent")
par(mar=c(0,0,0,0))
## Exactly the same plot as above, but no type="n" and instead specify points, color, etc.
plot(x = utr$LFC_plot, y = utr$sig_plot,
     pch=16, cex=0.9,
     xlim=c(-6, 6), ylim=c(0, 20),
     axes=FALSE,
     xlab="Log2 Fold Before vs After Guide tet Induction", 
     ylab = expression('Statistical Significance -log10(Q)'), 
     col = alpha("gray52", 1.0))
points(-RNA_metabolic_processutr$logFC,
       -log10(pmax(10^-20, RNA_metabolic_processutr$adj.P.Val)),
       pch=16, cex=0.9,
       col="goldenrod1"
)
points(-cellular_amino_acid_metabolic_processutr$logFC,
       -log10(pmax(10^-20, cellular_amino_acid_metabolic_processutr$adj.P.Val)),
       pch=16, cex=0.9,
       col = alpha("mediumorchid2", 1.0)
)
points(-tRNA_transcriptionutr$logFC,
       -log10(pmax(10^-20, tRNA_transcriptionutr$adj.P.Val)),
       pch=16, cex=0.9,
       col = alpha("palegreen4", 1.0)
)
points(-tRNA_metabolic_processutr$logFC,
       -log10(pmax(10^-20, tRNA_metabolic_processutr$adj.P.Val)),
       pch=16, cex=0.9,
       col = alpha("orange1", 1.0)
)
points(-tRNA_aminoacylationutr$logFC,
       -log10(pmax(10^-20, tRNA_aminoacylationutr$adj.P.Val)),
       pch=16, cex=0.9,
       col = alpha("orangered1", 1.0)
)
points(-regulation_of_translationutr$logFC,
       -log10(pmax(10^-20, regulation_of_translationutr$adj.P.Val)),
       pch=16, cex=0.9,
       col="royalblue1"
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
axis(side = 1, at = c(-6, -4, -2, 0, 2, 4, 6), labels=c("<1/64", "1/16", "1/4", "1", "4", "16", ">64"))
axis(side = 2, at = c(0, 5, 10, 15, 20), labels=c("0", "5", "10", "15", ">20"))
title(main="P(Gal) TF-UTR Guide Induction")
dev.off()
