#Line graphs of Tecan data
if (!requireNamespace("scales", quietly = TRUE))
  install.packages("scales")
library(scales)

bc1 <- read.csv("~/CiBER_github/CiBER_seq/scripts/R_figures/supplemental/fig1/bcodes_bc1.csv",
                       stringsAsFactors=FALSE, header = TRUE, sep = ",")
bc2 <- read.csv("~/CiBER_github/CiBER_seq/scripts/R_figures/supplemental/fig1/bcodes_bc2.csv",
                stringsAsFactors=FALSE, header = TRUE, sep = ",")
bc3 <- read.csv("~/CiBER_github/CiBER_seq/scripts/R_figures/supplemental/fig1/bcodes_bc3.csv",
                stringsAsFactors=FALSE, header = TRUE, sep = ",")

bc1 <- filter(bc1, bc1$PE.Tx.Red.YG.A > 0)
bc2 <- filter(bc2, bc2$PE.Tx.Red.YG.A > 0)
bc3 <- filter(bc3, bc3$PE.Tx.Red.YG.A > 0)

pdf("barcodes_estradiol.pdf", useDingbats=FALSE, width=5, height=5, colormodel="rgb")
d <- density(bc1$PE.Tx.Red.YG.A)
plot(d, xlim = range(0,1000), col="red", main = "")
lines(density(bc2$PE.Tx.Red.YG.A), col="green")
lines(density(bc3$PE.Tx.Red.YG.A), col="blue")
title(main = "Estradiol Titration FACS")
abline(v = mean(bc1$PE.Tx.Red.YG.A), col="red")
abline(v = mean(bc2$PE.Tx.Red.YG.A), col="green")
abline(v = mean(bc3$PE.Tx.Red.YG.A), col="blue")
dev.off()

median(bc1$PE.Tx.Red.YG.A)
median(bc2$PE.Tx.Red.YG.A)
median(bc3$PE.Tx.Red.YG.A)


