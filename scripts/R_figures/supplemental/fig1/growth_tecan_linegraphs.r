#Line graphs of Tecan data
if (!requireNamespace("scales", quietly = TRUE))
  install.packages("scales")
library(scales)

est_growth <- read.csv("~/CiBER_seq_package/scripts/R_figures/supplemental/fig1/citrine_est_titration.csv",
                       stringsAsFactors=FALSE, header = TRUE, sep = ",")

est_growth$minint15 <- c(1:96)
est_growth$hour <- est_growth$minint15/4
tail(est_growth)

pdf("tecan_est_titration.pdf", useDingbats=FALSE, width=5, height=5, colormodel="rgb")
plot(est_growth$hour, est_growth$citperOD_128,
     xlim = c(0,50/4),
     xlab="Time (hr)",
     ylab="Citrine per OD value",
     type="l", 
     col="#859900", 
     lwd=2)
lines(est_growth$hour, est_growth$citperOD_0,
      type="l", 
      col="#b58900", 
      lwd=2)
lines(est_growth$hour, est_growth$citperOD_2,
      type="l", 
      col="#cb4b16", 
      lwd=2)
lines(est_growth$hour, est_growth$citperOD_4,
      type="l", 
      col="#dc322f", 
      lwd=2)
lines(est_growth$hour, est_growth$citperOD_8,
      type="l", 
      col="#d33682", 
      lwd=2)
lines(est_growth$hour, est_growth$citperOD_16,
      type="l", 
      col="#6c71c4", 
      lwd=2)
lines(est_growth$hour, est_growth$citperOD_32,
      type="l", 
      col="#268bd2", 
      lwd=2)
lines(est_growth$hour, est_growth$citperOD_64,
      type="l", 
      col="#2aa198", 
      lwd=2)
segments(est_growth$hour, est_growth$citperOD_0 - est_growth$STD_0, 
         est_growth$hour, est_growth$citperOD_0 + est_growth$STD_0, col = alpha("#b58900", 0.4), lwd = 2)
segments(est_growth$hour, est_growth$citperOD_2 - est_growth$STD_2, 
         est_growth$hour, est_growth$citperOD_2 + est_growth$STD_2, col = alpha("#cb4b16", 0.4), lwd = 2)
segments(est_growth$hour, est_growth$citperOD_4 - est_growth$STD_4, 
         est_growth$hour, est_growth$citperOD_4 + est_growth$STD_4, col = alpha("#dc322f", 0.4), lwd = 2)
segments(est_growth$hour, est_growth$citperOD_8 - est_growth$STD_8, 
         est_growth$hour, est_growth$citperOD_8 + est_growth$STD_8, col = alpha("#d33682", 0.4), lwd = 2)
segments(est_growth$hour, est_growth$citperOD_16 - est_growth$STD_16, 
         est_growth$hour, est_growth$citperOD_16 + est_growth$STD_16, col = alpha("#6c71c4", 0.4), lwd = 2)
segments(est_growth$hour, est_growth$citperOD_32 - est_growth$STD_32, 
         est_growth$hour, est_growth$citperOD_32 + est_growth$STD_32, col = alpha("#268bd2", 0.4), lwd = 2)
segments(est_growth$hour, est_growth$citperOD_64 - est_growth$STD_64, 
         est_growth$hour, est_growth$citperOD_64 + est_growth$STD_64, col = alpha("#2aa198", 0.4), lwd = 2)
segments(est_growth$hour, est_growth$citperOD_128 - est_growth$STD_128, 
         est_growth$hour, est_growth$citperOD_128 + est_growth$STD_128, col = alpha("#859900", 0.4), lwd = 2)
## Draw axes, title, etc.
#axis(side = 1, at = c(-4, -2, 0, 2, 4), labels=c("<1/16", "1/4", "1", "4", ">16"))
#axis(side = 2, at = c(0, 12, 24, 36, 48), labels=c("0", "3", "6", "9", "12"))
title(main = "Estradiol Titration Time course")
dev.off()

#gRNA line plots_________________________________________________________________
est_growth <- read.csv("~/CiBER_seq_package/scripts/R_figures/supplemental/fig1/tet_growth_pgal.csv",
                       stringsAsFactors=FALSE, header = TRUE, sep = ",")

tet_growth_pgal <- est_growth
tet_growth_pgal$hour <- tet_growth_pgal$Cycle.Nr./4
tail(tet_growth_pgal)

pdf("tecan_tet_pgal_gRNAtitration.pdf", useDingbats=FALSE, width=5, height=5, colormodel="rgb")
plot(tet_growth_pgal$hour, tet_growth_pgal$CitperOD_0,
     xlim = c(0,50/4),
     xlab="Time (hr)",
     ylab="Citrine per OD value",
     type="l", 
     col="#859900", 
     lwd=2)
lines(tet_growth_pgal$hour, tet_growth_pgal$CitperOD_3.9,
      type="l", 
      col="#b58900", 
      lwd=2)
lines(tet_growth_pgal$hour, tet_growth_pgal$CitperOD_7.8,
      type="l", 
      col="#cb4b16", 
      lwd=2)
lines(tet_growth_pgal$hour, tet_growth_pgal$CitperOD_15.6,
      type="l", 
      col="#dc322f", 
      lwd=2)
lines(tet_growth_pgal$hour, tet_growth_pgal$CitperOD_31.3,
      type="l", 
      col="#d33682", 
      lwd=2)
lines(tet_growth_pgal$hour, tet_growth_pgal$CitperOD_62.5,
      type="l", 
      col="#6c71c4", 
      lwd=2)
lines(tet_growth_pgal$hour, tet_growth_pgal$CitperOD_125,
      type="l", 
      col="#268bd2", 
      lwd=2)
lines(tet_growth_pgal$hour, tet_growth_pgal$CitperOD_250,
      type="l", 
      col="#2aa198", 
      lwd=2)
segments(tet_growth_pgal$hour, tet_growth_pgal$CitperOD_0 - tet_growth_pgal$STD_0, 
         tet_growth_pgal$hour, tet_growth_pgal$CitperOD_0 + tet_growth_pgal$STD_0, col = alpha("#b58900", 0.4), lwd = 2)
segments(tet_growth_pgal$hour, tet_growth_pgal$CitperOD_3.9 - tet_growth_pgal$STD_3.9, 
         tet_growth_pgal$hour, tet_growth_pgal$CitperOD_3.9 + tet_growth_pgal$STD_3.9, col = alpha("#cb4b16", 0.4), lwd = 2)
segments(tet_growth_pgal$hour, tet_growth_pgal$CitperOD_7.8 - tet_growth_pgal$STD_7.8, 
         tet_growth_pgal$hour, tet_growth_pgal$CitperOD_7.8 + tet_growth_pgal$STD_7.8, col = alpha("#dc322f", 0.4), lwd = 2)
segments(tet_growth_pgal$hour, tet_growth_pgal$CitperOD_15.6 - tet_growth_pgal$STD_15.6, 
         tet_growth_pgal$hour, tet_growth_pgal$CitperOD_15.6 + tet_growth_pgal$STD_15.6, col = alpha("#d33682", 0.4), lwd = 2)
segments(tet_growth_pgal$hour, tet_growth_pgal$CitperOD_31.3 - tet_growth_pgal$STD_31.3, 
         tet_growth_pgal$hour, tet_growth_pgal$CitperOD_31.3 + tet_growth_pgal$STD_31.3, col = alpha("#6c71c4", 0.4), lwd = 2)
segments(tet_growth_pgal$hour, tet_growth_pgal$CitperOD_62.5 - tet_growth_pgal$STD_62.5, 
         tet_growth_pgal$hour, tet_growth_pgal$CitperOD_62.5 + tet_growth_pgal$STD_62.5, col = alpha("#268bd2", 0.4), lwd = 2)
segments(tet_growth_pgal$hour, tet_growth_pgal$CitperOD_125 - tet_growth_pgal$STD_125, 
         tet_growth_pgal$hour, tet_growth_pgal$CitperOD_125 + tet_growth_pgal$STD_125, col = alpha("#2aa198", 0.4), lwd = 2)
segments(tet_growth_pgal$hour, tet_growth_pgal$CitperOD_250 - tet_growth_pgal$STD_250, 
         tet_growth_pgal$hour, tet_growth_pgal$CitperOD_250 + tet_growth_pgal$STD_250, col = alpha("#859900", 0.4), lwd = 2)
## Draw axes, title, etc.
#axis(side = 1, at = c(-4, -2, 0, 2, 4), labels=c("<1/16", "1/4", "1", "4", ">16"))
#axis(side = 2, at = c(0, 12, 24, 36, 48), labels=c("0", "3", "6", "9", "12"))
title(main = "pGal gRNA tet Titration Time course")
dev.off()

est_growth <- read.csv("~/CiBER_seq_package/scripts/R_figures/supplemental/fig1/tet_growth_adh1.csv",
                       stringsAsFactors=FALSE, header = TRUE, sep = ",")

tet_growth_padh1 <- est_growth
tet_growth_padh1$hour <- tet_growth_padh1$Cycle.Nr./4
tail(tet_growth_padh1)

pdf("tecan_tet_adh1_gRNAtitration.pdf", useDingbats=FALSE, width=5, height=5, colormodel="rgb")
plot(tet_growth_padh1$hour, tet_growth_padh1$CitperOD_0,
     xlim = c(0,50/4),
     xlab="Time (hr)",
     ylab="Citrine per OD value",
     type="l", 
     col="#859900", 
     lwd=2)
lines(tet_growth_padh1$hour, tet_growth_padh1$CitperOD_3.9,
      type="l", 
      col="#b58900", 
      lwd=2)
lines(tet_growth_padh1$hour, tet_growth_padh1$CitperOD_7.8,
      type="l", 
      col="#cb4b16", 
      lwd=2)
lines(tet_growth_padh1$hour, tet_growth_padh1$CitperOD_15.6,
      type="l", 
      col="#dc322f", 
      lwd=2)
lines(tet_growth_padh1$hour, tet_growth_padh1$CitperOD_31.3,
      type="l", 
      col="#d33682", 
      lwd=2)
lines(tet_growth_padh1$hour, tet_growth_padh1$CitperOD_62.5,
      type="l", 
      col="#6c71c4", 
      lwd=2)
lines(tet_growth_padh1$hour, tet_growth_padh1$CitperOD_125,
      type="l", 
      col="#268bd2", 
      lwd=2)
lines(tet_growth_padh1$hour, tet_growth_padh1$CitperOD_250,
      type="l", 
      col="#2aa198", 
      lwd=2)
segments(tet_growth_padh1$hour, tet_growth_padh1$CitperOD_0 - tet_growth_padh1$STD_0, 
         tet_growth_padh1$hour, tet_growth_padh1$CitperOD_0 + tet_growth_padh1$STD_0, col = alpha("#b58900", 0.4), lwd = 2)
segments(tet_growth_padh1$hour, tet_growth_padh1$CitperOD_3.9 - tet_growth_padh1$STD_3.9, 
         tet_growth_padh1$hour, tet_growth_padh1$CitperOD_3.9 + tet_growth_padh1$STD_3.9, col = alpha("#cb4b16", 0.4), lwd = 2)
segments(tet_growth_padh1$hour, tet_growth_padh1$CitperOD_7.8 - tet_growth_padh1$STD_7.8, 
         tet_growth_padh1$hour, tet_growth_padh1$CitperOD_7.8 + tet_growth_padh1$STD_7.8, col = alpha("#dc322f", 0.4), lwd = 2)
segments(tet_growth_padh1$hour, tet_growth_padh1$CitperOD_15.6 - tet_growth_padh1$STD_15.6, 
         tet_growth_padh1$hour, tet_growth_padh1$CitperOD_15.6 + tet_growth_padh1$STD_15.6, col = alpha("#d33682", 0.4), lwd = 2)
segments(tet_growth_padh1$hour, tet_growth_padh1$CitperOD_31.3 - tet_growth_padh1$STD_31.3, 
         tet_growth_padh1$hour, tet_growth_padh1$CitperOD_31.3 + tet_growth_padh1$STD_31.3, col = alpha("#6c71c4", 0.4), lwd = 2)
segments(tet_growth_padh1$hour, tet_growth_padh1$CitperOD_62.5 - tet_growth_padh1$STD_62.5, 
         tet_growth_padh1$hour, tet_growth_padh1$CitperOD_62.5 + tet_growth_padh1$STD_62.5, col = alpha("#268bd2", 0.4), lwd = 2)
segments(tet_growth_padh1$hour, tet_growth_padh1$CitperOD_125 - tet_growth_padh1$STD_125, 
         tet_growth_padh1$hour, tet_growth_padh1$CitperOD_125 + tet_growth_padh1$STD_125, col = alpha("#2aa198", 0.4), lwd = 2)
segments(tet_growth_padh1$hour, tet_growth_padh1$CitperOD_250 - tet_growth_padh1$STD_250, 
         tet_growth_padh1$hour, tet_growth_padh1$CitperOD_250 + tet_growth_padh1$STD_250, col = alpha("#859900", 0.4), lwd = 2)
## Draw axes, title, etc.
#axis(side = 1, at = c(-4, -2, 0, 2, 4), labels=c("<1/16", "1/4", "1", "4", ">16"))
#axis(side = 2, at = c(0, 12, 24, 36, 48), labels=c("0", "3", "6", "9", "12"))
title(main = "pAdh1 gRNA tet Titration Time course")
dev.off()