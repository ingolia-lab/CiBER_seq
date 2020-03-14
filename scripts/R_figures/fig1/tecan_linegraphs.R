#Line graphs of Tecan data

est_growth <- see_CGA_estradiol_growth_stuff...citperOD
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
tet_growth_pgal <- `92118_Estra_Tet_46...p(gal)g2_final`
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


tet_growth_padh1 <- `92118_Estra_Tet_46...p(adh1)g2_final`
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





#extra experimenting
library(ggplot2)
e <- ggplot(est_growth, aes(hour))
e
h + geom_ribbon(aes(ymin=citperOD_64-STD_64, ymax=citperOD_64+STD_64))
h + geom_area(aes(y = level))
h + geom_ribbon(aes(ymin = citperOD_64 - 1, ymax = citperOD_64 + 1), fill = "grey70")

+
  geom_line(aes(y = level))


huron <- data.frame(year = 1875:1972, level = as.vector(LakeHuron))
h <- ggplot(huron, aes(year))

h + geom_ribbon(aes(ymin=0, ymax=level))
h +
  geom_ribbon(aes(ymin = level - 1, ymax = level + 1), fill = "grey70") +
  geom_line(aes(y = level))

ggplot(est_growth, aes(citperOD_0))
plot(est_growth$citperOD_0)




#Graph at max expression for the different estradiol concentrations
max_val <- c(272.2143115, 2386.113856, 4387.855756, 
             7603.326305, 9134.297129, 13779.35727, 
             15677.94506, 19800.59626	)
conc <- c(0, 2, 4, 8, 16, 32, 64, 128)
max_std <- c(	61.33724757, 79.48970905, 108.8775488, 
                 	298.5367273, 722.8422138, 458.6050947, 
                 	544.1148317, 1663.170019)
est_max <- data.frame(max_val, conc, max_std)

pdf("est_max.pdf", useDingbats=FALSE, width=5, height=5, colormodel="rgb")
plot(log10(est_max$conc + 1), log10(est_max$max_val),
     xlab="Concentration (nM)",
     ylab="Citrine per OD value",
     axes = FALSE, 
     type="l", 
     col="#b58900", 
     lwd=2)
segments(log10(est_max$conc + 1), log10(est_max$max_val - est_max$max_std), 
         log10(est_max$conc + 1), log10(est_max$max_val + est_max$max_std), col = alpha("#b58900", 0.8), lwd = 2)
## Draw axes, title, etc.
axis(side = 1, at = c(log10(1), log10(3), log10(5), log10(9), log10(17), log10(33), log10(65), log10(129)), 
     labels=c("0", "2", "4", "8", "16", "32", "64", "128"))
axis(side = 2, at = c(log10(1), log10(101), log10(201), log10(301), log10(401),
                      log10(501), log10(601), log10(701), log10(801), log10(901),
                      log10(1001), log10(2001), log10(3001), log10(4001), 
                      log10(5001), log10(6001), log10(7001), log10(8001), log10(9001), 
                      log10(10001), log10(11001), log10(12001), log10(13001), log10(14001), log10(15001), 
                      log10(16001), log10(17001), log10(18001), log10(19001), log10(20001)), 
     labels=c("0", "100", "200", "300", "400", 
              "500", "600", "700", "800", "900",
              "1000", "2000", "3000", "4000", 
              "5000", "6000", "7000", "8000", "9000", 
              "10000", "11000", "12000", "13000", "14000",
              "15000", "16000", "17000", "18000", "19000", "20000"))
title(main = "Citrine Expression vs Estradiol Concentration")
dev.off()


#Graph at max expression for the different tet concentrations
pgal_max_val <- c(10260.14284, 10717.84645, 9494.675883, 
                  6007.743148, 3589.354932, 3525.863596, 
                  3802.993432, 4051.52212)
pgal_conc <- c(0, 3.9, 7.8, 15.6, 31.3, 62.5, 125, 250)
pgal_max_std <- c(820.0955718, 340.1950584, 276.9320225, 
                  	335.4882575, 217.7146642, 28.12893389, 
                  	126.5493831, 411.3870831)
pgal_max <- data.frame(pgal_max_val, pgal_conc, pgal_max_std)

padh1_max_val <- c(15264.31562, 14510.31592, 13186.95844, 	
                   7192.598542, 4578.049562, 4026.019972, 
                   3385.416766, 2509.807193)
padh1_conc <- c(0, 3.9, 7.8, 15.6, 31.3, 62.5, 125, 250)
padh1_max_std <- c(	1398.837151, 679.9789641, 2424.579898, 
                   	2227.157657, 281.8380073, 328.4391969, 
                   	678.489613, 536.5066049)
padh1_max <- data.frame(padh1_max_val, padh1_conc, padh1_max_std)

pdf("tetgRNA_max.pdf", useDingbats=FALSE, width=5.5, height=5.5, colormodel="rgb")
plot(log10(pgal_max$pgal_conc + 1), pgal_max$pgal_max_val,
     xlab="Concentration (ng/mL)",
     ylab="Citrine per OD value",
     ylim = c(0,17000), 
     axes = FALSE, 
     type="l", 
     col="#6c71c4", 
     lwd=2)
lines(log10(padh1_max$padh1_conc + 1), padh1_max$padh1_max_val, 
      type="l", 
      col="#dc322f", 
      lwd=2)
segments(log10(pgal_max$pgal_conc + 1), pgal_max$pgal_max_val - pgal_max$pgal_max_std, 
         log10(pgal_max$pgal_conc + 1), pgal_max$pgal_max_val + pgal_max$pgal_max_std, col = alpha("#6c71c4", 0.8), lwd = 2)
segments(log10(padh1_max$padh1_conc + 1), padh1_max$padh1_max_val - padh1_max$padh1_max_std, 
         log10(padh1_max$padh1_conc + 1), padh1_max$padh1_max_val + padh1_max$padh1_max_std, col = alpha("#dc322f", 0.8), lwd = 2)
## Draw axes, title, etc.
axis(side = 1, at = c(log10(1), log10(4.9), log10(8.8), log10(16.6), log10(32.3), log10(63.5), log10(126), log10(251)), 
     labels=c("0", "3.9", "7.8", "15.6", "31.3", "62.5", "125", "250"))
axis(side = 2, at = c(0, 4000, 8000, 12000, 16000), labels=c("0", "4000", "8000", "12000", "16000"))
title(main = "Citrine Expression vs tet Concentration")
dev.off()
