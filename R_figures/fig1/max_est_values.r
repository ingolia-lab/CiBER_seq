#Graph at max expression for the different estradiol concentrations
max_val <- c(272.2143115, 2386.113856, 4387.855756, 
             7603.326305, 9134.297129, 13779.35727, 
             15677.94506, 19800.59626	)
conc <- c(0, 2, 4, 8, 16, 32, 64, 128)
max_std <- c(61.33724757, 79.48970905, 108.8775488, 
                 298.5367273, 722.8422138, 458.6050947, 
                 544.1148317, 1663.170019)
est_max <- data.frame(max_val, conc, max_std)

pdf("est_max.pdf", useDingbats=FALSE, width=5, height=5, colormodel="rgb")
plot(log10(est_max$conc + 1), est_max$max_val,
     xlab="Concentration (nM)",
     ylab="Citrine per OD value",
     axes = FALSE, 
     type="l", 
     col="#b58900", 
     lwd=2)
segments(log10(est_max$conc + 1), est_max$max_val - est_max$max_std, 
         log10(est_max$conc + 1), est_max$max_val + est_max$max_std, col = alpha("#b58900", 0.8), lwd = 2)
## Draw axes, title, etc.
axis(side = 1, at = c(log10(1), log10(3), log10(5), log10(9), log10(17), log10(33), log10(65), log10(129)), 
     labels=c("0", "2", "4", "8", "16", "32", "64", "128"))
axis(side = 2, at = c(0, 5000, 10000, 15000, 20000), labels=c("0", "5000", "10000", "15000", "20000"))
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
padh1_max_std <- c(1398.837151, 679.9789641, 2424.579898, 
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
