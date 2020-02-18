#Quantification of gcn4 western

hexo <- c(7758.841, 12698.154, 12782.154, 8144.447, 
          12080.962, 17121.104, 12142.497, 9738.619)
gcn4_tf <- c(8611.326, 9468.225, 13406.740, 10013.104, 6189.447, 9307.033, 6861.154, 4203.205)
samp <- c("Nup133_un", "Nup145_un", "UBC9_un", "Ulp1_un", 
          "Nup133_tet", "Nup145_tet", "UBC9_tet", "Ulp1_tet")
western <- data.frame(hexo, gcn4_tf, samp)
western$norm <- western$gcn4_tf/western$hexo

pdf("GCN4_western.pdf", useDingbats=FALSE, width=5, height=6)
par(mar=c(8,4,4,4)+.1)
barplot(western$norm, main="Western_gcn4_TF", las = 2)
dev.off()

