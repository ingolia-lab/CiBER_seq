#DNA replicates bc6__________________________________________
pdf("DNA_bc6_replicates.pdf", useDingbats=FALSE, width=5, height=5, colormodel="rgb")
par(pty="s")
plot(log10(bc6$V2+1), log10(bc6$V3+1), pch=20, 
     col = alpha("black", 0.1), axes = FALSE, 
     type="n")

## Capture the plot window coordinates
coords <- par("usr")
gx <- grconvertX(coords[1:2], "user", "inches")
gy <- grconvertY(coords[3:4], "user", "inches")
width <- max(gx) - min(gx)
height <- max(gy) - min(gy)

pointsfile <- "DNA_bc6.png"
png(pointsfile, width=width, height=height, units="in", res=1500, bg="transparent")
par(mar=c(0,0,0,0))
plot(log10(bc6$V2+1), log10(bc6$V3+1), pch=20, 
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
title(main = "DNA_replicates")
dev.off()