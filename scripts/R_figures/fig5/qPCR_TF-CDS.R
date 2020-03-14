#Barplots of qPCRs

cds_tf <- c(-0.000000001592368015, 2.300711998, 2.761421187, 2.547691513, 
            0.4620509363, 5.031665399, 4.243890226, 3.64003348)
names(cds_tf) <- c("NUP133_un", "NUP145_un", "UBC9_un", "ULP1_un", 
                      "NUP133_tet", "NUP145_tet", "UBC9_tet", "ULP1_tet")
pdf("cds_tf_qPCR.pdf", useDingbats=FALSE, width=5, height=6)
par(mar=c(8,4,4,4)+.1)
barplot(cds_tf, main="cds_tf qPCR", las = 2)
dev.off()

#re-normlaization
cds_tf <- c(0.9999999984, 3.300711998, 3.761421187, 3.547691513, 
            1.462050936, 6.031665399, 5.243890226, 4.64003348)
names(cds_tf) <- c("NUP133_un", "NUP145_un", "UBC9_un", "ULP1_un", 
                   "NUP133_tet", "NUP145_tet", "UBC9_tet", "ULP1_tet")
barplot(cds_tf, main="cds_tf qPCR", las = 2, axes = FALSE)
axis(side = 2, at = c(0, 1, 2, 3, 4, 5, 6), 
     labels=c("1/2", "0", "2", "4", "8", "16", "32"))

