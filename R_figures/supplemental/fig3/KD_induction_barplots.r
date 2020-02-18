# His4 reporter induction

KD <- c(2^0.00000009311479809,	2^2.39113289,	2^-0.2086242928,	
        2^1.092613453,	2^-1.71901124,	2^1.127188708)
names(KD) <- c("HTS1-null", "HTS1-KD", "RPC31-null", "RPC31-KD", "null", "3-at")

pdf("induction.pdf", useDingbats=FALSE, width=5, height=6)
par(mar=c(8,4,4,4)+.1)
barplot(KD, main="His4 reporter induction", las = 2)
dev.off()


# HTS1 and RPC31 KD

KD <- c(2^6.892104455, 2^4.097958944, 2^2.300242089, 2^-0.009419603206)
names(KD) <- c("HTS1un", "HTS1KD", "RPC31un", "RPC31KD")

pdf("KD.pdf", useDingbats=FALSE, width=5, height=6)
par(mar=c(8,4,4,4)+.1)
barplot(log2(KD + 1), main="HTS1 RPC31 KD", las = 2, axes = FALSE)
axis(side = 2, at = c(log2(1), log2(2), log2(3), log2(5), log2(9), log2(17), log2(33), log2(65)), 
     labels=c("0", "1", "2", "4", "8", "16", "32", "64"))
dev.off()

