#Barplots of qPCRs

epistasis <- c(0.9953896791, 2.417194113,	0.5638304635,	0.6431973347,      
               8.948628997, 4.573609948, 1.520978753, 
               8.774599838, 9.350217988, 0.9244496602)
names(epistasis) <- c("WT", "ILV2", "ILV2_GCN2KO", "ILV2_GCN4KO", 
                      "HTS1", "HTS1_GCN2KO", "HTS1_GCN4KO", 
                      "RPC31", "RPC31_GCN2KO", "RPC31_GCN4KO")

pdf("Epistasis_qPCR.pdf", useDingbats=FALSE, width=5, height=6)
par(mar=c(8,4,4,4)+.1)
barplot(epistasis, main="Epistasis qPCR", las = 2)
dev.off()

#_________________________________________________________________________

trna_over <- c(12.85418884,	8.850956263,	7.000703958,	6.912303994,	7.520524765,	7.660826246)
names(trna_over) <- c("HTS1",	"HTS1_tRNA",	"RPC31",	"RPC31_tRNA",	"EFB1",	"EFB1_tRNA")

pdf("trna_overqPCR.pdf", useDingbats=FALSE, width=5, height=6)
par(mar=c(8,4,4,4)+.1)
barplot(trna_over, main="tRNA Overexpression qPCR", las = 2)
dev.off()

#_________________________________________________________________________
#second trial tRNA overexpression, including ILV2
trna_over <- c(2.170205357, 25.57374826, 2.277448332, 18.38148954, 1.964243858, 16.14578671, 
               2.020696631, 23.58326022, 0.9999549397, 13.84015737, 1.723910968, 20.09962471)
names(trna_over) <- c("HTS1-0-JR", "HTS1-0-met", "RPC31-0-JR", 
                      "RPC31-0-met", "ILV2-0-JR", "ILV2-0-met", 	
                      "HTS1-tet-JR", "HTS1-tet-met", "RPC31-tet-JR", 
                      "RPC31-tet-met", "ILV2-tet-JR", "ILV2-tet-met")
pdf("trna_overqPCR.pdf", useDingbats=FALSE, width=5, height=6)
par(mar=c(8,4,4,4)+.1)
barplot(trna_over, main="tRNA Overexpression qPCR", las = 2)
dev.off()

#_________________________________________________________________________
#second trial tRNA overexpression, including ILV2, log scaled
trna_over <- c(1.964243858, 16.14578671, 2.170205357, 25.57374826, 2.277448332, 18.38148954,  
               1.723910968, 20.09962471, 2.020696631, 23.58326022, 0.9999549397, 13.84015737)
names(trna_over) <- c("ILV2-0-JR", "ILV2-0-met", "HTS1-0-JR", "HTS1-0-met", 
                      "RPC31-0-JR", "RPC31-0-met", "ILV2-tet-JR", "ILV2-tet-met", 	
                      "HTS1-tet-JR", "HTS1-tet-met", "RPC31-tet-JR", 
                      "RPC31-tet-met")
pdf("trna_overqPCR.pdf", useDingbats=FALSE, width=5, height=6)
par(mar=c(8,4,4,4)+.1)
barplot(log10(trna_over + 1), main="tRNA Overexpression qPCR", las = 2, axes = FALSE)
axis(side = 2, at = c(log10(1), log10(2), log10(3), log10(4), 
                      log10(5), log10(6), log10(7), log10(8), 
                      log10(9), log10(10), log10(11), log10(21), log10(31)), 
     labels=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9","10", 
              "20", "30"))
dev.off()

#_________________________________________________________________________

rough_densit <- c(1.2930653, 3.633240856, 0.7815894579,  3.312170639, 1.416073381,
                  0.9553005699,  0.6009188712, 1.1)
names(rough_densit) <- c("null", "HTS1",	"RPC31", "ILV2",
                         "null_GCN2KO", "HTS1_GCN2KO",	"RPC31_GCN2KO", "ILV2_GCN2KO")

pdf("rough_densit2.pdf", useDingbats=FALSE, width=5, height=6)
par(mar=c(8,4,4,4)+.1)
barplot(rough_densit, main="Rough Western Densitometry 2-alpha", las = 2)
dev.off()

#_________________________________________________________________________
# HTS1 and RPC31 KD

KD <- c(2^2.391132797, 2^1.301237746, 2^2.846199948, 2^-2.794145511, 2^-2.309661692)
names(KD) <- c("HTS1-cit", "RPC31-cit",	"3-at-cit", "HTS1KD",
                         "RPC31KD")

pdf("KD.pdf", useDingbats=FALSE, width=5, height=6)
par(mar=c(8,4,4,4)+.1)
barplot(KD, main="HTS1 RPC31 KD", las = 2)
dev.off()

#_________________________________________________________________________
# HTS1 and RPC31 GO terms barplot

HTS1_sub <- HTS1_sum_mpralm[c(1,2,6,9)] 
RPC31_sub <- RPC31_sum_mpralm[c(1,2,6,9)] 
his4_sub <- his4_pooled_sum_mpralm[c(1,2,6,9)]
names(HTS1_sub) <- c("Guide", "HTS1_logFC", "HTS1_adj.P.val", "HTS1_Yorf1")
names(RPC31_sub) <- c("Guide", "RPC31_logFC", "RPC31_adj.P.val", "RPC31_Yorf1")
names(his4_sub) <- c("Guide", "his4_logFC", "his4_adj.P.val", "his4_Yorf1")

head(HTS1_sub)
head(RPC31_sub)
head(his4_sub)

HTS1vRPC31 <- merge(HTS1_sub, RPC31_sub, by = "Guide")
head(HTS1vRPC31)

HTS1vRPC31vhis4 <- merge(HTS1vRPC31, his4_sub, by = "Guide")
head(HTS1vRPC31vhis4)

all_sig <- filter(HTS1vRPC31vhis4, HTS1vRPC31vhis4$RPC31_adj.P.val < 0.05 & HTS1vRPC31vhis4$HTS1_adj.P.val < 0.05 & HTS1vRPC31vhis4$his4_adj.P.val < 0.05)
HTS1RPC31nohis4 <- filter(HTS1vRPC31vhis4, HTS1vRPC31vhis4$RPC31_adj.P.val < 0.05 & HTS1vRPC31vhis4$HTS1_adj.P.val < 0.05)
HTS1RPC31nohis4 <- filter(HTS1RPC31nohis4, !HTS1RPC31nohis4$Guide %in% all_sig$Guide)

a <- HTS1RPC31nohis4[complete.cases(HTS1RPC31nohis4[,4]),]
b <- HTS1vRPC31vhis4[complete.cases(HTS1vRPC31vhis4[,4]),]

write.table(a$HTS1_Yorf1, "~/Downloads/HTS1RPC31nohis4_sigdwn.txt", row.names = FALSE, col.names = FALSE, quote = FALSE) 
write.table(b$HTS1_Yorf1, "~/Downloads/HTS1RPC31nohis4_backgrnd.txt", row.names = FALSE, col.names = FALSE, quote = FALSE) 


HTS1RPC31nohis4 <- c(1.2930653, 3.633240856, 0.7815894579,  3.312170639, 1.416073381,
                  0.9553005699,  0.6009188712, 1.1)
names(HTS1RPC31nohis4) <- c("null", "HTS1",	"RPC31", "ILV2",
                         "null_GCN2KO", "HTS1_GCN2KO",	"RPC31_GCN2KO", "ILV2_GCN2KO")

pdf("HTS1RPC31nohis4.pdf", useDingbats=FALSE, width=5, height=6)
par(mar=c(8,4,4,4)+.1)
barplot(HTS1RPC31nohis4, main="GO terms HTS1 RPC31 no His4", las = 2)
dev.off()






GO_his4_3AT_dwn <- c(8.01e-07, 1.00e-05, 2.72e-04, 4.85e-02)
names(GO_his4_3AT_dwn) <- c("tRNA transcription by RNA polymerase III", 
                            "tRNA aminoacylation for protein translation", 
                            "cellular amino acid metabolic process", 
                            "tRNA processing")
logplot <- -log10(GO_his4_3AT_dwn)

pdf("His4_3AT_dwn.pdf", useDingbats=FALSE, width=5, height=7)
par(mar=c(17,4,4,4)+.1)
barplot(logplot, main="P(His4) 3AT Treatment Transcription Repressors", 
        ylab = "-log10(FDR)", 
        col =c("palegreen3", "orangered1", "violet", "orange1"), las = 2)
dev.off()


GO_his4_pooled_up <- c(5.95e-06, 2.38e-08, 3.43e-07, 9.24e-05)
names(GO_his4_pooled_up) <- c("tRNA transcription by RNA polymerase III", 
                              "tRNA aminoacylation for protein translation", 
                              "cellular amino acid metabolic process", 
                              "tRNA processing")
logplot <- -log10(GO_his4_pooled_up)

pdf("His4_pooled_up.pdf", useDingbats=FALSE, width=5, height=7)
par(mar=c(17,4,4,4)+.1)
barplot(logplot, main="P(His4) Transcription Activators", 
        ylab = "-log10(FDR)", 
        col =c("palegreen3", "orangered1", "violet", "orange1"), las = 2)
dev.off()





GO_his4_pooled_dwn <- c(8.01e-07, 1.00e-05, 2.72e-04, 4.85e-02)
names(GO_his4_pooled_dwn) <- c("tRNA transcription by RNA polymerase III", 
                               "tRNA aminoacylation for protein translation", 
                               "cellular amino acid metabolic process", 
                               "tRNA processing")
logplot <- -log10(GO_his4_pooled_dwn)

pdf("His4_pooled_dwn.pdf", useDingbats=FALSE, width=5, height=7)
par(mar=c(17,4,4,4)+.1)
barplot(logplot, main="P(His4) pooled Treatment Transcription Repressors", 
        ylab = "-log10(FDR)", 
        col =c("palegreen3", "orangered1", "violet", "orange1"), las = 2)
dev.off()

#GO terms that are in activators or repressors (could also do pie charts or even overlapping density plots
#that take pAdj vs count)

x <- cellular_amino_acid_metabolic_processhis4_3AT
sig <- filter(x, x$adj.P.Val < 0.05)
d <- density(log10(sig$adj.P.Val)) # returns the density data
plot(d)
lines(density(MyData$Column2))

x <- his4_pooledsig_up
subset(x, !(x$Guide %in% a$Guide))


GO_his4_pooled_dwn <- c(8.01e-07, 1.00e-05, 2.72e-04, 4.85e-02)
names(GO_his4_pooled_dwn) <- c("tRNA transcription by RNA polymerase III", 
                               "tRNA aminoacylation for protein translation", 
                               "cellular amino acid metabolic process", 
                               "tRNA processing", "other")
logplot <- -log10(GO_his4_pooled_dwn)

pdf("His4_pooled_dwn_numbers.pdf", useDingbats=FALSE, width=5, height=7)
par(mar=c(17,4,4,4)+.1)
barplot(logplot, main="P(His4) pooled Treatment Transcription Repressors", 
        ylab = "-log10(FDR)", 
        col =c("palegreen3", "orangered1", "violet", "orange1"), las = 2)
dev.off()

GO_his4_pooled_up <- c(8.01e-07, 1.00e-05, 2.72e-04, 4.85e-02)
names(GO_his4_pooled_up) <- c("tRNA transcription by RNA polymerase III", 
                              "tRNA aminoacylation for protein translation", 
                              "cellular amino acid metabolic process", 
                              "tRNA processing")
logplot <- -log10(GO_his4_pooled_up)

pdf("His4_pooled_up_numbers.pdf", useDingbats=FALSE, width=5, height=7)
par(mar=c(17,4,4,4)+.1)
barplot(logplot, main="P(His4) pooled Treatment Transcription Repressors", 
        ylab = "-log10(FDR)", 
        col =c("palegreen3", "orangered1", "violet", "orange1"), las = 2)
dev.off()

#qPCR: Efficient KD, Epistasis, tRNA overexpression, 

GO_his4_pooled_up <- c(0, 0, 0, 0)
names(GO_his4_pooled_up) <- c("tRNA transcription by RNA polymerase III", 
                              "tRNA aminoacylation for protein translation", 
                              "cellular amino acid metabolic process", 
                              "tRNA processing")
logplot <- -log10(GO_his4_pooled_up)

pdf("His4_Activatotion_by_Guide_KD.pdf", useDingbats=FALSE, width=5, height=7)
par(mar=c(5,4,4,4)+.1)
barplot(logplot, main="P(His4) pooled Treatment Transcription Repressors", 
        ylab = "-log10(FDR)", 
        col =c("palegreen3", "orangered1", "violet", "orange1"), las = 2)
dev.off()



