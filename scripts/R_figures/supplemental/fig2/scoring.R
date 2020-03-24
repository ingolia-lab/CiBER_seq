## mpralm scoring results for P(HIS4)
mpralmFile <- '/mnt/ingolialab/rmuller1/CiBER_seq_package/all_raw_fastq/HIS4_PGK1_pooled/his4_pooled_sum_mpralm.txt'
mpralmAll <- read.delim(mpralmFile, stringsAsFactors=FALSE, row.names=1)
mpralm <- mpralmAll[,c("Guide", "logFC", "AveExpr", "P.Value", "adj.P.Val")]
neg <- mpralm[grepl("No_gRNA", mpralm$Guide),]
mpralm <- mpralm[!grepl("No_gRNA", mpralm$Guide),]

## Guide information -- targeting and scoring
seqGoodTargetFile <- '/mnt/ingolialab/ingolia/Prog/yeast-crispri/library_v1/sequence-good-targets.txt'
seqRedundantTargetFile <- '/mnt/ingolialab/ingolia/Prog/yeast-crispri/library_v1/sequence-redundant-targets.txt'
yorfCleanFile <- '/mnt/ingolialab/ingolia/Prog/yeast-crispri/library_v1/yorf-clean-target.txt'
rescoreFile <- '/mnt/ingolialab/ingolia/Prog/yeast-crispri/mcglincy_2019/phenotype/figures/rescore.txt'

targGood <- read.delim(seqGoodTargetFile, stringsAsFactors=FALSE)
targRed <- read.delim(seqRedundantTargetFile, stringsAsFactors=FALSE)
yorfClean <- read.delim(yorfCleanFile, stringsAsFactors=FALSE, header=FALSE)
rescore <- read.delim(rescoreFile, stringsAsFactors=FALSE)

targetOligoFile <- 'target-oligos.txt'
targetOligos <- read.delim(targetOligoFile, stringsAsFactors=FALSE, row.names=NULL)

## Link mpralm results with targeting through the oligo sequence
mpralm$Oligo <- targetOligos[match(mpralm$Guide, targetOligos$Yorf_GNo), "Oligo"]
mpralm <- cbind.data.frame(mpralm,
                           targGood[match(mpralm$Oligo, targGood$Oligo),
                                    c("TargetLoc", "Yorf1", "Offset1", "StartType1", "Yorfs")])
mpralm$Oligo <- NULL
mpralm$rescore <- rescore[match(mpralm$TargetLoc, rescore$TargetLoc), "score"]

## Pick out candidates: genes with a significant (q < 0.01), strong (|logFC| > 1) guide
sig <- mpralm[mpralm$adj.P.Val<0.01 & abs(mpralm$logFC) > 1,]
cat(sprintf("%d significant guides\n", nrow(sig)))
## Further filter on "clean" (non-divergent) promoters
sigclean <- sig[sig$Yorf1 %in% yorfClean$V1,]
cleantargs <- unique(sigclean$Yorf1)
cat(sprintf("%d hits at clean promoters\n", length(cleantargs)))

cands <- mpralm[mpralm$Yorf1 %in% cleantargs,]
cands <- cands[order(cands$Yorf1),]
cands$scorebin <- round(cands$rescore)

signifScore <- ecdf(cands[cands$adj.P.Val < 0.01,]$rescore)
nonsigScore <- ecdf(cands[cands$adj.P.Val > 0.5,]$rescore)
allScore <- ecdf(cands$rescore)

sbinneg2 <- ecdf(abs(cands[cands$scorebin <= -2,]$logFC))
sbinneg1 <- ecdf(abs(cands[cands$scorebin == -1,]$logFC))
sbinzero <- ecdf(abs(cands[cands$scorebin ==  0,]$logFC))
sbinpos1 <- ecdf(abs(cands[cands$scorebin ==  1,]$logFC))
sbinpos2 <- ecdf(abs(cands[cands$scorebin >=  2,]$logFC))

library("RColorBrewer")

rdbu <- brewer.pal(9, "RdBu")
puor <- brewer.pal(9, "PuOr")

pdf("score-by-signif.pdf", width=6, height=4, useDingbats=FALSE)
xs <- seq(-5,5,0.05)
plot(xs, allScore(xs),
     type="l", lwd=3, col="grey",
     xlim=c(-5,5), ylim=c(0,1), yaxt="n",
     xlab="logit activity score", ylab="Fract guides",
     xaxp=c(-4,4,4))
lines(xs, signifScore(xs),
      col=rdbu[[2]], lwd=2)
lines(xs, nonsigScore(xs),
      col=rdbu[[8]], lwd=2)
legend(x="topleft", bty="n",
       legend=c("q < 0.01", "all", "q > 0.5"),
       lwd=c(2,3,2), col=c(rdbu[[2]], "grey", rdbu[[8]]),
       text.col=c(rdbu[[2]], "grey", rdbu[[8]]))
dev.off()

pdf("change-by-score.pdf", width=6, height=4, useDingbats=FALSE)
xs <- seq(0, 3, 0.05)
plot(xs, sbinneg2(xs),
     type="l", lwd=2, col=puor[[9]],
     xlim=c(0,3), ylim=c(0,1),
     xlab="Absolute fold-change", ylab="Fract guides",
     xaxt="n", yaxt="n")
lines(xs, sbinneg1(xs),
      col=puor[[7]], lwd=2)
lines(xs, sbinzero(xs),
      col="grey", lwd=2)
lines(xs, sbinpos1(xs),
      col=puor[[3]], lwd=2)
lines(xs, sbinpos2(xs),
      col=puor[[1]], lwd=2)      
axis(side=1,
     at=log2(seq(1,8)),
     labels=seq(1,8))
legend(x="bottomright", bty="n",
       legend=c("-2", "-1", "0", "1", "2"),
       lwd=2, col=c(puor[[9]], puor[[7]], "grey", puor[[3]], puor[[1]]),
       text.col=c(puor[[9]], puor[[7]], "grey", puor[[3]], puor[[1]]))
dev.off()


