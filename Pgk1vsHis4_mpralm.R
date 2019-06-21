#start with trimmed and tabulated datasets "all_his4_counts.txt", "all_his4_counts2.txt", 
# "all_pgk1_counts.txt", and "all_pgk1_counts2.txt"


#loading packages and libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("mpra", version = "3.8")

library(mpra)

#creating dataframes that mpralm can use
#all_pooled_pgk1_32 is filtering the pooled dataset with 32 read cut-off for both pre IVT samples
dna_pgk1 <- all_pooled_pgk1_32
dna_pgk1 <- merge(dna_pgk1, grna.assign.barcode.grna.good, by="barcode")
dna_pgk1 <- data.frame(row.names=dna_pgk1$barcode, dna_pgk1$IVT_preL, dna_pgk1$IVT_postL, 
                       dna_pgk1$IVT_preR, dna_pgk1$IVT_postR)
names(dna_pgk1) <- c("preL", "postL", "preR", "postR")
head(dna_pgk1)

rna_pgk1 <- all_pooled_pgk1_32
rna_pgk1 <- merge(rna_pgk1, grna.assign.barcode.grna.good, by="barcode")
rna_pgk1 <- data.frame(row.names=rna_pgk1$barcode, rna_pgk1$RNA_preL, rna_pgk1$RNA_postL, 
                       rna_pgk1$RNA_preR, rna_pgk1$RNA_postR)
names(rna_pgk1) <- c("preL", "postL", "preR", "postR")
head(rna_pgk1)

eid_pgk1 <- all_pooled_pgk1_32
eid_pgk1 <- merge(eid_pgk1, grna.assign.barcode.grna.good, by="barcode")
eid_pgk1 <- as.character(eid_pgk1$guide)
head(eid_pgk1)

#mpra dataset construction
mpraset_PGK1 <- MPRASet(DNA = dna_pgk1, RNA = rna_pgk1, eid = eid_pgk1, eseq = NULL, barcode = NULL)
mpraset_PGK1

#mpra data set analysis
data(mpraset_PGK1)
design <- data.frame(intcpt = 1,
                     pre = grepl("pre", colnames(mpraset_PGK1)))
mpralm_fit <- mpralm(object = mpraset_PGK1, design = design,
                     aggregate = "mean", normalize = TRUE,
                     model_type = "indep_groups", plot = TRUE)
toptab <- topTable(mpralm_fit, coef = 2, number = Inf)
head(toptab)

#modifications to output analysis
pgk1_mpralm <- toptab 
install.packages("data.table")
library(data.table)
setDT(pgk1_mpralm, keep.rownames = TRUE)[]
names(pgk1_mpralm)[names(pgk1_mpralm) == "barcode"] <- "guide"
pgk1_mpralm$yorf <- pgk1_mpralm$guide
pgk1_mpralm$yorf <- sub("_[0-9]*$", "", pgk1_mpralm$yorf)
pgk1_mpralm$gene <- sgd[match(pgk1_mpralm$yorf, sgd$name), "gene"]
pgk1_mpralm$desc <- sgd[match(pgk1_mpralm$yorf, sgd$name), "desc"]
head(pgk1_mpralm)

#______________repeat analysis with His4

#creating dataframes that mpralm can use
dna_his4 <- all_pooled_his4_32
dna_his4 <- merge(dna_his4, grna.assign.barcode.grna.good, by="barcode")
dna_his4 <- data.frame(row.names=dna_his4$barcode, dna_his4$IVT_preL, dna_his4$IVT_postL, 
                       dna_his4$IVT_preR, dna_his4$IVT_postR)
names(dna_his4) <- c("preL", "postL", "preR", "postR")
head(dna_his4)

rna_his4 <- all_pooled_his4_32
rna_his4 <- merge(rna_his4, grna.assign.barcode.grna.good, by="barcode")
rna_his4 <- data.frame(row.names=rna_his4$barcode, rna_his4$RNA_preL, rna_his4$RNA_postL, 
                       rna_his4$RNA_preR, rna_his4$RNA_postR)
names(rna_his4) <- c("preL", "postL", "preR", "postR")
head(rna_his4)

eid_his4 <- all_pooled_his4_32
eid_his4 <- merge(eid_his4, grna.assign.barcode.grna.good, by="barcode")
eid_his4 <- as.character(eid_his4$guide)
head(eid_his4)

#mpra dataset construction
mpraset_his4 <- MPRASet(DNA = dna_his4, RNA = rna_his4, eid = eid_his4, eseq = NULL, barcode = NULL)
mpraset_his4

#mpra data set analysis
data(mpraset_his4)
design <- data.frame(intcpt = 1,
                     pre = grepl("pre", colnames(mpraset_his4)))
mpralm_fit <- mpralm(object = mpraset_his4, design = design,
                     aggregate = "mean", normalize = TRUE,
                     model_type = "indep_groups", plot = TRUE)
toptab <- topTable(mpralm_fit, coef = 2, number = Inf)
head(toptab)

#modifications to output analysis
his4_mpralm <- toptab 
install.packages("data.table")
library(data.table)
setDT(his4_mpralm, keep.rownames = TRUE)[]
names(his4_mpralm)[names(his4_mpralm) == "rn"] <- "guide"
his4_mpralm$yorf <- his4_mpralm$guide
his4_mpralm$yorf <- sub("_[0-9]*$", "", his4_mpralm$yorf)
his4_mpralm$gene <- sgd[match(his4_mpralm$yorf, sgd$name), "gene"]
his4_mpralm$desc <- sgd[match(his4_mpralm$yorf, sgd$name), "desc"]
head(his4_mpralm)

#saving mpralm analysis
write.table(pgk1_mpralm, "~/Documents/pgk1_mpralm.txt", sep="\t")
write.table(his4_mpralm, "~/Documents/his4_mpralm.txt", sep="\t")

