
#Load in mpralm datasets
options(stringsAsFactors=FALSE)
if (!exists("his4_3AT")) {
  his4_3AT <- read.delim("~/test_sandbox/gz_to_mpralm/rehearsal_dinner/all_raw_fasta_gz/HIS4_PGK1_3AT/his4_postv3AT_sum_mpralm.txt",
                                    stringsAsFactors=FALSE)
}

options(stringsAsFactors=FALSE)
if (!exists("pgk1_3AT")) {
  pgk1_3AT <- read.delim("~/test_sandbox/gz_to_mpralm/rehearsal_dinner/all_raw_fasta_gz/HIS4_PGK1_3AT/pgk1_postv3AT_sum_mpralm.txt",
                                    stringsAsFactors=FALSE)
}

options(stringsAsFactors=FALSE)
if (!exists("his4_pooled")) {
  his4_pooled <- read.delim("~/test_sandbox/gz_to_mpralm/rehearsal_dinner/all_raw_fasta_gz/HIS4_PGK1_pooled/his4_pooled_sum_mpralm.txt",
                                       stringsAsFactors=FALSE)
}

options(stringsAsFactors=FALSE)
if (!exists("pgk1_pooled")) {
  pgk1_pooled <- read.delim("~/test_sandbox/gz_to_mpralm/rehearsal_dinner/all_raw_fasta_gz/HIS4_PGK1_pooled/pgk1_pooled_sum_mpralm.txt",
                                       stringsAsFactors=FALSE)
}

options(stringsAsFactors=FALSE)
if (!exists("HTS1")) {
  HTS1 <- read.delim("~/test_sandbox/gz_to_mpralm/rehearsal_dinner/all_raw_fasta_gz/HTS1_RPC31/HTS1_sum_mpralm.txt",
                                stringsAsFactors=FALSE)
}

options(stringsAsFactors=FALSE)
if (!exists("RPC31")) {
  RPC31 <- read.delim("~/test_sandbox/gz_to_mpralm/rehearsal_dinner/all_raw_fasta_gz/HTS1_RPC31/RPC31_sum_mpralm.txt",
                                 stringsAsFactors=FALSE)
}

options(stringsAsFactors=FALSE)
if (!exists("TF_CDS")) {
  cds <- read.delim("~/test_sandbox/gz_to_mpralm/rehearsal_dinner/all_raw_fasta_gz/GCN4_CDS_UTR/cds_sum_mpralm.txt",
                               stringsAsFactors=FALSE)
}

options(stringsAsFactors=FALSE)
if (!exists("TF_UTR")) {
  utr <- read.delim("~/test_sandbox/gz_to_mpralm/rehearsal_dinner/all_raw_fasta_gz/GCN4_CDS_UTR/utr_sum_mpralm.txt",
                               stringsAsFactors=FALSE)
}

#GO analysis, generating His4 3AT background and significant up and down text files
his4_3AT_noempty <- his4_3AT
his4_3AT_noempty[his4_3AT_noempty == "Neg_ctrl"] <- NA
his4_3AT_noempty <- his4_3AT_noempty[complete.cases(his4_3AT_noempty[,9]),]

his4_3AT_background <- his4_3AT_noempty
his4_3AT_sig_up <- his4_3AT_noempty
his4_3AT_sig_down <- his4_3AT_noempty
his4_3AT_sig_up <- filter(his4_3AT_sig_up, his4_3AT_sig_up$logFC < -0.5 & his4_3AT_sig_up$adj.P.Val < 0.05)
his4_3AT_sig_down <- filter(his4_3AT_sig_down, his4_3AT_sig_down$logFC > 0.5 & his4_3AT_sig_down$adj.P.Val < 0.05)

write(his4_3AT_background$Yorf1, "/mnt/ingolialab/rmuller1/test_sandbox/gz_to_mpralm/rehearsal_dinner/scripts/GO_analysis_files/his4_3AT_background.txt") 
write(his4_3AT_sig_up$Yorf1, "/mnt/ingolialab/rmuller1/test_sandbox/gz_to_mpralm/rehearsal_dinner/scripts/GO_analysis_files/his4_3AT_sig_up.txt")
write(his4_3AT_sig_down$Yorf1, "/mnt/ingolialab/rmuller1/test_sandbox/gz_to_mpralm/rehearsal_dinner/scripts/GO_analysis_files/his4_3AT_sig_down.txt")

#GO analysis, generating Pgk1 3AT background and significant up and down text files
pgk1_3AT_noempty <- pgk1_3AT
pgk1_3AT_noempty[pgk1_3AT_noempty == "Neg_ctrl"] <- NA
pgk1_3AT_noempty <- pgk1_3AT_noempty[complete.cases(pgk1_3AT_noempty[,9]),]

pgk1_3AT_background <- pgk1_3AT_noempty
pgk1_3AT_sig_up <- pgk1_3AT_noempty
pgk1_3AT_sig_down <- pgk1_3AT_noempty
pgk1_3AT_sig_up <- filter(pgk1_3AT_sig_up, pgk1_3AT_sig_up$logFC < -0.5 & pgk1_3AT_sig_up$adj.P.Val < 0.05)
pgk1_3AT_sig_down <- filter(pgk1_3AT_sig_down, pgk1_3AT_sig_down$logFC > 0.5 & pgk1_3AT_sig_down$adj.P.Val < 0.05)

write(pgk1_3AT_background$Yorf1, "/mnt/ingolialab/rmuller1/test_sandbox/gz_to_mpralm/rehearsal_dinner/scripts/GO_analysis_files/pgk1_3AT_background.txt") 
write(pgk1_3AT_sig_up$Yorf1, "/mnt/ingolialab/rmuller1/test_sandbox/gz_to_mpralm/rehearsal_dinner/scripts/GO_analysis_files/pgk1_3AT_sig_up.txt")
write(pgk1_3AT_sig_down$Yorf1, "/mnt/ingolialab/rmuller1/test_sandbox/gz_to_mpralm/rehearsal_dinner/scripts/GO_analysis_files/pgk1_3AT_sig_down.txt")

#GO analysis, generating His4 pooled background and significant up and down text files
his4_pooled_noempty <- his4_pooled
his4_pooled_noempty[his4_pooled_noempty == "Neg_ctrl"] <- NA
his4_pooled_noempty <- his4_pooled_noempty[complete.cases(his4_pooled_noempty[,9]),]

his4_pooled_background <- his4_pooled_noempty
his4_pooled_sig_up <- his4_pooled_noempty
his4_pooled_sig_down <- his4_pooled_noempty
his4_pooled_sig_up <- filter(his4_pooled_sig_up, his4_pooled_sig_up$logFC < -0.5 & his4_pooled_sig_up$adj.P.Val < 0.05)
his4_pooled_sig_down <- filter(his4_pooled_sig_down, his4_pooled_sig_down$logFC > 0.5 & his4_pooled_sig_down$adj.P.Val < 0.05)

write(his4_pooled_background$Yorf1, "/mnt/ingolialab/rmuller1/test_sandbox/gz_to_mpralm/rehearsal_dinner/scripts/GO_analysis_files/his4_pooled_background.txt") 
write(his4_pooled_sig_up$Yorf1, "/mnt/ingolialab/rmuller1/test_sandbox/gz_to_mpralm/rehearsal_dinner/scripts/GO_analysis_files/his4_pooled_sig_up.txt")
write(his4_pooled_sig_down$Yorf1, "/mnt/ingolialab/rmuller1/test_sandbox/gz_to_mpralm/rehearsal_dinner/scripts/GO_analysis_files/his4_pooled_sig_down.txt")

#GO analysis, generating Pgk1 pooled background and significant up and down text files
pgk1_pooled_noempty <- pgk1_pooled
pgk1_pooled_noempty[pgk1_pooled_noempty == "Neg_ctrl"] <- NA
pgk1_pooled_noempty <- pgk1_pooled_noempty[complete.cases(pgk1_pooled_noempty[,9]),]

pgk1_pooled_background <- pgk1_pooled_noempty
pgk1_pooled_sig_up <- pgk1_pooled_noempty
pgk1_pooled_sig_down <- pgk1_pooled_noempty
pgk1_pooled_sig_up <- filter(pgk1_pooled_sig_up, pgk1_pooled_sig_up$logFC < -0.5 & pgk1_pooled_sig_up$adj.P.Val < 0.05)
pgk1_pooled_sig_down <- filter(pgk1_pooled_sig_down, pgk1_pooled_sig_down$logFC > 0.5 & pgk1_pooled_sig_down$adj.P.Val < 0.05)

write(pgk1_pooled_background$Yorf1, "/mnt/ingolialab/rmuller1/test_sandbox/gz_to_mpralm/rehearsal_dinner/scripts/GO_analysis_files/pgk1_pooled_background.txt") 
write(pgk1_pooled_sig_up$Yorf1, "/mnt/ingolialab/rmuller1/test_sandbox/gz_to_mpralm/rehearsal_dinner/scripts/GO_analysis_files/pgk1_pooled_sig_up.txt")
write(pgk1_pooled_sig_down$Yorf1, "/mnt/ingolialab/rmuller1/test_sandbox/gz_to_mpralm/rehearsal_dinner/scripts/GO_analysis_files/pgk1_pooled_sig_down.txt")


#GO analysis, generating HTS1 background and significant up and down text files
HTS1_noempty <- HTS1
HTS1_noempty[HTS1_noempty == "Neg_ctrl"] <- NA
HTS1_noempty <- HTS1_noempty[complete.cases(HTS1_noempty[,9]),]

HTS1_background <- HTS1_noempty
HTS1_sig_up <- HTS1_noempty
HTS1_sig_down <- HTS1_noempty
HTS1_sig_up <- filter(HTS1_sig_up, HTS1_sig_up$logFC < -0.5 & HTS1_sig_up$adj.P.Val < 0.05)
HTS1_sig_down <- filter(HTS1_sig_down, HTS1_sig_down$logFC > 0.5 & HTS1_sig_down$adj.P.Val < 0.05)

write(HTS1_background$Yorf1, "/mnt/ingolialab/rmuller1/test_sandbox/gz_to_mpralm/rehearsal_dinner/scripts/GO_analysis_files/HTS1_background.txt") 
write(HTS1_sig_up$Yorf1, "/mnt/ingolialab/rmuller1/test_sandbox/gz_to_mpralm/rehearsal_dinner/scripts/GO_analysis_files/HTS1_sig_up.txt")
write(HTS1_sig_down$Yorf1, "/mnt/ingolialab/rmuller1/test_sandbox/gz_to_mpralm/rehearsal_dinner/scripts/GO_analysis_files/HTS1_sig_down.txt")


#GO analysis, generating RPC31 background and significant up and down text files
RPC31_noempty <- RPC31
RPC31_noempty[RPC31_noempty == "Neg_ctrl"] <- NA
RPC31_noempty <- RPC31_noempty[complete.cases(RPC31_noempty[,9]),]

RPC31_background <- RPC31_noempty
RPC31_sig_up <- RPC31_noempty
RPC31_sig_down <- RPC31_noempty
RPC31_sig_up <- filter(RPC31_sig_up, RPC31_sig_up$logFC < -0.5 & RPC31_sig_up$adj.P.Val < 0.05)
RPC31_sig_down <- filter(RPC31_sig_down, RPC31_sig_down$logFC > 0.5 & RPC31_sig_down$adj.P.Val < 0.05)

write(RPC31_background$Yorf1, "/mnt/ingolialab/rmuller1/test_sandbox/gz_to_mpralm/rehearsal_dinner/scripts/GO_analysis_files/RPC31_background.txt") 
write(RPC31_sig_up$Yorf1, "/mnt/ingolialab/rmuller1/test_sandbox/gz_to_mpralm/rehearsal_dinner/scripts/GO_analysis_files/RPC31_sig_up.txt")
write(RPC31_sig_down$Yorf1, "/mnt/ingolialab/rmuller1/test_sandbox/gz_to_mpralm/rehearsal_dinner/scripts/GO_analysis_files/RPC31_sig_down.txt")


#GO analysis, generating TF-cds background and significant up and down text files
cds_noempty <- cds
cds_noempty[cds_noempty == "Neg_ctrl"] <- NA
cds_noempty <- cds_noempty[complete.cases(cds_noempty[,9]),]

cds_background <- cds_noempty
cds_sig_up <- cds_noempty
cds_sig_down <- cds_noempty
cds_sig_up <- filter(cds_sig_up, cds_sig_up$logFC < -0.5 & cds_sig_up$adj.P.Val < 0.05)
cds_sig_down <- filter(cds_sig_down, cds_sig_down$logFC > 0.5 & cds_sig_down$adj.P.Val < 0.05)

write(cds_background$Yorf1, "/mnt/ingolialab/rmuller1/test_sandbox/gz_to_mpralm/rehearsal_dinner/scripts/GO_analysis_files/cds_background.txt") 
write(cds_sig_up$Yorf1, "/mnt/ingolialab/rmuller1/test_sandbox/gz_to_mpralm/rehearsal_dinner/scripts/GO_analysis_files/cds_sig_up.txt")
write(cds_sig_down$Yorf1, "/mnt/ingolialab/rmuller1/test_sandbox/gz_to_mpralm/rehearsal_dinner/scripts/GO_analysis_files/cds_sig_down.txt")

#GO analysis, generating TF-utr background and significant up and down text files
utr_noempty <- utr
utr_noempty[utr_noempty == "Neg_ctrl"] <- NA
utr_noempty <- utr_noempty[complete.cases(utr_noempty[,9]),]

utr_background <- utr_noempty
utr_sig_up <- utr_noempty
utr_sig_down <- utr_noempty
utr_sig_up <- filter(utr_sig_up, utr_sig_up$logFC < -0.5 & utr_sig_up$adj.P.Val < 0.05)
utr_sig_down <- filter(utr_sig_down, utr_sig_down$logFC > 0.5 & utr_sig_down$adj.P.Val < 0.05)

write(utr_background$Yorf1, "/mnt/ingolialab/rmuller1/test_sandbox/gz_to_mpralm/rehearsal_dinner/scripts/GO_analysis_files/utr_background.txt") 
write(utr_sig_up$Yorf1, "/mnt/ingolialab/rmuller1/test_sandbox/gz_to_mpralm/rehearsal_dinner/scripts/GO_analysis_files/utr_sig_up.txt")
write(utr_sig_down$Yorf1, "/mnt/ingolialab/rmuller1/test_sandbox/gz_to_mpralm/rehearsal_dinner/scripts/GO_analysis_files/utr_sig_down.txt")

