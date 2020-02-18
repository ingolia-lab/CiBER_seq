#Stacked barplot for 3AT, HTS1, and RPC31 

#Saturated
saturated <- filter(HTS1_act_no_ind, HTS1_act_no_ind$RPC31_adj.P.Val < 0.05 & HTS1_act_no_ind$RPC31_logFC < -0.5)
saturated <- unique(saturated$Yorf1)

sat_aabiosyn <- subset(saturated, (saturated %in% aabiosynthesis_handcur$name))
sat_polIII <- subset(saturated, (saturated %in% polIIIsubunit_handcur$name))
sat_transl <- subset(saturated, (saturated %in% translationcontrol_handcur$name))
sat_tRNAsynth <- subset(saturated, (saturated %in% tRNAaasynthetase_handcur$name))
sat_tRNAproc <- subset(saturated, (saturated %in% tRNAprocessing_handcur$name))
'%notin%' <- Negate('%in%')
sat_other <- saturated
sat_other <- subset(sat_other, (sat_other %notin% aabiosynthesis_handcur$name))
sat_other <- subset(sat_other, (sat_other %notin% polIIIsubunit_handcur$name))
sat_other <- subset(sat_other, (sat_other %notin% translationcontrol_handcur$name))
sat_other <- subset(sat_other, (sat_other %notin% tRNAaasynthetase_handcur$name))
sat_other <- subset(sat_other, (sat_other %notin% tRNAprocessing_handcur$name))
sat_other <- subset(sat_other, (sat_other != "NA"))
  
categories <- c("amino acid biosyn", "tRNA_synth", "translation_ctrl", "RNA_pol_3", "tRNA_proc", "other")
3AT_sat <- c(16, 25, 4, 16, 9, 44)
HTS1_sat <- c(27, 39, 19, 26, 19, 123)
RPC31_sat <- c(26, 35, 15, 25, 16, 90)

#Blocked
categories <- c("ISR + pol3 + tRNA_proc", "actin", "ER_trafficking", 
                "ribosomal_subunits", "proteasome", "other")
3AT_block <- c()
HTS1_block <- c()
RPC31_block <- c()

# saturated stacked bar graph (by targeting guides)
library(ggplot2)
# create a dataset
specie <- c(rep("3AT" , 6) , rep("HTS1" , 6) , rep("RPC31" , 6))
categories <- rep(c("bamino acid biosyn", "btRNA_synth", "ctranslation_ctrl", 
                    "dRNA_pol_3", "etRNA_proc", "aother"), 3)
thrAT_sat <- c(16, 25, 4, 16, 9, 44)
HTS1_sat <- c(27, 39, 19, 26, 19, 123)
RPC31_sat <- c(26, 35, 15, 25, 16, 90)
value <- c(16, 25, 4, 16, 9, 44, 27, 39, 19, 26, 19, 123, 26, 35, 15, 25, 16, 90)
data <- data.frame(specie,categories,value)
# Stacked
ggplot(data, aes(fill=categories, y=value, x=specie)) + 
  geom_bar(position="stack", stat="identity")

# saturated stacked bar graph (by targeted genes)
library(ggplot2)
# create a dataset
specie <- c(rep("3AT" , 6) , rep("HTS1" , 6) , rep("RPC31" , 6))
categories <- rep(c("bamino acid biosyn", "btRNA_synth", "ctranslation_ctrl", 
                    "dRNA_pol_3", "etRNA_proc", "aother"), 3)
thrAT_sat <- c(9, 11, 3, 8, 6, 37)
HTS1_sat <- c(10, 16, 6, 9, 9, 95)
RPC31_sat <- c(11, 15, 6, 9, 9, 69)
value <- c(9, 11, 3, 8, 6, 37, 10, 16, 6, 9, 9, 95, 11, 15, 6, 9, 9, 69)
data <- data.frame(specie,categories,value)
# Stacked
ggplot(data, aes(fill=categories, y=value, x=specie)) + 
  geom_bar(position="stack", stat="identity")


# blocked stacked bar graph (by targeted genes)

blocked <- filter(HTS1_block, HTS1_block$HTS1_adj.P.Val < 0.05 & HTS1_block$HTS1_logFC < -0.5)
blocked <- unique(blocked$Yorf1)

actin_block <- subset(blocked, (blocked %in% actin_extra_handcur$name))
srp_block <- subset(blocked, (blocked %in% srp_handcur$name))
ribo_sub_block <- subset(blocked, (blocked %in% ribosomal_subs_handcur$name))
protea_block <- subset(blocked, (blocked %in% proteasome_handcur$name))
'%notin%' <- Negate('%in%')
block_other <- blocked
block_other <- subset(block_other, (block_other %notin% actin_extra_handcur$name))
block_other <- subset(block_other, (block_other %notin% srp_handcur$name))
block_other <- subset(block_other, (block_other %notin% ribosomal_subs_handcur$name))
block_other <- subset(block_other, (block_other %notin% proteasome_handcur$name))
block_other <- subset(block_other, (block_other != "NA"))

library(ggplot2)
# create a dataset
specie <- c(rep("3AT" , 6) , rep("HTS1" , 6) , rep("RPC31" , 6))
categories <- rep(c("bactin", "signal_rec_particle", 
                    "ribosomal_subunits", "proteasome", "aother"), 3)
thrAT_sat <- c(1, 0, 1, 0, 64)
HTS1_sat <- c(9, 0, 13, 9, 300)
#RPC31_sat <- c(11, 15, 6, 9, 9, 69)
#value <- c(9, 11, 3, 8, 6, 37, 10, 16, 6, 9, 9, 95, 11, 15, 6, 9, 9, 69)
data <- data.frame(specie,categories,value)
# Stacked
ggplot(data, aes(fill=categories, y=value, x=specie)) + 
  geom_bar(position="stack", stat="identity")

