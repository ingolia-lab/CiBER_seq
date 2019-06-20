# CiBER_seq
Code to make all the datasets and analysis for Amino-acid starvation CiBER_seq project

This repository contains code to generate the following datasets:
  -PGK1 vs His4 promoter
  -His4 promoter in background gRNA knockdown of HTS1 or RPC31
  -Synthetic Transcription Factor with promoter + 5'UTR of GCN4 and GCN4_CDS-fused TF
  -barcode counting validation sequencing experiment
  
The starting gzipped files are in rmuller1@compute1:/mnt/ingolialab/rmuller1/CiBER-seq_paper/all_raw_fasta_gz/

Each dataset has a bash script associated with it. The bash script will take the starting gzipped files, unzip them, use a combination of fastx_trimmer and cut_adapt to assign 5nt identifiers to the sample type and trim the barcode reads to 25nt. The bash script then uses bc-count.py to count the number of times a barcode appears and bc-tabulate.py to organize the separate bc-count files into a dataframe matrix.

The dataframe matrix can then be analyzed in R. Each dataset has an associated R script that will take the bc-tabulate matrix, filter cut-off for barcodes with an IVT-pre count >32 reads, and perform mpralm analysis (Kasper Hansen) to identify barcodes whose RNA:DNA ratio significantly changes in response to guide RNA induction. 

Separately, the IVT only samples from His4/Pgk1 dataset (treated as the WT condition) and HTS1/RPC31 gRNA KD dataset (treated as the KD condition) were analyzed in DESeq2 as a growth screen to look for synthetic lethal guides and generally guides that affect growth. 


