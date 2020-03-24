# CiBER_seq

#This repository is organized as follows:
  -Code for generating barcode count tables and violin plots (figure 1) is found in /CiBER_seq/scripts/est_barcodes
  -Barcode gRNA assignment code is located in /CiBER_seq/barcode_assign
  -Code that processes raw sequencing .fastq files into barcode count tables is found in /CiBER_seq/scripts/fastq_mpralm
  -Code that processes count tables to perform mpralm analysis is found in /CiBER_seq/scripts/R_mpralm
  -Code that uses mpralm output data tables to graph volcano plots, scatter plots, and other paper figures is located in /CiBER_seq/scripts/R_figures
  -Hand curated gene lists by biological function are found in /CiBER_seq/scripts/GO_analysis_files/GO_annotation_lists
  -Additional code and reference .txt files for SRA download, pipeline analysis, and count table organization is found in /CiBER_seq/scripts
  
#How the analysis pipeline works:
All script directory callings within the code assumes the CiBER_seq repository has been downloaded into a folder in the home directory called "CiBER_github". The analysis pipeline is performed by running the master.sh script from the /CiBER_seq/scripts/ directory. The master.sh script calls the sra_download.sh script which makes a new directory ~/CiBER_seq_package/all_raw_fastq/ and downloads all the sra datasets to this location. The master.sh script then calls each dataset analysis pipeline. These pipelines trim adapters from the raw .fastq files, generate barcode counts, merge the barcode counts for each sample into a data matrix by experiment, and generate: 
  -Barcode count tables for figure 1
  -mpralm data tables for:
    -P(HIS4) CiBER-seq
    -P(PGK1) CiBER-seq
    -P(HIS4) 3AT treatment CiBER-seq
    -P(HIS4) Epistasis during 3AT treatment, HTS1 knockdown, and RPC31 knockdown
    -P(GEM synthetic) with GCN4_CDS-transcription factor, GCN4 post-translational CiBER-seq
    -P(GEM synthetic) with GCN4_UTR-transcription factor, GCN4 translation control CiBER-seq
The barcode count tables and mpralm data tables are written to their corresponding directories labeled by experiment. They can be found in subdirectories of ~/CiBER_seq_package/all_raw_fastq/


