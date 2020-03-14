#!/bin/bash -x

HOME=$(echo ~)
DIR='/CiBER_seq_package/all_raw_fasta_gz/bc_validation/'
GZDIR=${HOME}${DIR}

cutadapt --discard-untrimmed --minimum-length 25 --maximum-length 25 -m 10 -a GGAGTTATCAACAGAGCGCGAGATCG \
	-o ${GZDIR}/DNA_rep1.trimmed.fastq.gz ${GZDIR}/DNA_rep1_S118_L004_R1_001.fastq.gz > ${GZDIR}/DNA_rep1.cutadapt.log
cutadapt --discard-untrimmed --minimum-length 25 --maximum-length 25 -m 10 -a GGAGTTATCAACAGAGCGCGAGATCG \
	-o ${GZDIR}/DNA_rep2.trimmed.fastq.gz ${GZDIR}/DNA_rep2_S119_L004_R1_001.fastq.gz > ${GZDIR}/DNA_rep2.cutadapt.log
cutadapt --discard-untrimmed --minimum-length 25 --maximum-length 25 -m 10 -a GGAGTTATCAACAGAGCGCGAGATCG \
	-o ${GZDIR}/RNA_rep1.trimmed.fastq.gz ${GZDIR}/RNA_rep1_S120_L004_R1_001.fastq.gz > ${GZDIR}/RNA_rep1.cutadapt.log
cutadapt --discard-untrimmed --minimum-length 25 --maximum-length 25 -m 10 -a GGAGTTATCAACAGAGCGCGAGATCG \
	-o ${GZDIR}/RNA_rep2.trimmed.fastq.gz ${GZDIR}/RNA_rep2_S121_L004_R1_001.fastq.gz > ${GZDIR}/RNA_rep2.cutadapt.log

gunzip ${GZDIR}/DNA_rep1.trimmed.fastq.gz
gunzip ${GZDIR}/DNA_rep2.trimmed.fastq.gz
gunzip ${GZDIR}/RNA_rep1.trimmed.fastq.gz
gunzip ${GZDIR}/RNA_rep2.trimmed.fastq.gz

python ~/CiBER_seq_package/scripts/est_barcodes/tabulate_bc_est.py


