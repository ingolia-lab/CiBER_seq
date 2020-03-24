#!/bin/bash -x

HOME=$(echo ~)
DIR='/CiBER_seq_package/all_raw_fastq/bc_validation/'
GZDIR=${HOME}${DIR}

cutadapt --discard-untrimmed --minimum-length 25 --maximum-length 25 -m 10 -a GGAGTTATCAACAGAGCGCGAGATCG \
	-o ${GZDIR}/DNA_rep1.trimmed.fastq ${GZDIR}/DNA_rep1.fastq > ${GZDIR}/DNA_rep1.cutadapt.log
cutadapt --discard-untrimmed --minimum-length 25 --maximum-length 25 -m 10 -a GGAGTTATCAACAGAGCGCGAGATCG \
	-o ${GZDIR}/DNA_rep2.trimmed.fastq ${GZDIR}/DNA_rep2.fastq > ${GZDIR}/DNA_rep2.cutadapt.log
cutadapt --discard-untrimmed --minimum-length 25 --maximum-length 25 -m 10 -a GGAGTTATCAACAGAGCGCGAGATCG \
	-o ${GZDIR}/RNA_rep1.trimmed.fastq ${GZDIR}/RNA_rep1.fastq > ${GZDIR}/RNA_rep1.cutadapt.log
cutadapt --discard-untrimmed --minimum-length 25 --maximum-length 25 -m 10 -a GGAGTTATCAACAGAGCGCGAGATCG \
	-o ${GZDIR}/RNA_rep2.trimmed.fastq ${GZDIR}/RNA_rep2.fastq > ${GZDIR}/RNA_rep2.cutadapt.log

pip install biopython
python ${HOME}/CiBER_github/CiBER_seq/scripts/est_barcodes/tabulate_bc_est.py


