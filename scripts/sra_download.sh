#!/bin/bash

#defining directories
HOME=$(echo ~)
DIR='/CiBER_github/CiBER_seq/scripts/'
SRADIR=${HOME}${DIR}
GZ='/CiBER_seq_package/all_raw_fastq/'
GZDIR=${HOME}${GZ}

mkdir ${GZDIR}/bc_validation/
mkdir ${GZDIR}/GCN4_CDS_UTR/
mkdir ${GZDIR}/HIS4_PGK1_pooled/
mkdir ${GZDIR}/PE_bc_gRNA_assignment/
#mkdir ${GZDIR}/HIS4_PGK1_3AT/
mkdir ${GZDIR}/HTS1_RPC31/
mkdir ${GZDIR}/epistasis_3AT_HTS1_RPC31/

names=${SRADIR}/sra_filenames.txt

cat ${names} | while read field1 field2;
do
	fastq-dump -O ${GZDIR} "$field1"
	mv ${GZDIR}/"$field1".fastq ${GZDIR}/"$field2".fastq
done

fastq-dump -O ${GZDIR}/PE_bc_gRNA_assignment/ --split-files SRR10327353
mv ${GZDIR}/SRR10327353.fastq ${GZDIR}/mod_bc_gRNA_R1.fastq
mv ${GZDIR}/SRR10327353_1.fastq ${GZDIR}/mod_bc_gRNA_R2.fastq

mv ${GZDIR}/*rep*.fastq ${GZDIR}/bc_validation/
mv ${GZDIR}/*CU*.fastq ${GZDIR}/GCN4_CDS_UTR/
mv ${GZDIR}/PH_*.fastq ${GZDIR}/HIS4_PGK1_pooled/
mv ${GZDIR}/mod*.fastq ${GZDIR}/PE_bc_gRNA_assignment/

cp -r ${GZDIR}/HIS4_PGK1_pooled/ ${GZDIR}/HIS4_PGK1_3AT/
mv ${GZDIR}/*3AT*.fastq ${GZDIR}/HIS4_PGK1_3AT/

mkdir ${GZDIR}/HTS1_RPC31/HTS1/
mkdir ${GZDIR}/HTS1_RPC31/RPC31/
mv ${GZDIR}/*HTS1*.fastq ${GZDIR}/HTS1_RPC31/HTS1
mv ${GZDIR}/*RPC31*.fastq ${GZDIR}/HTS1_RPC31/RPC31/
