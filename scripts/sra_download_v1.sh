#!/bin/bash

#defining directiories
HOME=$(echo ~)
DIR='/CiBER_seq/scripts/sra_files/'
SRADIR=${HOME}${DIR}
GZ='/CiBER_seq_package/all_raw_fasta_gz/'
GZDIR=${HOME}${GZ}

#loop over list of SRA numbers in experiment files
bc_validation=${SRADIR}/bc_val.txt
GCN4_CDS_UTR=${SRADIR}/CDS_UTR.txt
HIS4_PGK1_pooled=${SRADIR}/HIS4_PGK1.txt
PE_bc_gRNA_assignment=${SRADIR}/bc_gRNA.txt
HIS4_PGK1_3AT=${SRADIR}/HIS4_PGK1_3AT.txt
HTS1_RPC31=${SRADIR}/HTS1_RPC31.txt

#Put them into directory names that match with existing code

for i in $(cat ${bc_validation})
do
	fastq-dump -O ${GZDIR}/bc_validation/${i}
done

for i in $(cat ${HIS4_PGK1_pooled})
do
        fastq-dump -O ${GZDIR}/HIS4_PGK1_pooled/${i}
done

for i in $(cat ${GCN4_CDS_UTR})
do
        fastq-dump -O ${GZDIR}/GCN4_CDS_UTR/${i}
done

fastq-dump -O ${GZDIR}/PE_bc_gRNA_assignment/mod_bc_gRNA_R1.fastq.gz ${GZDIR}/PE_bc_gRNA_assignment/mod_bc_gRNA_R2.fastq.gz --split-files SRR10327353

for i in $(cat ${HIS4_PGK1_3AT})
do
        fastq-dump -O ${GZDIR}/HIS4_PGK1_3AT/${i}
done

for i in $(cat ${HTS1_RPC31})
do
        fastq-dump -O ${GZDIR}/HTS1_RPC31/${i}
done
