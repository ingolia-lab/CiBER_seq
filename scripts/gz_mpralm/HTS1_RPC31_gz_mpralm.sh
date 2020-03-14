#!/bin/bash -x

#Start of code, assuming you have the .gz files in a folder in your current working directory
HOME=$(echo ~)
DIR='/CiBER_seq_package/all_raw_fasta_gz/HTS1_RPC31/'
GZDIR1=${HOME}${DIR}/replicate_1/
GZDIR2=${HOME}${DIR}/replicate_3/

rm -f ${GZDIR1}/samples.txt
> ${GZDIR1}/samples.txt
echo ${GZDIR1}*.gz | cut -d ' ' -f1- --output-delimiter=$'\n' | cut -d '.' -f 1 >> ${GZDIR1}/samples.txt
SAMP=${GZDIR1}/samples.txt

for i in $(cat ${SAMP})
do
	if [[ ! -e "${GZDIR1}all_HTS1_rep1_counts.txt" ]]
	then
		gunzip $i.fastq.gz
		fastx_trimmer -f 1 -l 25 -i $i.fastq -o $i-trim25.fastq
		~/CiBER_seq_package/scripts/debug/bc-count --fastq $i-trim25.fastq --output ${i}-count.txt --neighborhood $i
	fi
done

#generating barcode count table from second run (rep3)

rm -f ${GZDIR2}/samples.txt
> ${GZDIR2}/samples.txt
echo ${GZDIR2}*.gz | cut -d ' ' -f1- --output-delimiter=$'\n' | cut -d '.' -f 1 >> ${GZDIR2}/samples.txt
SAMP2=${GZDIR2}/samples.txt

for i in $(cat ${SAMP2})
do
        if [[ ! -e "${GZDIR2}all_HTS1_rep1_counts.txt" ]]
        then
                gunzip $i.fastq.gz
                fastx_trimmer -f 1 -l 25 -i $i.fastq -o $i-trim25.fastq
                ~/CiBER_seq_package/scripts/debug/bc-count --fastq $i-trim25.fastq --output ${i}-count.txt --neighborhood $i
        fi
done

COUNTS1=$(echo ${GZDIR1}*-count.txt)
COUNTS2=$(echo ${GZDIR2}*-count.txt)
ALLCOUNT="${COUNTS1} ${COUNTS2}"
python bc-tabulate.py -o ${HOME}${DIR}all_HTS1_RPC31_counts.txt ${ALLCOUNT}

Rscript ./R_mpralm/HTS1_RPC31_mpralm.r
