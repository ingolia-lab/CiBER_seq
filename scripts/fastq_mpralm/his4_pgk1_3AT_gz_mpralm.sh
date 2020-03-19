#!/bin/bash -x

#Start of code, assuming you have the .gz files in a folder in your current working directory
HOME=$(echo ~)
DIR='/CiBER_seq_package/all_raw_fastq/HIS4_PGK1_3AT/'
GZDIR=${HOME}${DIR}

HIS4='GAGATCCAGTCACTCGGaactgTTAC'
PGK1='GAGATCCAGTCACTCGGtcgatTTAC'

rm -f ${GZDIR}/samples.txt
> ${GZDIR}/samples.txt
echo ${GZDIR}*.fastq | cut -d ' ' -f1- --output-delimiter=$'\n' | cut -d '.' -f 1 >> ${GZDIR}/samples.txt
SAMP=${GZDIR}/samples.txt

for i in $(cat ${SAMP})
do
	if [[ ! -e "${i}his4.fastq.gz" ]]
	then
		cutadapt -a his4="${HIS4}" -a pgk1="${PGK1}" --minimum-length 10 \
        	-o ${i}{name}.fastq ${i}.fastq
	fi
done

PROM='his4 pgk1'

for i in $(cat ${SAMP})
do
	for j in ${PROM}
	do
		if [[ ! -e "${i}${j}count.txt" ]]
        	then
			~/CiBER_seq/barcode_assign/bc-count --fastq $i$j.fastq --output $i${j}-count.txt --neighborhood $i${j}
		fi
	done
done

HISCOUNTS=$(echo ${GZDIR}*his4-count.txt)
PGKCOUNTS=$(echo ${GZDIR}*pgk1-count.txt)

if [[ ! -e "${GZDIR}all_his4_seq1_counts.txt" ]]
then
	python bc-tabulate.py -o ${GZDIR}all_his4_seq1_counts.txt ${HISCOUNTS}
	python bc-tabulate.py -o ${GZDIR}all_pgk1_seq1_counts.txt ${PGKCOUNTS}
fi

Rscript ./R_mpralm/his4_3AT_mpralm.r

