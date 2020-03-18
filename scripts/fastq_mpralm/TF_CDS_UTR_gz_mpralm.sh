#!/bin/bash -x

#Start of code, assuming you have the .gz files in a folder in your current working directory
HOME=$(echo ~)
DIR='/CiBER_seq_package/all_raw_fastq/GCN4_CDS_UTR/'
GZDIR=${HOME}${DIR}

CDS='GAGATCCAGTCACTCGGCCTTGCTAT'
UTR='GAGATCCAGTCACTCGGTCAGACTAT'

rm -f ${GZDIR}/samples.txt
> ${GZDIR}/samples.txt
echo ${GZDIR}*.gz | cut -d ' ' -f1- --output-delimiter=$'\n' | cut -d '.' -f 1 >> ${GZDIR}/samples.txt
SAMP=${GZDIR}/samples.txt

for i in $(cat ${SAMP})
do
	if [[ ! -e "${GZDIR}all_TF_CDS_counts.txt" ]]
	then
		cutadapt -a cds="${CDS}" -a utr="${UTR}" --minimum-length 10 \
        	-o ${i}{name}.fastq.gz ${i}.fastq.gz
	fi
done

PROM='cds utr'

for i in $(cat ${SAMP})
do
	for j in ${PROM}
	do
		gunzip $i$j.fastq.gz
		~/CiBER_seq_package/scripts/debug/bc-count --fastq $i$j.fastq --output $i${j}-count.txt --neighborhood $i${j}
	done
done

CDSCOUNTS=$(echo ${GZDIR}*cds-count.txt)
UTRCOUNTS=$(echo ${GZDIR}*utr-count.txt)

python bc-tabulate.py -o ${GZDIR}all_TF_CDS_counts.txt ${CDSCOUNTS}
python bc-tabulate.py -o ${GZDIR}all_TF_UTR_counts.txt ${UTRCOUNTS}

Rscript ./R_mpralm/TF_CDS_UTR_mpralm.r

