#!/bin/bash -x

#Start of code, assuming you have the .gz files in a folder in your current working directory
HOME=$(echo ~)
DIR='/CiBER_seq_package/all_raw_fasta_gz/HIS4_PGK1_pooled/'
GZDIR=${HOME}${DIR}/
GZDIR2=${HOME}${DIR}/more_seq/

HIS4='GAGATCCAGTCACTCGGaactgTTAC'
PGK1='GAGATCCAGTCACTCGGtcgatTTAC'

PROM='his4 pgk1'

rm -f ${GZDIR2}/samples.txt
> ${GZDIR2}/samples.txt
echo ${GZDIR2}*.gz | cut -d ' ' -f1- --output-delimiter=$'\n' | cut -d '.' -f 1 >> ${GZDIR2}/samples.txt
SAMP2=${GZDIR2}/samples.txt

for i in $(cat ${SAMP2})
do
        if [[ ! -e "${GZDIR}all_his4_moreseq_counts.txt" ]]
        then
                cutadapt -a his4="${HIS4}" -a pgk1="${PGK1}" --minimum-length 10 \
                -o ${i}{name}.fastq.gz ${i}.fastq.gz
        fi
done

PROM='his4 pgk1'

for i in $(cat ${SAMP2})
do
        for j in ${PROM}
        do
                gunzip $i$j.fastq.gz
                ~/CiBER_seq_package/scripts/debug/bc-count --fastq $i$j.fastq --output $i${j}-count.txt --neighborhood $i${j}
        done
done

HISCOUNTS2=$(echo ${GZDIR2}*his4-count.txt)
PGKCOUNTS2=$(echo ${GZDIR2}*pgk1-count.txt)

python bc-tabulate.py -o ${GZDIR2}all_his4_moreseq_counts.txt $HISCOUNTS2
python bc-tabulate.py -o ${GZDIR2}all_pgk1_moreseq_counts.txt $PGKCOUNTS2

Rscript ./R_mpralm/his4_pooled_mpralm.r
