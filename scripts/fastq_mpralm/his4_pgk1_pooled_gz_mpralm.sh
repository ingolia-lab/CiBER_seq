#!/bin/bash -x

#Start of code, assuming you have the .gz files in a folder in your current working directory
HOME=$(echo ~)
DIR='/CiBER_seq_package/all_raw_fastq/HIS4_PGK1_pooled/'
GZDIR=${HOME}${DIR}/
#GZDIR2=${HOME}${DIR}/more_seq/

HIS4='GAGATCCAGTCACTCGGaactgTTAC'
PGK1='GAGATCCAGTCACTCGGtcgatTTAC'

PROM='his4 pgk1'

rm -f ${GZDIR}/samples.txt
> ${GZDIR}/samples.txt
echo ${GZDIR}*.fastq | cut -d ' ' -f1- --output-delimiter=$'\n' | cut -d '.' -f 1 >> ${GZDIR}/samples.txt
SAMP=${GZDIR}/samples.txt

for i in $(cat ${SAMP})
do
        if [[ ! -e "${GZDIR}all_his4_moreseq_counts.txt" ]]
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
                ../barcode_assign/target/debug/bc-count --fastq $i$j.fastq --output $i${j}-count.txt --neighborhood $i${j}
        done
done

HISCOUNTS=$(echo ${GZDIR}*his4-count.txt)
PGKCOUNTS=$(echo ${GZDIR}*pgk1-count.txt)

python bc-tabulate.py -o ${GZDIR}all_his4_pooled_counts.txt $HISCOUNTS
python bc-tabulate.py -o ${GZDIR}all_pgk1_pooled_counts.txt $PGKCOUNTS

Rscript ./R_mpralm/his4_pooled_mpralm.r
