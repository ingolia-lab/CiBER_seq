#!/bin/bash -x

# Adapter sequences for barcode (3' adapter) and guide (both sides)
BCTRIM=GAGATCCAGTCACTCGGGATCCgatctgccaattgaacataacatggtagt
GTRIM3=aagttaaaataaggct

BCDIR=../barcode_assign/target/debug/
DATADIR=~/CiBER_seq_package/all_raw_fastq/PE_bc_gRNA_assignment/
#DATADIR=~/bc_transcript_CRISPRi/updated_bc_gRNA_assignment/PE_guide_alignments/

# wget http://hannonlab.cshl.edu/fastx_toolkit/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2
# tar -xjf fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2
# sudo cp ./bin/* /usr/local/bin

#gunzip ${DATADIR}/mod_bc_gRNA_R1.fastq.gz
#gunzip ${DATADIR}/mod_bc_gRNA_R2.fastq.gz
fastx_reverse_complement -i ${DATADIR}/mod_bc_gRNA_R2.fastq -o ${DATADIR}/revcomp_mod_bc_gRNA_R2.fastq

# Trim barcode adapter
if [[ ! -e "${DATADIR}/bctrim_R1.fq" ]]
then
    cutadapt -a "${BCTRIM}" --pair-filter=both \
	   --untrimmed-output="${DATADIR}/no_bctrim_R1.fq" \
	   --untrimmed-paired-output="${DATADIR}/no_bctrim_R2.fq" \
	   --interleaved \
	   ${DATADIR}/mod_bc_gRNA_R1.fastq \
	   ${DATADIR}/revcomp_mod_bc_gRNA_R2.fastq \
	   2> "${DATADIR}/bctrim-adapter.txt" \
        | cutadapt --interleaved \
	         --pair-filter=any --minimum-length=12 \
	         --too-short-output="${DATADIR}/short_bctrim_R1.fq" \
	         --too-short-paired-output="${DATADIR}/short_bctrim_R2.fq" \
	         -o "${DATADIR}/bctrim_R1.fq" -p "${DATADIR}/bctrim_R2.fq" \
	         - \
	         > "${DATADIR}/bctrim-length.txt"  2>&1
fi

# Use barcode (R1) sequence as name for guide (R2) sequence
if [[ ! -e "${DATADIR}/grna_barcode_barcoded.fq" ]]
then
    ${BCDIR}/bc-seqs \
        --barcodes "${DATADIR}/bctrim_R1.fq" \
        --sequences "${DATADIR}/bctrim_R2.fq" \
        --outbase "${DATADIR}/grna_barcode"
fi

# Trim barcoded guide sequences
if [[ ! -e "${DATADIR}/grna_barcode_trimmed.fq" ]]
then
    cutadapt -a "${GTRIM3}" \
	   --untrimmed-output="${DATADIR}/no_grnatrim.fq" \
	   --minimum-length=20 \
	   --too-short-output="${DATADIR}/short_grnatrim.fq" \
	   "${DATADIR}/grna_barcode_barcoded.fq" \
	   -o "${DATADIR}/grna_barcode_trimmed.fq" \
	   > "${DATADIR}/grnatrim-adapter.txt" 2>&1
fi

# Generate bowtie index of guide sequences
if [[ ! -e "${DATADIR}/target-oligos.1.bt2" ]]
then
    grep -v Yorf "./target-oligos.txt" \
        | sort -k3,3 | uniq -f 2 \
        | awk '{print ">" $1 "\n" $3}' \
	    > "${DATADIR}/target-oligos.fa"
    echo ">No_gRNA" >> "${DATADIR}/target-oligos.fa"
    echo "gctgggaacgaaactctgggagctgcgattggcagCCTAGgttttagagctagaaatagc" >> "${DATADIR}/target-oligos.fa"
    bowtie2-build "${DATADIR}/target-oligos.fa" "${DATADIR}/target-oligos"
fi

# Align guides against bowtie index
if [[ ! -e "${DATADIR}/grna_barcode.bam" ]]
then
    bowtie2 -p38 \
	  -x "${DATADIR}/target-oligos" \
	  -U "${DATADIR}/grna_barcode_trimmed.fq" \
	  2> "${DATADIR}/grna-bowtie.txt" \
        | samtools view -b -o "${DATADIR}/grna_barcode.bam" - \
	         2> "${DATADIR}/grna-bowtie-samtools.txt"
fi

# Sort guide sequence alignments by name
if [[ ! -e "${DATADIR}/grna_barcode_sorted.bam" ]]
then
    samtools sort -m 4G -n -o "${DATADIR}/grna_barcode_sorted.bam" "${DATADIR}/grna_barcode.bam"
fi

# Determine barcode-to-guide assignments
if [[ ! -e "${DATADIR}/grna-assign-barcode-fates.txt" ]]
then
    ${BCDIR}/bc-grna \
        --align-start 0 \
        --bam-by-name "${DATADIR}/grna_barcode_sorted.bam" \
        --outbase "${DATADIR}/grna-assign-"

    cut -f1 "${DATADIR}/grna-assign-barcode-grna-good.txt" \
        "${DATADIR}/grna-assign-barcode-bad-assign.txt" \
        > "${DATADIR}/barcodes-known.txt"

#    grep -F -f "${DATADIR}/barcodes-known.txt" \
#         "${DATADIR}/rym.txt" \
#         > rym-known.txt
#
#    grep -v -F -f "${DATADIR}/barcodes-known.txt" \
#         "${DATADIR}/rym.txt" \
#         > rym-unknown.txt
fi
