#!/bin/bash -x

#make fastq data directory
mkdir ~/CiBER_seq_package/
mkdir ~/CiBER_seq_package/all_raw_fastq/

DATADIR='~/CiBER_seq_package/all_raw_fastq/'

#downloads sra files and puts them into corresponding folders
nohup bash sra_download.sh

#Generate count tables for barcode validation experiment (fig1)
nohup bash est_barcodes/bc_est_count.sh

#Takes PE sequencing data to make the bc-gRNA assignments (also includes empty barcodes for FDR calculations)
PEDIR='~/CiBER_seq_package/all_raw_fastq/PE_bc_gRNA_assignment/'
if [[ ! -e "${PEDIR}/grna-assign-barcode-grna-good.txt" ]]
then
    ./gRNA_assign_plus_empty.sh
fi

#Analysis for His4/Pgk1 3AT dataset (fig3) Takes sra fastq and outputs mpralm analysis dataset
if [[ ! -e "${DATADIR}/HIS4_PGK1_3AT/his4_3AT_sum_mpralm.txt" ]]
then
    nohup bash fastq_mpralm/his4_pgk1_3AT_gz_mpralm.sh
fi

#Analysis for His4/Pgk1 pooled dataset (fig2) Takes sra fastq and outputs mpralm analysis dataset
if [[ ! -e "${DATADIR}/HIS4_PGK1_pooled/his4_pooled_sum_mpralm.txt" ]]
then
    nohup bash fastq_mpralm/his4_pgk1_pooled_gz_mpralm.sh
fi

#Analysis for HTS1/RPC31 dataset (fig4) Takes sra fastq and outputs mpralm analysis dataset
if [[ ! -e "${DATADIR}/HTS1_RPC31/HTS1_sum_mpralm.txt" ]]
then
    nohup bash fastq_mpralm/HTS1_RPC31_gz_mpralm.sh
fi

#Analysis for Epistasis: 3AT, HTS1, and RPC31 (fig4)
Rscript R_mpralm/epistasis_3AT_HTS1_RPC31.r

#Analysis for TF-CDS/TF-UTR dataset (fig5, fig6) Takes sra fastq and outputs mpralm analysis dataset
if [[ ! -e "${DATADIR}/GCN4_CDS.UTR/cds_sum_mpralm.txt" ]]
then
    nohup bash fastq_mpralm/TF_CDS_UTR_gz_mpralm.sh
fi
