#!/bin/bash -x

#From the PE sequencing data, make the bc-gRNA assignments (also includes empty barcodes for FDR calculations)
#Download PE sequencing .gz files from sra
PEDIR='~/CiBER_seq_package/all_raw_fasta_gz/PE_bc_gRNA_assignment/'
if [[ ! -e "${PEDIR}/grna-assign-barcode-grna-good.txt" ]]
then
    ./gRNA_assign_plus_empty.sh
fi

#Analysis for His4/Pgk1 3AT dataset
#Download from sra and place .gz files into a folder named his4_3AT
DATADIR='~/CiBER_seq_package/all_raw_fasta_gz/'

if [[ ! -e "${DATADIR}/HIS4_PGK1_3AT/his4_3AT_sum_mpralm.txt" ]]
then
    nohup bash gz_mpralm/his4_pgk1_3AT_gz_mpralm.sh
fi

#Analysis for His4/Pgk1 pooled dataset
if [[ ! -e "${DATADIR}/HIS4_PGK1_pooled/his4_pooled_sum_mpralm.txt" ]]
then
    nohup bash gz_mpralm/his4_pgk1_pooled_gz_mpralm.sh
fi

#Analysis for HTS1/RPC31 dataset
if [[ ! -e "${DATADIR}/HTS1_RPC31/HTS1_sum_mpralm.txt" ]]
then
    nohup bash gz_mpralm/HTS1_RPC31_gz_mpralm.sh
fi

#Analysis for TF-CDS/TF-UTR dataset
if [[ ! -e "${DATADIR}/GCN4_CDS.UTR/cds_sum_mpralm.txt" ]]
then
    nohup bash gz_mpralm/TF_CDS_UTR_gz_mpralm.sh
fi
