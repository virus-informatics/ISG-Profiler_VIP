#!/bin/sh
#$ -S /bin/bash
###a #$ -l lmem
#$ -l arm
#$ -o log
#$ -e log
#$ -cwd
#$ -l s_vmem=50G
#$ -pe def_slot 2

export OPENBLAS_NUM_THREADS=1
module load R/4.3.2
NCBI_EMAIL="your email address"

# input ISG & control gene symbol lists
Rscript 01_get_gene_ID_QC.R \
    isg_60.list \
    cntl_100.list \
    "$NCBI_EMAIL" \
    ISGcntl_geneID_260904.txt || exit 1

# get Isoform list
Rscript 02_get_ISGcntl.R \
    ISGcntl_geneID_260904.txt \
    gene2refseq_Amniota_ISGcntl_260904.list \
    Isoform_260904.list || exit 1

set -euo pipefail

# get Isoform fasta
python3 03_get_fasta.py \
    Isoform_260904.list \
    Isoform_260904.fasta \
    --email "$NCBI_EMAIL" \
    --batch-size 200 \
    2> get_fasta_260904.log || exit 1

echo "FASTA download completed."

echo -n "Input accessions: "
wc -l < Isoform_260904.list

echo -n "FASTA records: "
grep -c '^>' Isoform_260904.fasta

# salmon indexing
salmon_cmd="/usr/local/package/salmon/1.10.1/bin/salmon"

${salmon_cmd} index -p 8 -t Isoform_260904.fasta -i Isoform_260904_salmon &>> index.log || exit 1
