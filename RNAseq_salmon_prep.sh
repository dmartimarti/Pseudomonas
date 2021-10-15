## MAIN STEPS FOR RNA SEQUENCING

## Project: Pseudomonas project from Jen and Filipe

## First steps done in the HPC
# analysis directory: /rds/general/user/dmarti14/home/Pseudomonas

module load anaconda3/personal
source activate salmon_final

cd /rds/general/user/dmarti14/home/Pseudomonas

# this will donwload the cDNA data from C. elegans from the Ensemblgenomes website 
# this is the target transcriptome
# https://www.ensembl.org/Caenorhabditis_elegans/Info/Index
# release-104 
curl http://ftp.ensembl.org/pub/release-104/fasta/caenorhabditis_elegans/cdna/Caenorhabditis_elegans.WBcel235.cdna.all.fa.gz -o celegans.fa.gz

# generate the index of the target transcriptome
salmon index -t celegans.fa.gz -i celegans_index

### IMPORTANT
# to process the data, use the scripts available in the repo
# 1. trimm_seqs.pbs to quality filter
# 2. salmon_hpc.pbs to process the filtered seqs into counts with Salmon

## Copy the quants folder into my computer
cd /mnt/d/MRC_Postdoc/Pseudomonas/RNAseq_pseudomonas
# scp to my computer
scp -r dmarti14@login.hpc.imperial.ac.uk:/rds/general/user/dmarti14/home/Pseudomonas/quants ./

### THE NEXT PARTS OF THE ANALYSIS ARE DONE IN R WITH DESEQ2
