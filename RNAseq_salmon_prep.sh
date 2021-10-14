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

WORK=/rds/general/project/lms-cabreiro-raw/live/211013_VH00504_11_AAAMNH7M5/Unaligned

# bioF mutant control
salmon quant -i celegans_index -l A -1 $WORK/290921_JV_B0L1_S5_L001_R1_001.fastq.gz  -2 $WORK/290921_JV_B0L1_S5_L001_R2_001.fastq.gz  -p 8 --validateMappings --gcBias -o quants/B0_1
salmon quant -i celegans_index -l A -1 $WORK/290921_JV_B0L2_S12_L001_R1_001.fastq.gz -2 $WORK/290921_JV_B0L2_S12_L001_R2_001.fastq.gz -p 8 --validateMappings --gcBias -o quants/B0_2
salmon quant -i celegans_index -l A -1 $WORK/290921_JV_B0L3_S19_L001_R1_001.fastq.gz -2 $WORK/290921_JV_B0L3_S19_L001_R2_001.fastq.gz -p 8 --validateMappings --gcBias -o quants/B0_3
salmon quant -i celegans_index -l A -1 $WORK/290921_JV_B0L4_S26_L001_R1_001.fastq.gz -2 $WORK/290921_JV_B0L4_S26_L001_R2_001.fastq.gz -p 8 --validateMappings --gcBias -o quants/B0_4

# bioF mutant Metformin
salmon quant -i celegans_index -l A -1 $WORK/290921_JV_B50L1_S6_L001_R1_001.fastq.gz  -2 $WORK/290921_JV_B50L1_S6_L001_R2_001.fastq.gz  -p 8 --validateMappings --gcBias -o quants/B50_1
salmon quant -i celegans_index -l A -1 $WORK/290921_JV_B50L2_S13_L001_R1_001.fastq.gz -2 $WORK/290921_JV_B50L2_S13_L001_R2_001.fastq.gz -p 8 --validateMappings --gcBias -o quants/B50_2
salmon quant -i celegans_index -l A -1 $WORK/290921_JV_B50L3_S20_L001_R1_001.fastq.gz -2 $WORK/290921_JV_B50L3_S20_L001_R2_001.fastq.gz -p 8 --validateMappings --gcBias -o quants/B50_3
salmon quant -i celegans_index -l A -1 $WORK/290921_JV_B50L4_S27_L001_R1_001.fastq.gz -2 $WORK/290921_JV_B50L4_S27_L001_R2_001.fastq.gz -p 8 --validateMappings --gcBias -o quants/B50_4

# gacA mutant control
salmon quant -i celegans_index -l A -1 $WORK/290921_JV_G0L1_S7_L001_R1_001.fastq.gz  -2 $WORK/290921_JV_G0L1_S7_L001_R2_001.fastq.gz  -p 8 --validateMappings --gcBias -o quants/G0_1
salmon quant -i celegans_index -l A -1 $WORK/290921_JV_G0L2_S14_L001_R1_001.fastq.gz -2 $WORK/290921_JV_G0L2_S14_L001_R2_001.fastq.gz -p 8 --validateMappings --gcBias -o quants/G0_2
salmon quant -i celegans_index -l A -1 $WORK/290921_JV_G0L3_S21_L001_R1_001.fastq.gz -2 $WORK/290921_JV_G0L3_S21_L001_R3_001.fastq.gz -p 8 --validateMappings --gcBias -o quants/G0_3
salmon quant -i celegans_index -l A -1 $WORK/290921_JV_G0L4_S28_L001_R1_001.fastq.gz -2 $WORK/290921_JV_G0L4_S28_L001_R4_001.fastq.gz -p 8 --validateMappings --gcBias -o quants/G0_4






