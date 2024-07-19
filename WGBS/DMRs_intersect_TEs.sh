#!/bin/bash
#SBATCH -p jic-medium
#SBATCH -t 1-00:00
#SBATCH -c 4
#SBATCH --mem=200000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=delfi.dorussen@jic.ac.uk
#SBATCH -o /jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/DMR_overlap/slurm_output/overlap-slurm-%j.out
#SBATCH -e /jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/DMR_overlap/slurm_output/overlap-slurm-%j.err

source package /nbi/software/production/bin/bedtools-2.17.0

cd /jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/DMR_calling/

bedtools intersect -wa -wb \
-a /jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/DMR_overlap/iwgsc_refseqv1.0_TransposableElements_2017Mar13_complete.gff3 \
-b A_single_hyper_DMRs.bed #Change to sample of interest