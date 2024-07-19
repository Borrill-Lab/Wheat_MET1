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
-a /jic/scratch/groups/Philippa-Borrill/References/iwgsc_ref_seq_1.1/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706_genes1kb.gff3 \
-b A_single_hypo_DMRs.bed #Change to sample of interest