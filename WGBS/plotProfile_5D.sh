#!/bin/bash
#SBATCH -p jic-medium
#SBATCH -t 1-00:00
#SBATCH -c 4
#SBATCH --mem=200000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=delfi.dorussen@jic.ac.uk
#SBATCH -o /jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/filter_by_chrom/average_methylation/slurm_output/plotprofile-slurm-%j.out
#SBATCH -e /jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/filter_by_chrom/average_methylation/slurm_output/plotprofile-slurm-%j.err

source package 6daf0c37-1c5e-4cd6-9884-2d0b4d5f9d8f

cd /jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/filter_by_chrom/average_methylation/

plotProfile -m gene_matrix_5A.gz \ #Edit input matrix file here eg for different methylation context
-o gene_CGmethylation_5A.pdf \
--averageType "mean" \
--plotType "lines" \
--colors black yellow blue red \
--perGroup \
--plotHeight 10 \
--plotWidth 15 \
--yMin 0 \
--yMax 1
