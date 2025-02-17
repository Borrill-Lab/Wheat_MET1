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

plotProfile -m gene_matrix_all.gz \
-o gene_CGmethylation_all.pdf \
--averageType "mean" \
--plotType "lines" \
--colors black yellow blue red \
--perGroup \
--plotHeight 10 \
--plotWidth 15 \
--yMin 0 \
--yMax 1

plotProfile -m gene_matrix_all_CHG.gz \
-o gene_CHGmethylation_all.pdf \
--averageType "mean" \
--plotType "lines" \
--colors black yellow blue red \
--perGroup \
--plotHeight 10 \
--plotWidth 15 \
--yMin 0 \
--yMax 0.8

plotProfile -m gene_matrix_all_CHH.gz \
-o gene_CHHmethylation_all.pdf \
--averageType "mean" \
--plotType "lines" \
--colors black yellow blue red \
--perGroup \
--plotHeight 10 \
--plotWidth 15 \
--yMin 0 \
--yMax 0.4

plotProfile -m TE_matrix_all.gz \
-o TE_CGmethylation_all.pdf \
--averageType "mean" \
--plotType "lines" \
--colors black yellow blue red \
--perGroup \
--plotHeight 10 \
--plotWidth 15 \
--yMin 0 \
--yMax 1

plotProfile -m TE_matrix_all_CHG.gz \
-o TE_CHGmethylation_all.pdf \
--averageType "mean" \
--plotType "lines" \
--colors black yellow blue red \
--perGroup \
--plotHeight 10 \
--plotWidth 15 \
--yMin 0 \
--yMax 0.8

plotProfile -m TE_matrix_all_CHH.gz \
-o TE_CHHmethylation_all.pdf \
--averageType "mean" \
--plotType "lines" \
--colors black yellow blue red \
--perGroup \
--plotHeight 10 \
--plotWidth 15 \
--yMin 0 \
--yMax 0.4