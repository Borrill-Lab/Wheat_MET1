#!/bin/bash
#SBATCH -p jic-medium
#SBATCH -t 1-00:00
#SBATCH -c 4
#SBATCH --mem=200000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=delfi.dorussen@jic.ac.uk
#SBATCH -o /jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/filter_by_chrom/average_methylation/slurm_output/computematrix-slurm-%j.out
#SBATCH -e /jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/filter_by_chrom/average_methylation/slurm_output/computematrix-slurm-%j.err

source package 6daf0c37-1c5e-4cd6-9884-2d0b4d5f9d8f

cd /jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/filter_by_chrom/bigwig_files/

computeMatrix scale-regions \
-S Z0001_Chr5D_CG.bw Z0003_Chr5D_CG.bw Z0007_Chr5D_CG.bw Z0009_Chr5D_CG.bw \
-R /jic/scratch/groups/Philippa-Borrill/References/iwgsc_ref_seq_1.1/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706_v2_chr5D.gtf \
-o /jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/filter_by_chrom/average_methylation/gene_matrix_5D.gz \
-m 2000 \
-b 2000 \
-a 2000