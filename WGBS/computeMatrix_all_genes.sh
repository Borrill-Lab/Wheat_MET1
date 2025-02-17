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
-S Z0001_all_chrom.bw Z0003_all_chrom.bw Z0007_all_chrom.bw Z0009_all_chrom.bw \
-R /jic/scratch/groups/Philippa-Borrill/References/iwgsc_ref_seq_1.1/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706_v2.gtf \
-o /jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/filter_by_chrom/average_methylation/gene_matrix_all.gz \
-m 2000 \
-b 2000 \
-a 2000

computeMatrix scale-regions \
-S Z0001_all_chrom_CHG.bw Z0003_all_chrom_CHG.bw Z0007_all_chrom_CHG.bw Z0009_all_chrom_CHG.bw \
-R /jic/scratch/groups/Philippa-Borrill/References/iwgsc_ref_seq_1.1/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706_v2.gtf \
-o /jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/filter_by_chrom/average_methylation/gene_matrix_all_CHG.gz \
-m 2000 \
-b 2000 \
-a 2000

computeMatrix scale-regions \
-S Z0001_all_chrom_CHH.bw Z0003_all_chrom_CHH.bw Z0007_all_chrom_CHH.bw Z0009_all_chrom_CHH.bw \
-R /jic/scratch/groups/Philippa-Borrill/References/iwgsc_ref_seq_1.1/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706_v2.gtf \
-o /jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/filter_by_chrom/average_methylation/gene_matrix_all_CHH.gz \
-m 2000 \
-b 2000 \
-a 2000
